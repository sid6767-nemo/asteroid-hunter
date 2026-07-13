# exp_faint_recovery.py - can the "hunt by ear" score recover KNOWN asteroids
# that the detection pipeline is too shallow to see?
#
# THE EXPERIMENT
#   The pipeline finds a dot in a SINGLE frame first, then links moving dots,
#   so it can only find objects bright enough to show up in one exposure.
#   The ear-tool pools the light from ALL frames along a track before deciding,
#   which reaches fainter objects (the whole point of shift-and-stack). Test:
#     1. Ask SkyBoT for EVERY known asteroid in the set203 field, per frame.
#     2. Turn each one's sky position into a pixel track via the frame-0 WCS.
#     3. Run the exact M0 score at each: TUNE (its real velocity) and ROAM
#        (blind - the tool guesses the velocity, like a user would).
#     4. Compare against random empty sky (the "fake floor") and the 3
#        pipeline-found asteroids (the "real" reference).
#   A faint KNOWN asteroid the pipeline MISSED but the ear scores high on = a
#   real recovery, not a false positive.
#
# Imports scripts/exp_hunt_score.py; does NOT modify it. The SkyBoT call needs
# internet (IMCCE), so run on your laptop:
#   .venv/bin/python scripts/exp_faint_recovery.py

import os
import sys
import glob

import numpy as np

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import exp_hunt_score as hs          # __main__ guard => import won't run its main()

from astropy.io import fits
from astropy.time import Time
from astropy.coordinates import SkyCoord
import astropy.units as u

PIPELINE_FOUND = {'2015 RM287', '2004 RH62', '2002 GE56'}
OBS_CODE = 'F51'                     # Pan-STARRS 1 (Haleakala); set203 is PS1 data
FIELD_RADIUS_DEG = 0.16
SEARCH_BOX = 9                       # cursor nudge radius, binned px
SEARCH_STEP = 3


def build_wcs(header):
    from astropy.wcs import WCS
    w = WCS(naxis=2)
    w.wcs.crpix = [header['CRPIX1'], header['CRPIX2']]
    w.wcs.cdelt = [header['CDELT1'], header['CDELT2']]
    w.wcs.crval = [header['CRVAL1'], header['CRVAL2']]
    w.wcs.ctype = [header['CTYPE1'], header['CTYPE2']]
    crota = np.radians(header['CROTA2'])
    w.wcs.pc = [[np.cos(crota), -np.sin(crota)],
                [np.sin(crota),  np.cos(crota)]]
    return w


def frame_headers_and_times():
    files = sorted(glob.glob(os.path.join(hs.DATA_DIR, '*.fits')),
                   key=lambda f: fits.open(f)[0].header['MJD-OBS'])
    headers = [fits.open(f)[0].header for f in files]
    mjds = [h['MJD-OBS'] for h in headers]
    return headers[0], np.array(mjds)


def skybot_at(center, epoch_mjd):
    """Every known solar-system object in the field at one epoch. Mirrors the
    pipeline's working call, including the units handling that keeps getting
    reverted: build the SkyCoord from the WHOLE RA/DEC column."""
    from astroquery.imcce import Skybot
    epoch = Time(epoch_mjd, format='mjd')
    try:
        tbl = Skybot.cone_search(center, FIELD_RADIUS_DEG * u.deg, epoch, location=OBS_CODE)
    except TypeError:
        tbl = Skybot.cone_search(center, FIELD_RADIUS_DEG * u.deg, epoch)
    if tbl is None or len(tbl) == 0:
        return []
    try:
        coords = SkyCoord(ra=tbl['RA'], dec=tbl['DEC'])
    except Exception:
        coords = SkyCoord(ra=tbl['_raj2000'], dec=tbl['_decj2000'], unit='deg')
    out = []
    for i in range(len(tbl)):
        try:
            vmag = float(tbl['V'][i])
        except Exception:
            vmag = np.nan
        out.append({'name': str(tbl['Name'][i]).strip(),
                    'ra': float(coords[i].ra.deg), 'dec': float(coords[i].dec.deg),
                    'V': vmag})
    return out


def collect_objects(center, mjds, w, binf):
    """Query every epoch, match by name, build each object's per-frame BINNED
    pixel track (linear fit over whatever frames SkyBoT reports it in)."""
    per_epoch = [skybot_at(center, m) for m in mjds]
    names = {o['name'] for ep in per_epoch for o in ep}
    t_days = mjds - mjds[0]
    objs = {}
    for name in names:
        ts, xs, ys, vmag = [], [], [], np.nan
        for k, ep in enumerate(per_epoch):
            for o in ep:
                if o['name'] != name:
                    continue
                px, py = w.all_world2pix(o['ra'], o['dec'], 0)
                ts.append(t_days[k]); xs.append(float(px) / binf); ys.append(float(py) / binf)
                if not np.isnan(o['V']):
                    vmag = o['V']
        if len(ts) < 2:
            continue
        ts = np.array(ts)
        ax, bx = np.polyfit(ts, xs, 1)
        ay, by = np.polyfit(ts, ys, 1)
        objs[name] = {'V': vmag, 'vx': ax, 'vy': ay, 'x0': bx, 'y0': by,
                      'n_epochs': len(ts)}
    return objs


def best_scores(frames, sigmas, t_days, x0, y0, vx, vy, box=SEARCH_BOX, step=SEARCH_STEP,
                roam_only=False):
    """Score on a grid of cursor nudges around the predicted spot, to absorb
    ephemeris/parallax error (fast objects predict least accurately). Returns
    best tune, best roam, the roam-recovered (speed, angle), and how far the
    winning cursor sat from the prediction (binned px). Pass roam_only=True for
    empty-sky floor sites so blank sky gets the SAME best-of-grid search the
    real objects get - otherwise the comparison is unfair."""
    offs = list(range(-box, box + 1, step))
    best_tune, best_roam, best_v, best_off = 0.0, 0.0, (0.0, 0.0), 0.0
    for dx in offs:
        for dy in offs:
            x, y = x0 + dx, y0 + dy
            if not roam_only:
                t = hs.score(frames, sigmas, t_days, x, y, vx, vy)
                if t > best_tune:
                    best_tune = t
            r, v = hs.roam_verified(frames, sigmas, t_days, x, y, detail=True)
            if r > best_roam:
                best_roam, best_v, best_off = r, v, float(np.hypot(dx, dy))
    return best_tune, best_roam, best_v, best_off


def main():
    print('Loading + aligning set203 ...')
    aligned, t_days = hs.load_and_align()
    frames, sigmas, _ = hs.prepare(aligned)
    H2, W2 = frames[0].shape
    Horig, Worig = aligned[0].shape
    binf = hs.BIN

    hdr0, mjds = frame_headers_and_times()
    w = build_wcs(hdr0)
    center = SkyCoord(*[float(v) for v in w.all_pix2world(Worig / 2, Horig / 2, 0)], unit='deg')
    print(f'  binned {W2}x{H2}, arc {t_days[-1]*24*60:.1f} min, '
          f'field center RA={center.ra.deg:.4f} Dec={center.dec.deg:.4f}')

    print('\nQuerying SkyBoT for every known asteroid in the field (4 epochs)...')
    objs = collect_objects(center, mjds, w, binf)
    print(f'  SkyBoT reports {len(objs)} known object(s) in/near the field.')

    MARGIN = 20
    rows = []
    for name, o in objs.items():
        x0, y0 = o['x0'], o['y0']
        if not (-MARGIN <= x0 < W2 + MARGIN and -MARGIN <= y0 < H2 + MARGIN):
            rows.append({'name': name, 'V': o['V'], 'note': 'off frame'})
            continue
        true_sp = np.hypot(o['vx'], o['vy']) * hs.PIXEL_SCALE * binf
        true_ang = np.degrees(np.arctan2(o['vy'], o['vx'])) % 360
        tune, roam, (rsp, rang), off = best_scores(frames, sigmas, t_days, x0, y0, o['vx'], o['vy'])
        rows.append({'name': name, 'V': o['V'], 'true_sp': true_sp, 'true_ang': true_ang,
                     'tune': tune, 'roam': roam, 'rec_sp': rsp, 'rec_ang': rang, 'off': off,
                     'found': name in PIPELINE_FOUND, 'note': ''})

    print('Scoring random empty sky (same best-of-grid search) for a FAIR floor ...')
    rng = np.random.default_rng(0)
    m = 60
    fake = np.array([best_scores(frames, sigmas, t_days,
                                 rng.uniform(m, W2 - m), rng.uniform(m, H2 - m),
                                 0.0, 0.0, roam_only=True)[1]
                     for _ in range(12)])
    floor = float(fake.max())

    def verdict(r):
        if r['note']:
            return r['note']
        if r['found']:
            return 'pipeline control'
        if r['roam'] >= 3 * floor:
            return 'RECOVERED (well above floor)'
        if r['roam'] >= 2 * floor:
            return 'candidate (above floor)'
        return 'not above floor'

    rows.sort(key=lambda r: (r['V'] if r['V'] is not None and not np.isnan(r['V']) else 99))
    print('\n' + '=' * 96)
    print(f'{"asteroid":13s} {"V":>5s} {"speed":>7s} {"tune":>7s} {"roam":>7s} '
          f'{"off":>5s} {"pipeline":>8s}   verdict')
    print('-' * 96)
    for r in rows:
        if r.get('note'):
            print(f'{r["name"]:13s} {r["V"]:5.1f} {"":>7s} {"":>7s} {"":>7s} {"":>5s} '
                  f'{"":>8s}   {r["note"]}')
            continue
        pf = 'FOUND' if r['found'] else 'missed'
        print(f'{r["name"]:13s} {r["V"]:5.1f} {r["true_sp"]:7.0f} {r["tune"]:7.3f} '
              f'{r["roam"]:7.3f} {r["off"]:5.0f} {pf:>8s}   {verdict(r)}')
    print('-' * 96)
    print(f'FAIR fake floor (best-of-grid roam, 12 empty sites) = {floor:.3f}   '
          f'median = {np.median(fake):.3f}')
    print('=' * 96)

    # velocity-match diagnostic: does the ear's recovered track match the truth?
    print('\nTRACK MATCH (does the ear\'s recovered velocity match SkyBoT truth?):')
    for r in rows:
        if r.get('note'):
            continue
        dang = abs((r['rec_ang'] - r['true_ang'] + 180) % 360 - 180)
        ratio = r['rec_sp'] / r['true_sp'] if r['true_sp'] else 0
        ok = 'MATCH' if (dang < 15 and 0.7 < ratio < 1.4) else 'no match'
        tag = 'ctrl  ' if r['found'] else 'missed'
        print(f'  [{tag}] {r["name"]:13s} recovered {r["rec_sp"]:4.0f}"/d @ {r["rec_ang"]:3.0f}deg  '
              f'vs true {r["true_sp"]:4.0f}"/d @ {r["true_ang"]:3.0f}deg   '
              f'(dAngle {dang:.0f}deg, speed x{ratio:.2f})  -> {ok}')
    print('\nA "missed" object that is above the floor AND whose recovered track MATCHES '
          'the\ntrue motion is a genuine faint recovery - the ear seeing what the pipeline '
          'cannot.')


if __name__ == '__main__':
    main()
