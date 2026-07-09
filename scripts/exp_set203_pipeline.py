# exp_set203_pipeline.py - moving-object detection for ANY number of frames (3+)
#
# Run on the sample data:        python scripts/exp_set203_pipeline.py
# Run on your own folder:        python scripts/exp_set203_pipeline.py --data data/myfield
# See all options:               python scripts/exp_set203_pipeline.py --help

import os
import glob
import argparse
import warnings
import numpy as np
from astropy.io import fits
import astroalign as aa
from astropy.stats import sigma_clipped_stats
from scipy.spatial import KDTree
from scipy import ndimage
from photutils.segmentation import detect_sources, deblend_sources, SourceCatalog
import matplotlib.pyplot as plt

warnings.filterwarnings('ignore', message='.*deblending mode.*')

PIXEL_SCALE      = 0.256
MIN_MOVE         = 3
MAX_MOVE         = 100
CONFIRM_RADIUS   = 5
GIANT_AREA_MIN   = 150
GIANT_RADIUS_MIN = 180
CV_MAX           = 0.35
ELONG_MAX        = 2.3    # max stretched-ness; rejects diffraction-spike artifacts
CONC_MIN         = 0.20   # min concentration for FAINT candidates
CONC_BRIGHT      = 0.32   # stricter concentration required for BRIGHT candidates
SNR_BRIGHT       = 15.0   # above this SNR a candidate counts as 'bright'
MAX_FRAMES       = 10   # safety cap


def get_args():
    p = argparse.ArgumentParser(
        description="Detect moving asteroids across a set of telescope frames (3 or more).")
    p.add_argument('--data', default='data/set203',
                   help="folder containing the .fits frames (default: data/set203)")
    p.add_argument('--output', default='outputs',
                   help="folder to save result images (default: outputs)")
    p.add_argument('--threshold', type=float, default=3.0,
                   help="detection threshold in sigma above background (default: 3.0)")
    p.add_argument('--giant-radius-k', type=float, default=10.0,
                   help="size of rejection zones around bright stars (default: 10.0)")
    p.add_argument('--min-frames', type=int, default=None,
                   help="how many frames an object must appear in "
                        "(default: all of them)")
    p.add_argument('--rms-max', type=float, default=1.0,
                   help="max straight-line fit error in px (default: 1.0)")
    p.add_argument('--no-skybot', action='store_true',
                   help="skip the SkyBoT online cross-match")
    p.add_argument('--save-frames', action='store_true',
                   help="also save each aligned frame as its own image (for the web blink viewer)")
    return p.parse_args()


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


def classify(rate):
    if rate < 1:    return "too slow -> likely star or artifact"
    elif rate < 50: return "slow -> TNO-like"
    elif rate < 500:return "main-belt asteroid range"
    else:           return "fast -> NEO-like"


def main():
    args = get_args()
    DATA_DIR, OUTPUT_DIR = args.data, args.output
    THRESH_SIGMA, GIANT_RADIUS_K = args.threshold, args.giant_radius_k
    RMS_MAX, DO_SKYBOT = args.rms_max, not args.no_skybot
    SAVE_FRAMES = args.save_frames

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    dataset_name = os.path.basename(os.path.normpath(DATA_DIR))

    # --- load ALL frames, sorted by time ---
    files = sorted(glob.glob(os.path.join(DATA_DIR, '*.fits')),
                   key=lambda f: fits.open(f)[0].header['MJD-OBS'])
    N = len(files)
    if N < 3:
        print(f"ERROR: found {N} .fits files in '{DATA_DIR}', need at least 3.")
        return
    if N > MAX_FRAMES:
        print(f"Note: found {N} frames, using the first {MAX_FRAMES} by time.")
        files = files[:MAX_FRAMES]; N = MAX_FRAMES

    MIN_FRAMES = args.min_frames if args.min_frames is not None else N
    MIN_FRAMES = max(3, min(MIN_FRAMES, N))
    print(f"Loading {N} frames from '{DATA_DIR}' (chronological order)...")

    images, mjds, headers = [], [], []
    for f in files:
        with fits.open(f) as hdul:
            images.append(hdul[0].data.astype(float))
            mjds.append(hdul[0].header['MJD-OBS'])
            headers.append(hdul[0].header)
            print(f"  {os.path.basename(f)[:24]}...  MJD={mjds[-1]:.6f}")

    t_days = [mjds[i] - mjds[0] for i in range(N)]
    gaps_min = [(mjds[i+1]-mjds[i])*24*60 for i in range(N-1)]
    print(f"\n{N} frames, time gaps (min): " + ", ".join(f"{g:.1f}" for g in gaps_min))
    print(f"Requiring an object to appear in at least {MIN_FRAMES} of {N} frames.")

    # --- align all frames to frame 0 ---
    print("\nAligning frames to frame 1...")
    aligned = [images[0]]
    for i in range(1, N):
        try:
            reg, _ = aa.register(images[i], images[0], detection_sigma=5)
            aligned.append(reg); print(f"  Frame {i+1} aligned OK")
        except Exception as e:
            print(f"  Frame {i+1} FAILED: {e}"); aligned.append(images[i])

    # --- background / threshold / saturation (frame 1) ---
    mask0 = (aligned[0] == 0)
    _, ref_median, ref_std = sigma_clipped_stats(aligned[0], sigma=3.0, mask=mask0)
    THRESHOLD = THRESH_SIGMA * ref_std
    img_median = np.median(aligned[0][aligned[0] > 0])
    SATURATION_LEVEL = img_median + 200 * ref_std

    # --- giant-star rejection zones ---
    lbl, nblob = ndimage.label(aligned[0] > SATURATION_LEVEL)
    sizes = ndimage.sum(np.ones_like(lbl), lbl, range(1, nblob + 1))
    giant_zones = []
    for i in range(1, nblob + 1):
        if sizes[i-1] >= GIANT_AREA_MIN:
            ys, xs = np.where(lbl == i)
            giant_zones.append((xs.mean(), ys.mean(),
                                max(GIANT_RADIUS_MIN, GIANT_RADIUS_K*np.sqrt(sizes[i-1]))))
    def in_giant_zone(x, y):
        return any(np.hypot(x-cx, y-cy) < r for cx, cy, r in giant_zones)

    # --- detect sources per frame ---
    print("\nDetecting sources (segmentation + deblending)...")
    all_pos, all_peak, all_elong = [], [], []
    for i, img in enumerate(aligned):
        mask = (img == 0)
        _, median, _ = sigma_clipped_stats(img, sigma=3.0, mask=mask)
        data = img - median
        segm = detect_sources(data, THRESHOLD, n_pixels=5, mask=mask)
        if segm is None:
            all_pos.append(np.empty((0,2))); all_peak.append(np.empty(0)); all_elong.append(np.empty(0))
            print(f"  Frame {i+1}: 0 sources"); continue
        segm = deblend_sources(data, segm, n_pixels=5, n_levels=32, contrast=0.0001)
        cat = SourceCatalog(data, segm, mask=mask)
        all_pos.append(np.array(list(zip(cat.x_centroid, cat.y_centroid))))
        all_peak.append(np.array(cat.max_value, dtype=float))
        all_elong.append(np.array(cat.elongation, dtype=float))
        print(f"  Frame {i+1}: {len(all_pos[-1])} sources")

    trees = [KDTree(p) if len(p) else None for p in all_pos]

    def mutual_pairs(iA, iB):
        A, B = all_pos[iA], all_pos[iB]
        if len(A) == 0 or len(B) == 0: return []
        tA, tB = KDTree(A), KDTree(B)
        df, jf = tB.query(A); _, jb = tA.query(B)
        return [(A[i], B[j], d) for i,(d,j) in enumerate(zip(df,jf))
                if MIN_MOVE < d < MAX_MOVE and jb[j] == i]

    confirmed = []

    def track_quality(fpos, fpk):
        ts = np.array(t_days)
        xs = np.array([p[0] for p in fpos]); ys = np.array([p[1] for p in fpos])
        px = np.polyfit(ts, xs, 1); py = np.polyfit(ts, ys, 1)
        rms = np.sqrt(np.mean((xs-np.polyval(px,ts))**2 + (ys-np.polyval(py,ts))**2))
        speed = np.hypot(px[0], py[0]) * PIXEL_SCALE
        pk = np.array(fpk); pk = pk[~np.isnan(pk)]
        cv = (pk.std()/pk.mean()) if len(pk) >= 2 and pk.mean() > 0 else 9.9
        return rms, cv, speed

    def stack_concentration(fpos):
        # Stack cutouts centered on the tracked position in each frame; measure
        # how concentrated the light is at the center. Real objects pile up into
        # a central point (high); halo/galaxy/noise fakes spread out (low).
        Z = 16; cuts = []
        for fi in range(N):
            x, y = int(round(fpos[fi][0])), int(round(fpos[fi][1]))
            if y-Z < 0 or x-Z < 0 or y+Z > aligned[fi].shape[0] or x+Z > aligned[fi].shape[1]:
                continue
            c = aligned[fi][y-Z:y+Z, x-Z:x+Z]
            _, md, _ = sigma_clipped_stats(c)
            cuts.append(c - md)
        if not cuts:
            return 0.0
        st = np.clip(np.mean(cuts, axis=0), 0, None)
        yy, xx = np.mgrid[0:2*Z, 0:2*Z]; r = np.hypot(xx - Z, yy - Z)
        outer = st[r < 12].sum()
        return float(st[r < 3].sum() / outer) if outer > 0 else 0.0

    def line_check_and_add(pairs, start_idx):
        # start_idx = index of the first frame of the seed pair (step is per-frame)
        for A, B, dist in pairs:
            step = B - A
            frame_pred = [A + step*(k - start_idx) for k in range(N)]
            fpos, fpk, fel, hits = [], [], [], 0
            for fi, pred in enumerate(frame_pred):
                if trees[fi] is None:
                    fpos.append(pred); fpk.append(np.nan); fel.append(np.nan); continue
                d, idx = trees[fi].query(pred)
                if d < CONFIRM_RADIUS:
                    hits += 1; fpos.append(all_pos[fi][idx]); fpk.append(all_peak[fi][idx]); fel.append(all_elong[fi][idx])
                else:
                    fpos.append(pred); fpk.append(np.nan); fel.append(np.nan)
            if hits < MIN_FRAMES:
                continue
            f1 = fpos[0]
            if any(in_giant_zone(p[0], p[1]) for p in fpos):
                continue
            rms, cv, speed = track_quality(fpos, fpk)
            if rms > RMS_MAX or cv > CV_MAX:
                continue
            mean_elong = float(np.nanmean(fel))
            if mean_elong > ELONG_MAX:      # reject stretched-out spike artifacts
                continue
            conc = stack_concentration(fpos)
            mean_snr = float(np.nanmean(fpk)) / ref_std
            conc_needed = CONC_BRIGHT if mean_snr > SNR_BRIGHT else CONC_MIN
            if conc < conc_needed:          # bright objects must be tightly concentrated; faint get a pass
                continue
            if any(np.hypot(f1[0]-e['f1'][0], f1[1]-e['f1'][1]) < 30 for e in confirmed):
                continue
            confirmed.append({'f1': f1, 'fpos': fpos, 'fpk': fpk, 'hits': hits,
                              'rate': speed, 'rms': rms, 'cv': cv, 'conc': conc})
            print(f"  CONFIRMED #{len(confirmed)}: {speed:.0f}\"/day  "
                  f"({f1[0]:.0f},{f1[1]:.0f})  {hits}/{N}  RMS={rms:.2f} CV={cv:.2f} conc={conc:.2f}")

    print("\nLinking + filtering tracks...")
    for i in range(N - 1):
        line_check_and_add(mutual_pairs(i, i+1), i)
    print(f"\n{len(confirmed)} confirmed candidate(s) after filtering")

    print("\n=== FINAL CONFIRMED CANDIDATES ===")
    for n, c in enumerate(confirmed):
        print(f"\nCandidate #{n+1}")
        print(f"  Speed:     {c['rate']:.0f} arcsec/day  ->  {classify(c['rate'])}")
        print(f"  Frame-1:   ({c['f1'][0]:.0f},{c['f1'][1]:.0f})")
        print(f"  Quality:   {c['hits']}/{N} frames, linRMS={c['rms']:.2f}px, brightnessCV={c['cv']:.2f}")

    # --- Stage 1: write sky-coordinate astrometry (time, RA, Dec) per detection ---
    if confirmed:
        from astropy.wcs import WCS as _WCS
        from astropy.time import Time as _Time
        from astropy.coordinates import SkyCoord as _SkyCoord
        import astropy.units as _u
        wq = build_wcs(headers[0])   # frame-0 WCS maps the aligned grid -> sky

        def _mpc_line(desig, mjd, ra_deg, dec_deg, obscode='F51'):
            t = _Time(mjd, format='mjd')
            y, mo, d = t.datetime.year, t.datetime.month, t.datetime.day
            dayfrac = (t.mjd - _Time(f'{y}-{mo:02d}-{d:02d}', format='iso').mjd) + d
            c = _SkyCoord(ra_deg*_u.deg, dec_deg*_u.deg)
            ra_hms = c.ra.to_string(unit=_u.hour, sep=' ', precision=2, pad=True)
            dec_dms = c.dec.to_string(unit=_u.deg, sep=' ', precision=1, alwayssign=True, pad=True)
            return f"     {desig:<7} C{y} {mo:02d} {dayfrac:08.5f} {ra_hms} {dec_dms}                   {obscode}"

        astro_path = os.path.join(OUTPUT_DIR, f'{dataset_name}_astrometry.txt')
        with open(astro_path, 'w') as af:
            af.write(f"# Astrometry for {dataset_name}: sky position of each detection at each frame time\n")
            af.write(f"# Feed the MPC-format lines to an orbit tool (find_orb / OpenOrb) to fit an orbit.\n\n")
            for n, c in enumerate(confirmed):
                desig = f"CAND{n+1:03d}"
                af.write(f"# Candidate #{n+1}  ({c['rate']:.0f} arcsec/day)\n")
                af.write(f"#  frame  {'UTC time':23}  {'RA(deg)':>10}  {'Dec(deg)':>10}\n")
                mpc_lines = []
                for fi, pos in enumerate(c['fpos']):
                    ra, dec = wq.all_pix2world(pos[0], pos[1], 0)
                    ra, dec = float(ra), float(dec)
                    utc = _Time(mjds[fi], format='mjd').iso
                    af.write(f"#  {fi+1:<5}  {utc:23}  {ra:10.5f}  {dec:10.5f}\n")
                    mpc_lines.append(_mpc_line(desig, mjds[fi], ra, dec))
                af.write("# MPC 80-column observations:\n")
                for ln in mpc_lines:
                    af.write(ln + "\n")
                af.write("\n")
        print(f"Saved astrometry -> {astro_path}")

    # --- Sonification: sound-field (from real tracks) + backdrop image for the web "Listen" mode ---
    if confirmed:
        try:
            from scipy.spatial import KDTree as _KDTree
            from scipy.ndimage import gaussian_filter as _gf
            from PIL import Image as _Image
            H0, W0 = aligned[0].shape
            NW = 700; _scale = NW / W0; NH = int(H0 * _scale)

            # PER-FRAME sound-fields: one beacon per asteroid AT ITS POSITION IN THAT FRAME.
            # The user scans each frame by ear and locks where they think the object is;
            # because the beacon MOVES frame to frame, they discover the motion themselves.
            import json as _json2
            ys, xs = np.mgrid[0:NH, 0:NW].astype(float)
            R = 90.0                                    # reach in display px
            nfr = len(aligned)
            for fi in range(nfr):
                field = np.zeros((NH, NW))
                for c in confirmed:
                    cx, cy = np.array(c['fpos'])[fi] * _scale
                    field = np.maximum(field, np.exp(-((xs-cx)**2 + (ys-cy)**2) / (2*(R/2.5)**2)))
                _Image.fromarray((field * 255).astype(np.uint8)).save(
                    os.path.join(OUTPUT_DIR, f'{dataset_name}_soundfield_f{fi+1}.png'))
            # also keep a combined field (frame-1 beacons) for the simple single-frame listen mode
            field1 = np.zeros((NH, NW))
            for c in confirmed:
                cx, cy = np.array(c['fpos'])[0] * _scale
                field1 = np.maximum(field1, np.exp(-((xs-cx)**2 + (ys-cy)**2) / (2*(R/2.5)**2)))
            _Image.fromarray((field1 * 255).astype(np.uint8)).save(
                os.path.join(OUTPUT_DIR, f'{dataset_name}_soundfield.png'))
            # truth positions (display px) per asteroid per frame, for scoring the user's locks
            truth = {'nframes': nfr, 'asteroids': [
                {'id': n+1, 'positions': [[round(float(p[0]*_scale),1), round(float(p[1]*_scale),1)]
                                          for p in np.array(c['fpos'])]}
                for n, c in enumerate(confirmed)]}
            with open(os.path.join(OUTPUT_DIR, f'{dataset_name}_truth.json'), 'w') as _tf:
                _json2.dump(truth, _tf)

            # backdrop: stretched frame-1 so a sighted tester can also see the field
            ref0 = aligned[0]
            m0 = (ref0 == 0)
            from astropy.stats import sigma_clipped_stats as _scs
            _, med0, sd0 = _scs(ref0, sigma=3.0, mask=m0)
            disp = np.arcsinh(np.clip((ref0 - med0) / (sd0 * 4 if sd0 > 0 else 1), 0, None))
            disp = np.clip(disp / np.percentile(disp, 99.6), 0, 1)
            _Image.fromarray((disp * 255).astype(np.uint8)).resize((NW, NH)).save(
                os.path.join(OUTPUT_DIR, f'{dataset_name}_backdrop.jpg'), quality=72)
            print(f"Saved sonification -> {dataset_name}_soundfield.png + _backdrop.jpg")
        except Exception as _e:
            print(f"  (sonification skipped: {_e})")

    # --- visualization (one panel per frame) ---
    print("\nGenerating visualization...")
    fig, axes = plt.subplots(1, N, figsize=(5*N, 6))
    if N == 1: axes = [axes]
    fig.suptitle(f'Confirmed moving objects: {len(confirmed)}', fontsize=14)
    colors = ['cyan','orange','red','lime','magenta','yellow']
    for fi, (ax, img) in enumerate(zip(axes, aligned)):
        _, med, std = sigma_clipped_stats(img, sigma=3.0)
        ax.imshow(img, cmap='gray', vmin=med-2*std, vmax=med+4*std, origin='upper')
        ax.set_title(f'Frame {fi+1}')
        for n, c in enumerate(confirmed):
            pos = c['fpos'][fi]; col = colors[n % len(colors)]
            ax.add_patch(plt.Circle((pos[0], pos[1]), 15, color=col, fill=False, lw=2))
            ax.text(pos[0]+18, pos[1]+18, f'#{n+1}', color=col, fontsize=8, fontweight='bold')
    plt.tight_layout()
    out_png = os.path.join(OUTPUT_DIR, f'{dataset_name}_tracks.png')
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    print(f"Saved {out_png}")

    # individual aligned frames, all identical size, for the web blink viewer
    if SAVE_FRAMES:
        H, W = aligned[0].shape
        for fi, img in enumerate(aligned):
            f2 = plt.figure(figsize=(8, 8 * H / W))
            a2 = f2.add_axes([0, 0, 1, 1])   # axes fill the whole figure -> identical dims
            _, med, std = sigma_clipped_stats(img, sigma=3.0)
            a2.imshow(img, cmap='gray', vmin=med-2*std, vmax=med+4*std, origin='upper')
            a2.axis('off')
            for n, c in enumerate(confirmed):
                pos = c['fpos'][fi]; col = colors[n % len(colors)]
                a2.add_patch(plt.Circle((pos[0], pos[1]), 15, color=col, fill=False, lw=1.8))
                a2.text(pos[0]+20, pos[1]+20, f'#{n+1}', color=col, fontsize=11, fontweight='bold')
            a2.set_xlim(0, W); a2.set_ylim(H, 0)
            f2.savefig(os.path.join(OUTPUT_DIR, f'{dataset_name}_frame{fi+1}.png'), dpi=130)
            plt.close(f2)
        print(f"Saved {N} individual frames for the viewer")

    plt.show()

    # --- write a results file (candidate info + SkyBoT name if matched) for the web app ---
    import json as _json
    from astropy.coordinates import SkyCoord as _SC
    from astropy.time import Time as _T
    import astropy.units as _u2
    w = build_wcs(headers[0]); epoch = _T(mjds[0], format='mjd')
    results = []
    for n, c in enumerate(confirmed):
        ra, dec = w.all_pix2world(c['f1'][0], c['f1'][1], 0)
        _H0, _W0 = aligned[0].shape
        _frac = [[round(float(p[0])/_W0, 5), round(float(p[1])/_H0, 5)] for p in c['fpos']]
        results.append({'id': n+1, 'rate': round(float(c['rate'])),
                        'x': round(float(c['f1'][0])), 'y': round(float(c['f1'][1])),
                        'rms': round(float(c['rms']),2), 'cv': round(float(c['cv']),2),
                        'conc': round(float(c.get('conc',0)),2),
                        'ra': round(float(ra),5), 'dec': round(float(dec),5),
                        'fpos': _frac,                      # per-frame position, as fractions of the image
                        'name': None, 'sep_arcsec': None})

    if DO_SKYBOT:
        try:
            from astroquery.imcce import Skybot
        except Exception:
            Skybot = None
        if Skybot is not None:
            print("\nCross-matching candidates against SkyBoT...")
            for n, c in enumerate(confirmed):
                ra, dec = results[n]['ra'], results[n]['dec']
                print(f"\n--- Candidate #{n+1}: RA={ra:.5f} Dec={dec:.5f} ({c['rate']:.0f}\"/day) ---")
                try:
                    tbl = Skybot.cone_search(_SC(ra=ra, dec=dec, unit='deg'), 0.05*_u2.deg, epoch)
                    print(tbl)
                    if tbl is not None and len(tbl) > 0:
                        cand = _SC(ra=ra, dec=dec, unit='deg')
                        # SkyBoT's RA/DEC columns carry units; build the SkyCoord from the
                        # whole column at once rather than float()-ing each element.
                        try:
                            objs = _SC(ra=tbl['RA'], dec=tbl['DEC'])
                        except Exception:
                            objs = _SC(ra=tbl['_raj2000'], dec=tbl['_decj2000'])
                        seps = cand.separation(objs).arcsec
                        import numpy as _np2
                        best_i = int(_np2.argmin(seps)); best_sep = float(seps[best_i])
                        if best_sep < 15:
                            results[n]['name'] = str(tbl['Name'][best_i]).strip()
                            results[n]['sep_arcsec'] = round(best_sep,1)
                except Exception as e:
                    print(f"  No SkyBoT match / query failed: {e}")

    with open(os.path.join(OUTPUT_DIR, f'{dataset_name}_results.json'), 'w') as jf:
        _json.dump(results, jf)
    print(f"\nSaved results -> {os.path.join(OUTPUT_DIR, dataset_name + '_results.json')}")
    print("Done.")


if __name__ == '__main__':
    main()
