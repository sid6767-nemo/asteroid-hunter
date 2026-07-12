# exp_hunt_score.py - Milestone 0 feasibility test for "Hunt by Ear"
#
# Question this script answers BEFORE we build any web/audio code:
#   Can a shift-and-stack score, computed from raw pixels + a guessed velocity,
#   separate the 3 real set203 asteroids from stars and random sky?
#
# The score for a guess (x, y, speed, angle):
#   - sample a small patch in each frame at the position the object WOULD be
#     if it moved at the guessed velocity (this is shift-and-stack)
#   - c_k  = core flux above local background, in noise (sigma) units, frame k
#   - s    = min over frames of c_k   ("present in EVERY frame along the track")
#   - conc = fraction of the stacked light inside the tight core ("a point, not a smear")
#   - score = s * conc
#
# Pass criterion: worst asteroid score > 5x best star/random score.
#
# Run:  .venv/bin/python scripts/exp_hunt_score.py

import glob
import json
import os

import numpy as np
import astroalign as aa
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from scipy import ndimage

PIXEL_SCALE = 0.256          # arcsec per ORIGINAL pixel
BIN = 2                      # 2x2 binning, same as the future web export
SPEED_FLOOR = 70.0           # arcsec/day; slower is indistinguishable from a star in 63 min
CORE_R = 2.5                 # binned px, core aperture radius
ANN_IN, ANN_OUT = 6.0, 9.0   # binned px, local-background annulus
CONC_R = 7.0                 # binned px, outer radius for the concentration ratio
PATCH = 9                    # patch half-size: samples u,v in [-9, +9]

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DATA_DIR = os.path.join(ROOT, 'data', 'set203')
RESULTS_JSON = os.path.join(ROOT, 'outputs', 'set203_results.json')
WEB_RESULTS = os.path.join(ROOT, 'web', 'static', 'results')


def load_and_align():
    files = sorted(glob.glob(os.path.join(DATA_DIR, '*.fits')),
                   key=lambda f: fits.open(f)[0].header['MJD-OBS'])
    images, mjds = [], []
    for f in files:
        with fits.open(f) as hdul:
            images.append(hdul[0].data.astype(float))
            mjds.append(hdul[0].header['MJD-OBS'])
    t_days = np.array([m - mjds[0] for m in mjds])
    aligned = [images[0]]
    for i in range(1, len(images)):
        reg, _ = aa.register(images[i], images[0], detection_sigma=5)
        aligned.append(reg)
    return aligned, t_days


def bin2(img):
    h, w = img.shape
    h2, w2 = h // BIN * BIN, w // BIN * BIN
    return img[:h2, :w2].reshape(h2 // BIN, BIN, w2 // BIN, BIN).mean(axis=(1, 3))


def prepare(aligned):
    """Bin, background-subtract, then REMOVE THE STATIC SKY: the per-pixel
    median across the unshifted frames is made of things that sit still
    (stars), so subtracting it leaves only things that move (plus noise).
    An asteroid touches any given pixel in only 1 of 4 frames, so the median
    ignores it and its flux survives the subtraction. No-data pixels -> NaN."""
    frames, sat = [], None
    for img in aligned:
        b = bin2(img)
        nodata = bin2((img == 0).astype(float)) > 0
        _, med, sig = sigma_clipped_stats(b[~nodata], sigma=3.0)
        f = b - med
        f[nodata] = np.nan
        frames.append(f.astype(np.float32))
        # saturated pixels hold no valid measurement in ANY frame; mask them
        s = f > 200 * sig
        sat = s if sat is None else (sat | s)
    sat = ndimage.binary_dilation(sat, iterations=3)
    for f in frames:
        f[sat] = np.nan
    with np.errstate(invalid='ignore'):
        static_sky = np.nanmedian(frames, axis=0)
    diffs, sigmas = [], []
    for f in frames:
        d = f - static_sky
        _, _, sig = sigma_clipped_stats(d[~np.isnan(d)], sigma=3.0)
        diffs.append(d.astype(np.float32))
        sigmas.append(sig)
    return diffs, np.array(sigmas), frames[0]      # frames[0]: stars still present


# --- the score engine (this is what gets ported to JS) -----------------------

_uv = np.mgrid[-PATCH:PATCH + 1, -PATCH:PATCH + 1]          # (2, 19, 19)
_r = np.hypot(_uv[0], _uv[1])


def make_masks(core_r, ann_in, ann_out, conc_out):
    core = _r <= core_r
    return core, (_r > ann_in) & (_r <= ann_out), _r <= conc_out, int(core.sum())


SHARP = make_masks(CORE_R, ANN_IN, ANN_OUT, CONC_R)   # tuning: tight, precise
BROAD = make_masks(5.0, 7.0, 9.0, 9.0)                # roaming: wide, forgiving


def sample_patch(frame, x, y):
    """Bilinear-sample a (2*PATCH+1)^2 patch centered on (x, y) in binned px."""
    coords = np.array([_uv[0] + y, _uv[1] + x])              # row=y, col=x
    return ndimage.map_coordinates(frame, coords.reshape(2, -1), order=1,
                                   mode='constant', cval=np.nan).reshape(_r.shape)


def raw_s(frames, sigmas, t_days, x, y, vx, vy, masks):
    """min-over-frames core signal along a track, in LOCAL noise units,
    plus the stacked patches (for the concentration term)."""
    m_core, m_ann, m_out, n_core = masks
    cs, patches = [], []
    for k, f in enumerate(frames):
        p = sample_patch(f, x + vx * t_days[k], y + vy * t_days[k])
        if np.isnan(p[m_core]).any():
            continue                       # off-image or saturated: drop frame
        ann = p[m_ann]
        b = np.nanmedian(ann)
        # local noise: residual-infested neighborhoods are noisy; judging the
        # core against ITS OWN surroundings crushes them without touching
        # clean sky. (1.4826 * MAD = robust sigma estimate)
        loc = 1.4826 * np.nanmedian(np.abs(ann - b))
        noise = max(loc, 0.6 * sigmas[k])
        cs.append(np.sum(p[m_core] - b) / (noise * n_core))
        patches.append(p - b)
    if len(cs) < 3:
        return 0.0, None                                      # never invent data
    s = max(0.0, min(cs))
    if s > 0:
        # consistency: an asteroid's brightness is ~constant over one hour, so
        # its per-frame signals agree; junk tracks stitched from unrelated
        # residuals are wildly inconsistent (same physics as the pipeline's
        # CV_MAX filter). cv = scatter / mean of the per-frame signals.
        cv = float(np.std(cs) / np.mean(cs))
        s *= max(0.0, 1.0 - cv)
    return s, patches


def score(frames, sigmas, t_days, x, y, vx, vy, masks=SHARP, s_static=None):
    """(x, y) binned px in frame 0; (vx, vy) binned px per day.
    Zero-velocity veto: light piling up along the track only counts if it
    CANNOT be explained by something sitting still at (x, y) -- a static
    residual blob scores high at v=0 too, a real mover scores ~0 there."""
    if s_static is None:
        s_static, _ = raw_s(frames, sigmas, t_days, x, y, 0.0, 0.0, masks)
    s, patches = raw_s(frames, sigmas, t_days, x, y, vx, vy, masks)
    s = max(0.0, s - s_static)
    if s == 0.0 or patches is None:
        return 0.0
    m_core, _, m_out, _ = masks
    with np.errstate(invalid='ignore'):
        M = np.clip(np.nanmedian(patches, axis=0), 0, None)
    outer = np.nansum(M[m_out])
    conc = float(np.nansum(M[m_core]) / outer) if outer > 0 else 0.0
    return s * conc


_CORE_OFF = np.argwhere(_r <= CORE_R) - PATCH               # (21, 2) dy,dx offsets


_BANK_SP = np.arange(70.0, 1251.0, 25.0)
_BANK_AN = np.radians(np.arange(0.0, 360.0, 2.0))
_BANK_VX = ((_BANK_SP[:, None] / (PIXEL_SCALE * BIN)) * np.cos(_BANK_AN)[None, :]).ravel()
_BANK_VY = ((_BANK_SP[:, None] / (PIXEL_SCALE * BIN)) * np.sin(_BANK_AN)[None, :]).ravel()


def roam_bank(frames, sigmas, t_days, x, y):
    """Fast matched-filter bank: try EVERY velocity on a fine grid (25 "/day
    x 2 deg = 8640 tracks), sampling only the 21 sharp-core pixels per track
    (median subtraction already removed the background, so b=0 and global
    noise are good enough for this quick-look stage). Vectorized: one
    map_coordinates call per frame covers all 8640 tracks."""
    n = _BANK_VX.size
    cs = np.empty((len(frames), n))
    for k, f in enumerate(frames):
        yy = (y + _BANK_VY * t_days[k])[:, None] + _CORE_OFF[:, 0]
        xx = (x + _BANK_VX * t_days[k])[:, None] + _CORE_OFF[:, 1]
        vals = ndimage.map_coordinates(f, [yy.ravel(), xx.ravel()], order=1,
                                       mode='constant', cval=np.nan).reshape(n, -1)
        with np.errstate(invalid='ignore'):
            cs[k] = np.nanmean(vals, axis=1) / sigmas[k]
        cs[k][np.isnan(vals).any(axis=1)] = np.nan
    valid = (~np.isnan(cs)).sum(axis=0) >= 3
    with np.errstate(invalid='ignore', divide='ignore'):
        s = np.clip(np.nanmin(cs, axis=0), 0, None)
        consist = np.clip(1.0 - np.nanstd(cs, axis=0)
                          / np.maximum(np.nanmean(cs, axis=0), 1e-9), 0, 1)
    s = s * consist
    s[~valid] = 0.0
    return s


def roam_sharp(frames, sigmas, t_days, x, y, detail=False):
    """Raw bank maximum (sensitive but not selective at junk sites)."""
    s = roam_bank(frames, sigmas, t_days, x, y)
    i = int(np.argmax(s))
    if detail:
        return float(s[i]), (float(_BANK_SP[i // len(_BANK_AN)]),
                             float(np.degrees(_BANK_AN[i % len(_BANK_AN)])))
    return float(s[i])


def roam_verified(frames, sigmas, t_days, x, y, topk=12, detail=False):
    """Two-stage roam: the fast bank NOMINATES the best few tracks, the full
    careful score (local noise + background + zero-velocity veto +
    concentration + consistency) VERIFIES them, then a small local refinement
    polishes the winner (the bank grid is 25"/d x 2deg, coarser than the
    score's basin). Only the verified value is played as audio -- junk sites
    that fool the quick look cannot fool the full treatment."""
    s = roam_bank(frames, sigmas, t_days, x, y)
    idx = np.argsort(s)[::-1][:topk]
    s_static, _ = raw_s(frames, sigmas, t_days, x, y, 0.0, 0.0, SHARP)
    best, best_v = 0.0, (0.0, 0.0)
    for i in idx:
        if s[i] <= 0:
            break
        sc = score(frames, sigmas, t_days, x, y, _BANK_VX[i], _BANK_VY[i],
                   SHARP, s_static)
        if sc > best:
            best = sc
            best_v = (float(_BANK_SP[i // len(_BANK_AN)]),
                      float(np.degrees(_BANK_AN[i % len(_BANK_AN)])))
    if best > 0:
        for dsp in (-12.5, -6.0, 0.0, 6.0, 12.5):
            for dan in (-1.0, -0.5, 0.0, 0.5, 1.0):
                vx, vy = vel_from_speed_angle(best_v[0] + dsp, best_v[1] + dan)
                sc = score(frames, sigmas, t_days, x, y, vx, vy, SHARP, s_static)
                if sc > best:
                    best = sc
    return (best, best_v) if detail else best


def vel_from_speed_angle(speed, angle_deg):
    """arcsec/day + degrees -> binned px/day."""
    v = speed / PIXEL_SCALE / BIN
    a = np.radians(angle_deg)
    return v * np.cos(a), v * np.sin(a)


ROAM_SPEEDS = list(range(70, 1251, 75))                       # 16 speeds
ROAM_ANGLES = np.arange(0, 360, 5.0)                          # 72 angles


def grid_max(frames, sigmas, t_days, x, y, masks=BROAD, detail=False):
    """Roam-mode signal: best BROAD score over a velocity grid at (x, y).
    The broad core widens the basin to ~+-58 "/day and ~+-4.5 deg so the
    5 deg x 75 "/day grid cannot fall between two basins."""
    s_static, _ = raw_s(frames, sigmas, t_days, x, y, 0.0, 0.0, masks)
    best, best_v = 0.0, (0, 0)
    for sp in ROAM_SPEEDS:
        for ang in ROAM_ANGLES:
            vx, vy = vel_from_speed_angle(sp, ang)
            sc = score(frames, sigmas, t_days, x, y, vx, vy, masks, s_static)
            if sc > best:
                best, best_v = sc, (sp, ang)
    return (best, best_v) if detail else best


# --- ground truth from the pipeline's per-frame positions --------------------

def load_truth(shape_orig):
    H, W = shape_orig
    fp = None
    for f in sorted(glob.glob(os.path.join(WEB_RESULTS, '*_results.json'))):
        d = json.load(open(f))
        if d and 'fpos' in d[0] and d[0].get('name') == '2015 RM287':
            fp = d
            break
    if fp is None:
        raise SystemExit('No web results JSON with fpos found for set203.')
    truths = []
    for c in fp:
        pos = np.array(c['fpos']) * [W, H] / BIN               # binned px, per frame
        truths.append({'name': c['name'], 'rate': c['rate'], 'pos': pos})
    return truths


def main():
    print('Loading + aligning set203 (same astroalign call as the pipeline)...')
    aligned, t_days = load_and_align()
    frames, sigmas, starframe = prepare(aligned)
    H2, W2 = frames[0].shape
    print(f'  {len(frames)} frames, binned {W2}x{H2}, arc {t_days[-1]*24*60:.1f} min, '
          f'sigma(binned) = {np.round(sigmas, 2)}')

    truths = load_truth(aligned[0].shape)
    margin = int(1250 / PIXEL_SCALE / BIN * t_days[-1]) + PATCH + 2

    # 1) True asteroids at their true position+velocity
    print('\n=== 1. TRUE ASTEROIDS (position + velocity from pipeline fpos) ===')
    ast_scores = []
    for tr in truths:
        p0 = tr['pos'][0]
        v = (tr['pos'][-1] - tr['pos'][0]) / t_days[-1]        # binned px/day
        sp = np.hypot(*v) * PIXEL_SCALE * BIN
        ang = np.degrees(np.arctan2(v[1], v[0])) % 360
        sc = score(frames, sigmas, t_days, p0[0], p0[1], v[0], v[1])
        ast_scores.append(sc)
        print(f'  {tr["name"]:12s} rate={tr["rate"]:4d}"/d  (check: {sp:.0f}"/d @ {ang:.0f}deg)'
              f'  score = {sc:.3f}')

    # 2) 50 brightest (unmasked) stars, each given every chance to cheat
    print('\n=== 2. 50 BRIGHTEST STARS x velocity grid (should all be ~0) ===')
    # find stars in the PRE-subtraction frame (stars still present there),
    # and stay away from the known asteroid tracks
    f0 = np.nan_to_num(starframe)
    peaks = (ndimage.maximum_filter(f0, 9) == f0) & (f0 > 20 * sigmas[0])
    peaks[:margin] = peaks[-margin:] = False
    peaks[:, :margin] = peaks[:, -margin:] = False
    ys, xs = np.where(peaks)
    ast_all = np.vstack([tr['pos'] for tr in truths])
    far = np.array([np.min(np.hypot(ast_all[:, 0] - x, ast_all[:, 1] - y)) > 30
                    for x, y in zip(xs, ys)])
    ys, xs = ys[far], xs[far]
    order = np.argsort(f0[ys, xs])[::-1][:50]
    star_broad = []
    for i in order:
        sc, v = grid_max(frames, sigmas, t_days, xs[i], ys[i], detail=True)
        star_broad.append((sc, v, f0[ys[i], xs[i]], xs[i], ys[i]))
    star_broad.sort(key=lambda e: e[0], reverse=True)
    print('  top 5 (broad roam score / winning velocity / star brightness / binned x,y):')
    for sc, (sp, ang), pk, x, y in star_broad[:5]:
        print(f'    score {sc:7.3f}   at {sp:4d}"/d @ {ang:3.0f}deg   peak {pk:7.0f}'
              f'   ({x}, {y})')
    bs = np.array([e[0] for e in star_broad])
    print(f'  broad: max = {bs.max():.3f}  median = {np.median(bs):.3f}  >0.05: {(bs > 0.05).sum()}/50')
    star_sharp = [grid_max(frames, sigmas, t_days, e[3], e[4], masks=SHARP)
                  for e in star_broad[:10]]
    star_sharp_max = max(star_sharp) if star_sharp else 0.0
    print(f'  sharp grid on the 10 worst leakers: max = {star_sharp_max:.3f}')

    # 3) Random sky: 200 single sharp guesses + 30 roam-mode grid-maxes
    print('\n=== 3. RANDOM SKY ===')
    rng = np.random.default_rng(42)
    rx = rng.uniform(margin, W2 - margin, 200)
    ry = rng.uniform(margin, H2 - margin, 200)
    rsc = []
    for i in range(200):
        vx, vy = vel_from_speed_angle(rng.uniform(SPEED_FLOOR, 1250), rng.uniform(0, 360))
        rsc.append(score(frames, sigmas, t_days, rx[i], ry[i], vx, vy))
    rsc = np.array(rsc)
    rgm = np.array([grid_max(frames, sigmas, t_days, rx[i], ry[i]) for i in range(30)])
    print(f'  200 single random sharp guesses: max = {rsc.max():.3f}  median = {np.median(rsc):.3f}')
    print(f'  30 roam-mode grid maxes:         max = {rgm.max():.3f}  median = {np.median(rgm):.3f}')

    # 3b) Roam-mode detectability: does the broad grid CATCH the asteroids?
    print('\n=== 3b. ROAM GRID AT THE TRUE ASTEROID POSITIONS ===')
    ast_roam = []
    for tr in truths:
        p0 = tr['pos'][0]
        sc, (sp, ang) = grid_max(frames, sigmas, t_days, p0[0], p0[1], detail=True)
        ast_roam.append(sc)
        print(f'  {tr["name"]:12s} roam grid-max = {sc:.3f}  (grid hit {sp}"/d @ {ang:.0f}deg,'
              f' true {tr["rate"]}"/d)')

    # 3c) Roam v2: fine-grid sharp matched filter (asteroids vs stars vs random)
    print('\n=== 3c. ROAM v2 (sharp core, fine 25"/d x 2deg grid) ===')
    ast_roam2 = []
    for tr in truths:
        p0 = tr['pos'][0]
        sc, (sp, ang) = roam_sharp(frames, sigmas, t_days, p0[0], p0[1], detail=True)
        ast_roam2.append(sc)
        print(f'  {tr["name"]:12s} roam2 = {sc:.3f}  (hit {sp:.0f}"/d @ {ang:.0f}deg,'
              f' true {tr["rate"]}"/d)')
    star_roam2 = np.array([roam_sharp(frames, sigmas, t_days, e[3], e[4])
                           for e in star_broad[:10]])
    rand_roam2 = np.array([roam_sharp(frames, sigmas, t_days, rx[i], ry[i])
                           for i in range(30)])
    print(f'  10 worst star sites: max = {star_roam2.max():.3f}'
          f'   30 random sites: max = {rand_roam2.max():.3f}'
          f'  median = {np.median(rand_roam2):.3f}')

    # 3d) Roam VERIFIED: bank nominates, full score verifies
    print('\n=== 3d. ROAM VERIFIED (bank nominates top 12, full score verifies + refines) ===')
    ast_roamv = []
    for tr in truths:
        p0 = tr['pos'][0]
        sc, (sp, ang) = roam_verified(frames, sigmas, t_days, p0[0], p0[1], detail=True)
        ast_roamv.append(sc)
        print(f'  {tr["name"]:12s} verified roam = {sc:.3f}  (hit {sp:.0f}"/d @ {ang:.0f}deg,'
              f' true {tr["rate"]}"/d)')
    star_roamv = np.array([roam_verified(frames, sigmas, t_days, e[3], e[4])
                           for e in star_broad[:10]])
    rand_roamv = np.array([roam_verified(frames, sigmas, t_days, rx[i], ry[i])
                           for i in range(30)])
    print(f'  10 worst star sites: max = {star_roamv.max():.3f}'
          f'   30 random sites: max = {rand_roamv.max():.3f}'
          f'  median = {np.median(rand_roamv):.3f}')

    # 4) Basin widths around asteroid #1 (how forgiving is the tuning?)
    print('\n=== 4. BASIN WIDTH around', truths[0]['name'], '===')
    tr = truths[0]
    p0 = tr['pos'][0]
    v = (tr['pos'][-1] - tr['pos'][0]) / t_days[-1]
    sp0 = np.hypot(*v) * PIXEL_SCALE * BIN
    ang0 = np.degrees(np.arctan2(v[1], v[0]))
    print('  speed offset ("/day) -> score:')
    for d in [-100, -50, -25, -10, 0, 10, 25, 50, 100]:
        vx, vy = vel_from_speed_angle(sp0 + d, ang0)
        print(f'    {d:+4d}  {score(frames, sigmas, t_days, p0[0], p0[1], vx, vy):.3f}')
    print('  angle offset (deg) -> score:')
    for d in [-10, -5, -2, 0, 2, 5, 10]:
        vx, vy = vel_from_speed_angle(sp0, ang0 + d)
        print(f'    {d:+4d}  {score(frames, sigmas, t_days, p0[0], p0[1], vx, vy):.3f}')
    print('  position offset (binned px, x) -> score at true velocity:')
    for d in [-8, -4, -2, 0, 2, 4, 8]:
        print(f'    {d:+4d}  {score(frames, sigmas, t_days, p0[0] + d, p0[1], v[0], v[1]):.3f}')

    # 4b) How close must the CURSOR get before verified roam sings on the
    # faintest asteroid? This sets the UI step size / scan line spacing.
    print('\n=== 4b. VERIFIED-ROAM REACH vs cursor offset (faintest asteroid) ===')
    p0 = truths[0]['pos'][0]
    for d in [0, 2, 4, 6, 10]:
        sc = roam_verified(frames, sigmas, t_days, p0[0] + d, p0[1])
        print(f'  cursor {d:+3d} binned px off -> verified roam = {sc:.3f}')

    # Verdict: both channels must separate cleanly.
    print('\n=== VERDICT ===')
    roam_ast = min(ast_roamv)
    roam_fake = max(star_roamv.max(), rand_roamv.max())
    r_roam = roam_ast / roam_fake if roam_fake > 0 else np.inf
    print(f'  ROAM channel (verified): worst asteroid = {roam_ast:.3f}   best star/random = '
          f'{roam_fake:.3f}   ratio = {r_roam:.1f}x')
    tune_ast, tune_fake = min(ast_scores), max(star_sharp_max, rsc.max())
    r_tune = tune_ast / tune_fake if tune_fake > 0 else np.inf
    print(f'  TUNE channel:  worst asteroid = {tune_ast:.3f}   best star/random = '
          f'{tune_fake:.3f}   ratio = {r_tune:.1f}x')
    ok = r_roam > 5 and r_tune > 5
    print('  PASS (both > 5x)' if ok else
          '  FAIL -- tune the formula here before any web code')


if __name__ == '__main__':
    main()
