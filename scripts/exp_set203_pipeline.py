# exp_set203_pipeline.py - 4-image asteroid detection pipeline
#
# Run it on the included sample data with no arguments:
#     python scripts/exp_set203_pipeline.py
#
# Or point it at your own data and tweak settings from the terminal:
#     python scripts/exp_set203_pipeline.py --data data/myfield --threshold 3.5
#
# See all options with:
#     python scripts/exp_set203_pipeline.py --help

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

# --- fixed settings (rarely need changing) ---
PIXEL_SCALE      = 0.256   # arcsec/px
MIN_MOVE         = 3       # px: ignore sub-pixel jitter between frames
MAX_MOVE         = 100     # px: ignore absurdly fast matches
CONFIRM_RADIUS   = 5       # px: how close a prediction must land to count as a hit
GIANT_AREA_MIN   = 150     # saturated-blob area to count as a spike-making giant star
GIANT_RADIUS_MIN = 180     # px: smallest giant-star rejection zone
CV_MAX           = 0.35    # max brightness variation across frames (std/mean)


def get_args():
    p = argparse.ArgumentParser(
        description="Detect moving asteroids in a set of telescope images.",
        epilog="Example: python scripts/exp_set203_pipeline.py --data data/set203 --threshold 3.0")
    p.add_argument('--data', default='data/set203',
                   help="folder containing the .fits frames (default: data/set203)")
    p.add_argument('--output', default='outputs',
                   help="folder to save result images (default: outputs)")
    p.add_argument('--threshold', type=float, default=3.0,
                   help="detection threshold in sigma above background (default: 3.0)")
    p.add_argument('--giant-radius-k', type=float, default=10.0,
                   help="size of rejection zones around bright stars; raise to cut "
                        "more star junk, lower to keep more sky (default: 10.0)")
    p.add_argument('--min-frames', type=int, default=4, choices=[3, 4],
                   help="how many of the 4 frames an object must appear in (default: 4)")
    p.add_argument('--rms-max', type=float, default=1.0,
                   help="max straight-line fit error in px; lower = stricter (default: 1.0)")
    p.add_argument('--no-skybot', action='store_true',
                   help="skip the SkyBoT online cross-match (useful offline)")
    return p.parse_args()


def build_wcs(header):
    """Manual WCS from CRPIX/CDELT/CRVAL/CROTA2 (the header has a bad CROTA1
    that breaks astropy's automatic WCS, so we build it ourselves)."""
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
    DATA_DIR     = args.data
    OUTPUT_DIR   = args.output
    THRESH_SIGMA = args.threshold
    GIANT_RADIUS_K = args.giant_radius_k
    MIN_FRAMES   = args.min_frames
    RMS_MAX      = args.rms_max
    DO_SKYBOT    = not args.no_skybot

    os.makedirs(OUTPUT_DIR, exist_ok=True)
    dataset_name = os.path.basename(os.path.normpath(DATA_DIR))  # e.g. "set203"

    # --- load all frames, SORTED BY TIME (not filename) ---
    pattern = os.path.join(DATA_DIR, '*.fits')
    files = sorted(glob.glob(pattern),
                   key=lambda f: fits.open(f)[0].header['MJD-OBS'])
    if len(files) < 4:
        print(f"ERROR: found {len(files)} .fits files in '{DATA_DIR}', need at least 4.")
        return
    files = files[:4]
    print(f"Loading {len(files)} frames from '{DATA_DIR}' (chronological order)...")

    images, mjds, headers = [], [], []
    for f in files:
        with fits.open(f) as hdul:
            images.append(hdul[0].data.astype(float))
            mjds.append(hdul[0].header['MJD-OBS'])
            headers.append(hdul[0].header)
            print(f"  {os.path.basename(f)[:24]}...  MJD={mjds[-1]:.6f}")

    gaps = [(mjds[i+1] - mjds[i]) * 24 * 60 for i in range(3)]
    t_days = [mjds[i] - mjds[0] for i in range(4)]
    print(f"\nTime gaps: {gaps[0]:.1f}, {gaps[1]:.1f}, {gaps[2]:.1f} min")

    # --- align frames 2,3,4 to frame 1 ---
    print("\nAligning frames to frame 1...")
    aligned = [images[0]]
    for i in range(1, 4):
        try:
            registered, _ = aa.register(images[i], images[0], detection_sigma=5)
            aligned.append(registered)
            print(f"  Frame {i+1} aligned OK")
        except Exception as e:
            print(f"  Frame {i+1} FAILED: {e}")
            aligned.append(images[i])

    # --- background, threshold, saturation level (from frame 1) ---
    mask0 = (aligned[0] == 0)
    _, ref_median, ref_std = sigma_clipped_stats(aligned[0], sigma=3.0, mask=mask0)
    THRESHOLD = THRESH_SIGMA * ref_std
    valid_pixels = aligned[0][aligned[0] > 0]
    img_median = np.median(valid_pixels)
    SATURATION_LEVEL = img_median + 200 * ref_std
    print(f"\nref_std={ref_std:.1f}  threshold={THRESHOLD:.1f}  saturation={SATURATION_LEVEL:.0f}")

    # --- build rejection zones ONLY around giant (spike-making) stars ---
    lbl, nblob = ndimage.label(aligned[0] > SATURATION_LEVEL)
    sizes = ndimage.sum(np.ones_like(lbl), lbl, range(1, nblob + 1))
    giant_zones = []
    for i in range(1, nblob + 1):
        area = sizes[i - 1]
        if area >= GIANT_AREA_MIN:
            ys, xs = np.where(lbl == i)
            radius = max(GIANT_RADIUS_MIN, GIANT_RADIUS_K * np.sqrt(area))
            giant_zones.append((xs.mean(), ys.mean(), radius))
    print(f"{len(giant_zones)} giant-star rejection zone(s)")

    def in_giant_zone(x, y):
        return any(np.hypot(x - cx, y - cy) < r for cx, cy, r in giant_zones)

    # --- detect sources in each aligned frame (keep centroid AND peak) ---
    print("\nDetecting sources (segmentation + deblending)...")
    all_pos, all_peak = [], []
    for i, img in enumerate(aligned):
        mask = (img == 0)
        _, median, _ = sigma_clipped_stats(img, sigma=3.0, mask=mask)
        data = img - median
        segm = detect_sources(data, THRESHOLD, n_pixels=5, mask=mask)
        if segm is None:
            all_pos.append(np.empty((0, 2))); all_peak.append(np.empty(0))
            print(f"  Frame {i+1}: 0 sources"); continue
        segm = deblend_sources(data, segm, n_pixels=5, n_levels=32, contrast=0.0001)
        cat = SourceCatalog(data, segm, mask=mask).to_table()
        pos  = np.array([[float(r['x_centroid']), float(r['y_centroid'])] for r in cat])
        peak = np.array([float(r['max_value']) for r in cat])
        all_pos.append(pos); all_peak.append(peak)
        print(f"  Frame {i+1}: {len(pos)} sources")

    trees = [KDTree(p) if len(p) else None for p in all_pos]

    def mutual_pairs(iA, iB):
        A, B = all_pos[iA], all_pos[iB]
        if len(A) == 0 or len(B) == 0: return []
        tA, tB = KDTree(A), KDTree(B)
        df, jf = tB.query(A); _, jb = tA.query(B)
        out = []
        for i, (d, j) in enumerate(zip(df, jf)):
            if MIN_MOVE < d < MAX_MOVE and jb[j] == i:
                out.append((A[i], B[j], d))
        return out

    confirmed = []

    def track_quality(fpos, fpk):
        ts = np.array(t_days)
        xs = np.array([p[0] for p in fpos]); ys = np.array([p[1] for p in fpos])
        px = np.polyfit(ts, xs, 1); py = np.polyfit(ts, ys, 1)
        rms = np.sqrt(np.mean((xs - np.polyval(px, ts))**2 + (ys - np.polyval(py, ts))**2))
        speed = np.hypot(px[0], py[0]) * PIXEL_SCALE
        pk = np.array(fpk); pk = pk[~np.isnan(pk)]
        cv = (pk.std() / pk.mean()) if len(pk) >= 2 and pk.mean() > 0 else 9.9
        return rms, cv, speed

    def line_check_and_add(pairs, start_gap):
        for A, B, dist in pairs:
            step = B - A
            frame_pred = [A - step * start_gap + step * i for i in range(4)]
            fpos, fpk, hits = [], [], 0
            for fi, (tree, pred) in enumerate(zip(trees, frame_pred)):
                if tree is None:
                    fpos.append(pred); fpk.append(np.nan); continue
                d, idx = tree.query(pred)
                if d < CONFIRM_RADIUS:
                    hits += 1
                    fpos.append(all_pos[fi][idx])
                    fpk.append(all_peak[fi][idx])
                else:
                    fpos.append(pred); fpk.append(np.nan)
            if hits < 3:
                continue
            f1 = fpos[0]
            if any(in_giant_zone(p[0], p[1]) for p in fpos):
                continue
            if hits < MIN_FRAMES:
                continue
            rms, cv, speed = track_quality(fpos, fpk)
            if rms > RMS_MAX or cv > CV_MAX:
                continue
            if any(np.hypot(f1[0]-e['f1'][0], f1[1]-e['f1'][1]) < 30 for e in confirmed):
                continue
            confirmed.append({'f1': f1, 'fpos': fpos, 'fpk': fpk, 'hits': hits,
                              'rate': speed, 'rms': rms, 'cv': cv})
            print(f"  CONFIRMED #{len(confirmed)}: {speed:.0f}\"/day  "
                  f"({f1[0]:.0f},{f1[1]:.0f})  {hits}/4  RMS={rms:.2f} CV={cv:.2f}")

    print("\nLinking + filtering tracks...")
    line_check_and_add(mutual_pairs(0, 1), 0)
    line_check_and_add(mutual_pairs(1, 2), 1)
    line_check_and_add(mutual_pairs(2, 3), 2)
    print(f"\n{len(confirmed)} confirmed candidate(s) after filtering")

    print("\n=== FINAL CONFIRMED CANDIDATES ===")
    for n, c in enumerate(confirmed):
        print(f"\nCandidate #{n+1}")
        print(f"  Speed:     {c['rate']:.0f} arcsec/day  ->  {classify(c['rate'])}")
        print(f"  Frame-1:   ({c['f1'][0]:.0f},{c['f1'][1]:.0f})")
        print(f"  Quality:   {c['hits']}/4 frames, linRMS={c['rms']:.2f}px, brightnessCV={c['cv']:.2f}")

    # --- visualization ---
    print("\nGenerating visualization...")
    fig, axes = plt.subplots(1, 4, figsize=(20, 6))
    fig.suptitle(f'Confirmed moving objects: {len(confirmed)}', fontsize=14)
    colors = ['cyan', 'orange', 'red', 'lime', 'magenta', 'yellow']
    for frame_idx, (ax, img) in enumerate(zip(axes, aligned)):
        _, med, std = sigma_clipped_stats(img, sigma=3.0)
        ax.imshow(img, cmap='gray', vmin=med-2*std, vmax=med+4*std, origin='upper')
        ax.set_title(f'Frame {frame_idx+1}')
        for n, c in enumerate(confirmed):
            pos = c['fpos'][frame_idx]
            col = colors[n % len(colors)]
            ax.add_patch(plt.Circle((pos[0], pos[1]), 15, color=col, fill=False, lw=2))
            ax.text(pos[0]+18, pos[1]+18, f'#{n+1}', color=col, fontsize=8, fontweight='bold')
    plt.tight_layout()
    out_png = os.path.join(OUTPUT_DIR, f'{dataset_name}_tracks.png')
    plt.savefig(out_png, dpi=150, bbox_inches='tight')
    print(f"Saved {out_png}")
    plt.show()
    # --- SkyBoT cross-match ---
    if not DO_SKYBOT:
        print("\n(SkyBoT cross-match skipped via --no-skybot)")
        print("\nDone.")
        return

    from astropy.coordinates import SkyCoord
    from astropy.time import Time
    import astropy.units as u
    try:
        from astroquery.imcce import Skybot
    except Exception:
        print("\n(astroquery not available - skipping SkyBoT)")
        print("\nDone.")
        return

    w = build_wcs(headers[0])
    epoch = Time(mjds[0], format='mjd')
    print("\nCross-matching candidates against SkyBoT...")
    for n, c in enumerate(confirmed):
        ra, dec = w.all_pix2world(c['f1'][0], c['f1'][1], 0)
        ra, dec = float(ra), float(dec)
        print(f"\n--- Candidate #{n+1}: RA={ra:.5f} Dec={dec:.5f} ({c['rate']:.0f}\"/day) ---")
        try:
            res = Skybot.cone_search(SkyCoord(ra=ra, dec=dec, unit='deg'), 0.05*u.deg, epoch)
            print(res)
        except Exception as e:
            print(f"  No SkyBoT match / query failed: {e}")

    print("\nDone.")


if __name__ == '__main__':
    main()