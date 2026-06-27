# exp_set203_pipeline.py - 4-image asteroid detection on IASC set203
# Full pipeline: align -> detect -> link tracks -> filter -> cross-match.
#
# Everything is wired in. To use: replace your old file with this and run.
# Tunable knobs are in the CONFIG block right below.

import glob
import numpy as np
from astropy.io import fits
import astroalign as aa
from astropy.stats import sigma_clipped_stats
from scipy.spatial import KDTree
from scipy import ndimage
from photutils.segmentation import detect_sources, deblend_sources, SourceCatalog
import matplotlib.pyplot as plt

# ============================ CONFIG ========================================
PIXEL_SCALE      = 0.256   # arcsec/px
THRESH_SIGMA     = 3.0     # detection threshold (sigma above background)
MIN_MOVE         = 3       # px: ignore sub-pixel jitter between frames
MAX_MOVE         = 100     # px: ignore absurdly fast matches
CONFIRM_RADIUS   = 5       # px: how close a prediction must land to count as a hit

# --- giant-star spike rejection (replaces the old 200px blanket filter) ---
GIANT_AREA_MIN   = 150     # saturated-blob area to count as a spike-making giant
GIANT_RADIUS_K   = 10.0     # zone radius = K*sqrt(area)...
GIANT_RADIUS_MIN = 180     # ...floored at this many px

# --- candidate quality gates ---
RMS_MAX          = 1.0     # max linear-fit residual (px). real tracks: ~0.1-0.6
CV_MAX           = 0.35    # max brightness variation (std/mean of peak across frames)
REQUIRE_ALL_FOUR = True    # True = require 4/4 frames (clean). False = keep 3/4 too.
# ============================================================================


def build_wcs(header):
    """Manual WCS from CRPIX/CDELT/CRVAL/CROTA2 (header has a bad CROTA1 that
    breaks astropy's automatic WCS, so we build it ourselves)."""
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


# --- load all 4 frames, SORTED BY TIME (not filename) ---
files = sorted(glob.glob('data/set203/*.fits'),
               key=lambda f: fits.open(f)[0].header['MJD-OBS'])
print(f"Loading {len(files)} frames (chronological order)...")

images, mjds, headers = [], [], []
for f in files:
    with fits.open(f) as hdul:
        images.append(hdul[0].data.astype(float))
        mjds.append(hdul[0].header['MJD-OBS'])
        headers.append(hdul[0].header)
        print(f"  {f.split('/')[-1][:24]}...  MJD={mjds[-1]:.6f}")

gaps = [(mjds[i+1] - mjds[i]) * 24 * 60 for i in range(3)]   # minutes
t_days = [mjds[i] - mjds[0] for i in range(4)]               # offset from frame 1, days
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
giant_zones = []   # (cx, cy, radius)
for i in range(1, nblob + 1):
    area = sizes[i - 1]
    if area >= GIANT_AREA_MIN:
        ys, xs = np.where(lbl == i)
        radius = max(GIANT_RADIUS_MIN, GIANT_RADIUS_K * np.sqrt(area))
        giant_zones.append((xs.mean(), ys.mean(), radius))
print(f"{len(giant_zones)} giant-star rejection zone(s):")
for cx, cy, r in giant_zones:
    print(f"  ({cx:.0f},{cy:.0f}) r={r:.0f}px")

def in_giant_zone(x, y):
    return any(np.hypot(x - cx, y - cy) < r for cx, cy, r in giant_zones)

# --- detect sources in each aligned frame (keep centroid AND peak) ---
print("\nDetecting sources (segmentation + deblending)...")
all_pos, all_peak = [], []
for i, img in enumerate(aligned):
    mask = (img == 0)
    _, median, _ = sigma_clipped_stats(img, sigma=3.0, mask=mask)
    data = img - median
    segm = detect_sources(data, THRESHOLD, npixels=5, mask=mask)
    if segm is None:
        all_pos.append(np.empty((0, 2))); all_peak.append(np.empty(0))
        print(f"  Frame {i+1}: 0 sources"); continue
    segm = deblend_sources(data, segm, npixels=5, nlevels=32, contrast=0.0001)
    cat = SourceCatalog(data, segm, mask=mask).to_table()
    pos  = np.array([[float(r['x_centroid']), float(r['y_centroid'])] for r in cat])
    peak = np.array([float(r['max_value']) for r in cat])
    all_pos.append(pos); all_peak.append(peak)
    print(f"  Frame {i+1}: {len(pos)} sources")

trees = [KDTree(p) if len(p) else None for p in all_pos]

# --- mutual nearest-neighbour candidate pairs for a given frame pair ---
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

# --- linear-track confirmation + quality filtering ---
confirmed = []

def track_quality(fpos, fpk):
    ts = np.array(t_days)
    xs = np.array([p[0] for p in fpos]); ys = np.array([p[1] for p in fpos])
    px = np.polyfit(ts, xs, 1); py = np.polyfit(ts, ys, 1)
    rms = np.sqrt(np.mean((xs - np.polyval(px, ts))**2 + (ys - np.polyval(py, ts))**2))
    speed = np.hypot(px[0], py[0]) * PIXEL_SCALE        # arcsec/day
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
                fpos.append(all_pos[fi][idx])      # use the REAL detected position
                fpk.append(all_peak[fi][idx])
            else:
                fpos.append(pred); fpk.append(np.nan)
        if hits < 3:
            continue

        f1 = fpos[0]
        # --- quality gates (replaces old bright-star filter) ---
        if any(in_giant_zone(p[0], p[1]) for p in fpos):
            continue
        if REQUIRE_ALL_FOUR and hits < 4:
            continue
        rms, cv, speed = track_quality(fpos, fpk)
        if rms > RMS_MAX or cv > CV_MAX:
            continue
        # dedup against already-confirmed tracks
        if any(np.hypot(f1[0]-e['f1'][0], f1[1]-e['f1'][1]) < 30 for e in confirmed):
            continue

        confirmed.append({'f1': f1, 'fpos': fpos, 'fpk': fpk, 'hits': hits,
                          'rate': speed, 'rms': rms, 'cv': cv, 'start_gap': start_gap})
        print(f"  CONFIRMED #{len(confirmed)}: {speed:.0f}\"/day  "
              f"({f1[0]:.0f},{f1[1]:.0f})  {hits}/4  RMS={rms:.2f} CV={cv:.2f}")

print("\nLinking + filtering tracks...")
line_check_and_add(mutual_pairs(0, 1), 0)
line_check_and_add(mutual_pairs(1, 2), 1)
line_check_and_add(mutual_pairs(2, 3), 2)
print(f"\n{len(confirmed)} confirmed candidate(s) after filtering")

# --- summary ---
print("\n=== FINAL CONFIRMED CANDIDATES ===")
for n, c in enumerate(confirmed):
    print(f"\nCandidate #{n+1}")
    print(f"  Speed:     {c['rate']:.0f} arcsec/day  ->  {classify(c['rate'])}")
    print(f"  Frame-1:   ({c['f1'][0]:.0f},{c['f1'][1]:.0f})")
    print(f"  Quality:   {c['hits']}/4 frames, linRMS={c['rms']:.2f}px, brightnessCV={c['cv']:.2f}")

# --- visualization (uses start_gap so circles land correctly for every track) ---
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
plt.savefig('outputs/set203_tracks.png', dpi=150, bbox_inches='tight')
print("Saved outputs/set203_tracks.png")
plt.show()
# --- convert to RA/Dec + SkyBoT cross-match ---
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
try:
    from astroquery.imcce import Skybot
    have_skybot = True
except Exception:
    have_skybot = False
    print("\n(astroquery not available - skipping SkyBoT)")

w = build_wcs(headers[0])
epoch = Time(mjds[0], format='mjd')
print("\nCross-matching candidates against SkyBoT...")
for n, c in enumerate(confirmed):
    ra, dec = w.all_pix2world(c['f1'][0], c['f1'][1], 0)
    ra, dec = float(ra), float(dec)
    print(f"\n--- Candidate #{n+1}: RA={ra:.5f} Dec={dec:.5f} ({c['rate']:.0f}\"/day) ---")
    if not have_skybot:
        continue
    try:
        res = Skybot.cone_search(SkyCoord(ra=ra, dec=dec, unit='deg'), 0.05*u.deg, epoch)
        print(res)
    except Exception as e:
        print(f"  No SkyBoT match / query failed: {e}")

print("\nDone.")