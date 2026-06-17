# exp_set203_pipeline.py - 4-image asteroid detection on IASC set203
import glob
import numpy as np
from astropy.io import fits
import astroalign as aa
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from scipy.spatial import KDTree

PIXEL_SCALE = 0.256  # arcsec/px

# --- load all 4 frames ---
files = sorted(glob.glob('data/set203/*.fits'))
print(f"Loading {len(files)} frames...")

images = []
mjds = []
for f in files:
    with fits.open(f) as hdul:
        data = hdul[0].data.astype(float)
        mjd = hdul[0].header['MJD-OBS']
        images.append(data)
        mjds.append(mjd)
        print(f"  {f.split('/')[-1][:20]}...  MJD={mjd:.6f}")

# time gaps in minutes
gaps = [(mjds[i+1] - mjds[i]) * 24 * 60 for i in range(3)]
print(f"\nTime gaps: {gaps[0]:.1f} min, {gaps[1]:.1f} min, {gaps[2]:.1f} min")

# --- align frames 2,3,4 to frame 1 ---
print("\nAligning frames to frame 1...")
aligned = [images[0]]
for i in range(1, 4):
    try:
        registered, _ = aa.register(images[i], images[0])
        aligned.append(registered)
        print(f"  Frame {i+1} aligned OK")
    except Exception as e:
        print(f"  Frame {i+1} FAILED: {e}")

print(f"\nAligned {len(aligned)} of 4 frames successfully.")

# --- detect sources in each aligned frame ---
mask0 = (aligned[0] == 0)
_, ref_median, ref_std = sigma_clipped_stats(aligned[0], sigma=3.0, mask=mask0)
THRESHOLD = 4.2 * ref_std
print(f"\nReference std={ref_std:.1f}, global threshold={THRESHOLD:.1f} counts")

print("\nDetecting sources in each frame...")
all_sources = []
for i, img in enumerate(aligned):
    mask = (img == 0)
    _, median, _ = sigma_clipped_stats(img, sigma=3.0, mask=mask)
    dao = DAOStarFinder(fwhm=5.0, threshold=THRESHOLD)
    sources = dao(img - median, mask=mask)
    if sources is None:
        sources = []
        print(f"  Frame {i+1}: 0 sources")
    else:
        print(f"  Frame {i+1}: {len(sources)} sources")
    all_sources.append(sources)

print("\nColumn names:", all_sources[0].colnames)

# --- find candidates between frame 1 and frame 2 ---
MIN_MOVE = 3    # px
MAX_MOVE = 100  # px

print("\nFinding candidates (frame 1 -> frame 2)...")

pos1 = np.array([[s['x_centroid'], s['y_centroid']] for s in all_sources[0]])
pos2 = np.array([[s['x_centroid'], s['y_centroid']] for s in all_sources[1]])

tree2 = KDTree(pos2)
dists, idxs = tree2.query(pos1)

candidates = []
for i, (dist, j) in enumerate(zip(dists, idxs)):
    if MIN_MOVE < dist < MAX_MOVE:
        candidates.append({'pos1': pos1[i], 'pos2': pos2[j], 'dist': dist})

print(f"  {len(candidates)} candidate(s) between frame 1 and 2")
for n, c in enumerate(candidates):
    rate = c['dist'] * PIXEL_SCALE / (gaps[0] / (24*60))
    print(f"  #{n+1}: {c['dist']:.1f} px = {rate:.0f} arcsec/day  "
          f"({c['pos1'][0]:.0f},{c['pos1'][1]:.0f}) -> ({c['pos2'][0]:.0f},{c['pos2'][1]:.0f})")
    # --- line check: confirm candidates in frames 3 and 4 ---
CONFIRM_RADIUS = 5  # px - how close the source needs to be to predicted position

pos3 = np.array([[s['x_centroid'], s['y_centroid']] for s in all_sources[2]])
pos4 = np.array([[s['x_centroid'], s['y_centroid']] for s in all_sources[3]])

tree3 = KDTree(pos3)
tree4 = KDTree(pos4)

print("\nRunning line-check on candidates...")
confirmed = []
for n, c in enumerate(candidates):
    A = c['pos1']
    B = c['pos2']
    step = B - A

    # predict positions in frames 3 and 4
    pred3 = B + step
    pred4 = B + 2 * step

    # check if a source exists near predicted position
    dist3, _ = tree3.query(pred3)
    dist4, _ = tree4.query(pred4)

    hits = 0
    if dist3 < CONFIRM_RADIUS:
        hits += 1
    if dist4 < CONFIRM_RADIUS:
        hits += 1

    if hits >= 1:  # at least 3 out of 4 frames
        rate = c['dist'] * PIXEL_SCALE / (gaps[0] / (24*60))
        confirmed.append({**c, 'pred3': pred3, 'pred4': pred4,
                         'dist3': dist3, 'dist4': dist4, 'hits': hits})
        print(f"  CONFIRMED #{len(confirmed)}: {c['dist']:.1f} px = {rate:.0f} arcsec/day "
              f"({A[0]:.0f},{A[1]:.0f})->({B[0]:.0f},{B[1]:.0f}) "
              f"| frame3: {dist3:.1f}px off | frame4: {dist4:.1f}px off "
              f"| {hits+2}/4 frames")

print(f"\n{len(confirmed)} candidate(s) survived line-check out of {len(candidates)}")