# exp_set203_pipeline.py - 4-image asteroid detection on IASC set203
import glob
import numpy as np
from astropy.io import fits
import astroalign as aa
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from scipy.spatial import KDTree
import matplotlib.pyplot as plt

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

# --- find candidates between frame 1 and frame 2 ---
MIN_MOVE = 3
MAX_MOVE = 100

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

# --- line check ---
CONFIRM_RADIUS = 5

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
    pred3 = B + step
    pred4 = B + 2 * step
    dist3, _ = tree3.query(pred3)
    dist4, _ = tree4.query(pred4)
    hits = 0
    if dist3 < CONFIRM_RADIUS:
        hits += 1
    if dist4 < CONFIRM_RADIUS:
        hits += 1
    if hits >= 1:
        rate = c['dist'] * PIXEL_SCALE / (gaps[0] / (24*60))
        confirmed.append({**c, 'pred3': pred3, 'pred4': pred4,
                         'dist3': dist3, 'dist4': dist4, 'hits': hits})
        print(f"  CONFIRMED #{len(confirmed)}: {c['dist']:.1f} px = {rate:.0f} arcsec/day "
              f"({A[0]:.0f},{A[1]:.0f})->({B[0]:.0f},{B[1]:.0f}) "
              f"| frame3: {dist3:.1f}px off | frame4: {dist4:.1f}px off "
              f"| {hits+2}/4 frames")

def classify(rate):
    if rate < 1:
        return "too slow -> likely star or artifact"
    elif rate < 50:
        return "slow -> TNO-like"
    elif rate < 500:
        return "main-belt asteroid range"
    else:
        return "fast -> NEO-like"

print(f"\n{len(confirmed)} candidate(s) survived line-check out of {len(candidates)}")
print("\n=== FINAL CONFIRMED CANDIDATES ===")
for n, c in enumerate(confirmed):
    rate = c['dist'] * PIXEL_SCALE / (gaps[0] / (24*60))
    A, B = c['pos1'], c['pos2']
    print(f"\nCandidate #{n+1}")
    print(f"  Speed:     {rate:.0f} arcsec/day  ->  {classify(rate)}")
    print(f"  Track:     ({A[0]:.0f},{A[1]:.0f}) -> ({B[0]:.0f},{B[1]:.0f})")
    print(f"  Frame 3:   {c['dist3']:.1f}px from predicted position")
    print(f"  Frame 4:   {c['dist4']:.1f}px from predicted position")
    print(f"  Confirmed: {c['hits']+2}/4 frames")

# --- visualize all 4 frames with confirmed candidates ---
print("\nGenerating visualization...")
fig, axes = plt.subplots(1, 4, figsize=(20, 6))
fig.suptitle(f'Confirmed moving objects: {len(confirmed)}', fontsize=14)

colors = ['cyan', 'orange', 'red']
for frame_idx, (ax, img) in enumerate(zip(axes, aligned)):
    mean, med, std = sigma_clipped_stats(img, sigma=3.0)
    ax.imshow(img, cmap='gray', vmin=med-2*std, vmax=med+4*std, origin='lower')
    ax.set_title(f'Frame {frame_idx+1}')
    for n, c in enumerate(confirmed):
        A, B = c['pos1'], c['pos2']
        step = B - A
        color = colors[n % len(colors)]
        if frame_idx == 0:
            pos = A
        elif frame_idx == 1:
            pos = B
        elif frame_idx == 2:
            pos = B + step
        else:
            pos = B + 2 * step
        circle = plt.Circle((pos[0], pos[1]), 15,
                            color=color, fill=False, linewidth=2)
        ax.add_patch(circle)
        ax.text(pos[0]+18, pos[1]+18, f'#{n+1}',
               color=color, fontsize=8, fontweight='bold')

plt.tight_layout()
plt.savefig('outputs/set203_tracks.png', dpi=150, bbox_inches='tight')
print("Saved to outputs/set203_tracks.png")
plt.show()

# --- zoom in on candidate #3 ---
fig2, axes2 = plt.subplots(1, 4, figsize=(20, 5))
fig2.suptitle('Candidate #3 zoom (739 arcsec/day - NEO-like) 4/4 frames', fontsize=13)

c = confirmed[2]
A = c['pos1']
B = c['pos2']
step = B - A
positions = [A, B, B+step, B+2*step]
ZOOM = 80

for frame_idx, (ax, img, pos) in enumerate(zip(axes2, aligned, positions)):
    mean, med, std = sigma_clipped_stats(img, sigma=3.0)
    x, y = int(pos[0]), int(pos[1])
    x1, x2 = max(0, x-ZOOM), min(img.shape[1], x+ZOOM)
    y1, y2 = max(0, y-ZOOM), min(img.shape[0], y+ZOOM)
    cutout = img[y1:y2, x1:x2]
    ax.imshow(cutout, cmap='gray', vmin=med-2*std, vmax=med+4*std, origin='lower')
    ax.set_title(f'Frame {frame_idx+1}  ({x},{y})')
    circle = plt.Circle((pos[0]-x1, pos[1]-y1), 10,
                        color='red', fill=False, linewidth=2)
    ax.add_patch(circle)

plt.tight_layout()
plt.savefig('outputs/set203_candidate3_zoom.png', dpi=150, bbox_inches='tight')
print("Saved to outputs/set203_candidate3_zoom.png")
plt.show()