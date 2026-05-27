import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder

# Load difference image
with fits.open('data/warp_difference.fits') as hdul:
    diff = hdul[0].data.astype(np.float64)
diff = np.nan_to_num(diff)

# Detect sources
mean, median, std = sigma_clipped_stats(diff, sigma=3.0)
finder = DAOStarFinder(fwhm=3.0, threshold=5.0 * std)
sources_pos = finder(diff - median)
sources_neg = finder(-(diff - median))

pos_xy = np.array([sources_pos['xcentroid'], sources_pos['ycentroid']]).T
neg_xy = np.array([sources_neg['xcentroid'], sources_neg['ycentroid']]).T
print(f"Positives: {len(pos_xy)}, Negatives: {len(neg_xy)}")

MIN_MOVE = 3
MAX_MOVE = 40

# --- Helper: find the nearest neighbor of a point in a set ---
def nearest(point, candidates):
    """Return (index, distance) of the closest candidate to point."""
    dists = np.hypot(candidates[:, 0] - point[0], candidates[:, 1] - point[1])
    idx = np.argmin(dists)
    return idx, dists[idx]

# --- Mutual nearest neighbor matching ---
# A real candidate: pos's nearest neg is N, AND neg N's nearest pos is that same pos.
# Both must also be within the MIN_MOVE..MAX_MOVE window.

confirmed = []
for i, p in enumerate(pos_xy):
    # nearest negative to this positive
    j, d_pn = nearest(p, neg_xy)
    if not (MIN_MOVE <= d_pn <= MAX_MOVE):
        continue
    # now check: is this positive ALSO the nearest positive to that negative?
    k, d_np = nearest(neg_xy[j], pos_xy)
    if k == i:
        # mutual! they point at each other
        confirmed.append({
            'pos': p,
            'neg': neg_xy[j],
            'dist': d_pn,
        })

print(f"\nMutual-pair candidates (after filter): {len(confirmed)}")
for n, c in enumerate(confirmed):
    print(f"  Candidate {n+1}: moved {c['dist']:.1f} px, "
          f"from ({c['neg'][0]:.0f},{c['neg'][1]:.0f}) "
          f"to ({c['pos'][0]:.0f},{c['pos'][1]:.0f})")

# --- Visualize ---
fig, ax = plt.subplots(figsize=(10, 10))
diff_range = np.percentile(np.abs(diff[diff != 0]), 99)
ax.imshow(diff, cmap='gray', origin='lower', vmin=-diff_range, vmax=diff_range)

for c in confirmed:
    ax.annotate('', xy=c['pos'], xytext=c['neg'],
                arrowprops=dict(arrowstyle='->', color='lime', lw=2))
    ax.scatter(*c['neg'], s=120, edgecolor='red', facecolor='none', linewidth=1.8)
    ax.scatter(*c['pos'], s=120, edgecolor='cyan', facecolor='none', linewidth=1.8)

ax.set_title(f"Confirmed mutual-pair candidates: {len(confirmed)}")
plt.tight_layout()
plt.savefig('outputs/filtered_candidates.png', dpi=100)
plt.show()

print("\nSaved to outputs/filtered_candidates.png")