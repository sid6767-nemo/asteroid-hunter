import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder

# Load difference image
with fits.open('data/warp_difference.fits') as hdul:
    diff = hdul[0].data.astype(np.float64)
diff = np.nan_to_num(diff)

# Background stats
mean, median, std = sigma_clipped_stats(diff, sigma=3.0)
print(f"Background std: {std:.2f}")

# Detect positive and negative sources
finder = DAOStarFinder(fwhm=3.0, threshold=5.0 * std)
sources_pos = finder(diff - median)
sources_neg = finder(-(diff - median))

n_pos = len(sources_pos) if sources_pos is not None else 0
n_neg = len(sources_neg) if sources_neg is not None else 0
print(f"Positive: {n_pos}, Negative: {n_neg}")

# Build coordinate arrays
pos_xy = np.array([sources_pos['xcentroid'], sources_pos['ycentroid']]).T
neg_xy = np.array([sources_neg['xcentroid'], sources_neg['ycentroid']]).T

# --- The matching logic ---
# A moving object = a positive spot with a negative spot nearby,
# but NOT on top of it (that would just be a star that changed brightness).
# We look for pos/neg pairs separated by MIN_MOVE to MAX_MOVE pixels.

MIN_MOVE = 3      # must move at least this far (else it's a static star)
MAX_MOVE = 40     # but not absurdly far (else it's two unrelated objects)

candidates = []
for i, p in enumerate(pos_xy):
    for j, n in enumerate(neg_xy):
        dist = np.hypot(p[0] - n[0], p[1] - n[1])
        if MIN_MOVE <= dist <= MAX_MOVE:
            candidates.append({
                'pos': p,
                'neg': n,
                'dist': dist,
                'midpoint': ((p[0]+n[0])/2, (p[1]+n[1])/2)
            })

print(f"\nCandidate moving objects found: {len(candidates)}")
for k, c in enumerate(candidates):
    print(f"  Candidate {k+1}: moved {c['dist']:.1f} px, "
          f"from ({c['neg'][0]:.0f},{c['neg'][1]:.0f}) "
          f"to ({c['pos'][0]:.0f},{c['pos'][1]:.0f})")

# --- Visualize ---
fig, ax = plt.subplots(figsize=(10, 10))
diff_range = np.percentile(np.abs(diff[diff != 0]), 99)
ax.imshow(diff, cmap='gray', origin='lower', vmin=-diff_range, vmax=diff_range)

# Draw all candidates: line connecting neg->pos, with arrow
for c in candidates:
    ax.annotate('', xy=c['pos'], xytext=c['neg'],
                arrowprops=dict(arrowstyle='->', color='lime', lw=1.5))
    ax.scatter(*c['neg'], s=100, edgecolor='red', facecolor='none', linewidth=1.5)
    ax.scatter(*c['pos'], s=100, edgecolor='cyan', facecolor='none', linewidth=1.5)

ax.set_title(f"Candidate moving objects: {len(candidates)}")
plt.tight_layout()
plt.savefig('outputs/candidates.png', dpi=100)
plt.show()

print("\nSaved to outputs/candidates.png")