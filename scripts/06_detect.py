import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder

# Load the difference image from Day 8
with fits.open('data/warp_difference.fits') as hdul:
    diff = hdul[0].data.astype(np.float64)

diff = np.nan_to_num(diff)

# Calculate background statistics (robust to outliers)
mean, median, std = sigma_clipped_stats(diff, sigma=3.0)
print(f"Background mean: {mean:.2f}")
print(f"Background median: {median:.2f}")
print(f"Background std: {std:.2f}")

# DAOStarFinder looks for point sources brighter than a threshold
# We detect on BOTH positive sources (brighter in warp2) and
# negative sources (brighter in warp1) since asteroids appear in both

# Positive detections (object appeared / brighter in warp 2)
finder_pos = DAOStarFinder(fwhm=3.0, threshold=5.0 * std)
sources_pos = finder_pos(diff - median)

# Negative detections (object disappeared / brighter in warp 1)
finder_neg = DAOStarFinder(fwhm=3.0, threshold=5.0 * std)
sources_neg = finder_neg(-(diff - median))

n_pos = len(sources_pos) if sources_pos is not None else 0
n_neg = len(sources_neg) if sources_neg is not None else 0

print(f"\nPositive sources detected: {n_pos}")
print(f"Negative sources detected: {n_neg}")

# Plot the difference image with detected sources marked
fig, ax = plt.subplots(figsize=(10, 10))

diff_range = np.percentile(np.abs(diff[diff != 0]), 99)
ax.imshow(diff, cmap='gray', origin='lower', vmin=-diff_range, vmax=diff_range)

# Mark positive sources in blue circles
if sources_pos is not None:
    ax.scatter(sources_pos['xcentroid'], sources_pos['ycentroid'],
               s=80, edgecolor='cyan', facecolor='none', linewidth=1.2,
               label=f'Positive ({n_pos})')

# Mark negative sources in red circles
if sources_neg is not None:
    ax.scatter(sources_neg['xcentroid'], sources_neg['ycentroid'],
               s=80, edgecolor='red', facecolor='none', linewidth=1.2,
               label=f'Negative ({n_neg})')

ax.set_title("Detected sources in difference image")
ax.legend()
plt.tight_layout()
plt.savefig('outputs/detections.png', dpi=100)
plt.show()

print("\nSaved detections to outputs/detections.png")