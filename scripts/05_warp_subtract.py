import numpy as np
import astroalign as aa
import matplotlib.pyplot as plt
from astropy.io import fits

# Load both warps
with fits.open('data/warp1.fits') as hdul:
    warp1 = hdul[0].data.astype(np.float64)
    hdr1 = hdul[0].header

with fits.open('data/warp2.fits') as hdul:
    warp2 = hdul[0].data.astype(np.float64)
    hdr2 = hdul[0].header

# Get MJD (observation date) from headers if available
mjd1 = hdr1.get('MJD-OBS', 'unknown')
mjd2 = hdr2.get('MJD-OBS', 'unknown')
print(f"Warp 1 MJD: {mjd1}")
print(f"Warp 2 MJD: {mjd2}")

# Replace NaN with 0 (warps have lots of masked pixels)
warp1 = np.nan_to_num(warp1)
warp2 = np.nan_to_num(warp2)

# Align warp 2 onto warp 1's pixel grid
print("\nAligning warp2 onto warp1... may take 30 seconds")
try:
    aligned_warp2, footprint = aa.register(warp2, warp1)
    print("Alignment done!")
except Exception as e:
    print(f"Alignment failed: {e}")
    print("This can happen if the warps don't share enough bright stars.")
    print("Trying with relaxed parameters...")
    aligned_warp2, footprint = aa.register(warp2, warp1, max_control_points=100, detection_sigma=3)
    print("Alignment done with relaxed params!")

# Normalize both: subtract median, scale by std
warp1_norm = warp1 - np.median(warp1[warp1 > 0])
aligned_warp2_norm = aligned_warp2 - np.median(aligned_warp2[aligned_warp2 > 0])

# Scale warp 2 to match warp 1's brightness range
scale = np.std(warp1_norm[warp1_norm != 0]) / np.std(aligned_warp2_norm[aligned_warp2_norm != 0])
aligned_warp2_norm = aligned_warp2_norm * scale

# Subtract
diff = aligned_warp2_norm - warp1_norm

print("\nDifference image statistics:")
print(f"  Min: {diff.min():.2f}")
print(f"  Max: {diff.max():.2f}")
print(f"  Mean: {diff.mean():.2f}")
print(f"  Std: {diff.std():.2f}")

# Visualize all four panels
fig, axes = plt.subplots(2, 2, figsize=(12, 12))

v1 = np.percentile(warp1[warp1 > 0], [1, 99])
v2 = np.percentile(aligned_warp2[aligned_warp2 > 0], [1, 99])
diff_range = np.percentile(np.abs(diff[diff != 0]), 99)

axes[0, 0].imshow(warp1, cmap='gray', origin='lower', vmin=v1[0], vmax=v1[1])
axes[0, 0].set_title(f"Warp 1 (MJD {mjd1})")

axes[0, 1].imshow(aligned_warp2, cmap='gray', origin='lower', vmin=v2[0], vmax=v2[1])
axes[0, 1].set_title(f"Warp 2 aligned (MJD {mjd2})")

axes[1, 0].imshow(diff, cmap='RdBu', origin='lower', vmin=-diff_range, vmax=diff_range)
axes[1, 0].set_title("Difference (warp2 - warp1)")

# The footprint shows which pixels in aligned warp2 came from valid data
axes[1, 1].imshow(footprint, cmap='gray', origin='lower')
axes[1, 1].set_title("Footprint (white = valid pixels)")

plt.tight_layout()
plt.savefig('outputs/warp_difference.png', dpi=100)
plt.show()

# Save the difference
fits.writeto('data/warp_difference.fits', diff.astype(np.float32), overwrite=True)
print("\nSaved difference to data/warp_difference.fits")