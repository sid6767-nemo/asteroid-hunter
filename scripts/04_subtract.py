import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits

# Load the reference image and the aligned image
with fits.open('data/image1.fits') as hdul:
    img1 = hdul[0].data.astype(np.float64)

with fits.open('data/image2_aligned.fits') as hdul:
    img2_aligned = hdul[0].data.astype(np.float64)

# Replace any NaN with 0
img1 = np.nan_to_num(img1)
img2_aligned = np.nan_to_num(img2_aligned)

# Normalize both images to a common scale using their median values
# This handles the fact that g and r filters have different brightness
img1_norm = img1 - np.median(img1)
img2_norm = img2_aligned - np.median(img2_aligned)

# Scale image 2 to match image 1's brightness range
scale = np.std(img1_norm) / np.std(img2_norm)
img2_norm = img2_norm * scale

# Subtract
diff = img2_norm - img1_norm

print("Difference image statistics:")
print(f"  Min: {diff.min():.2f}")
print(f"  Max: {diff.max():.2f}")
print(f"  Mean: {diff.mean():.2f}")
print(f"  Std: {diff.std():.2f}")

# Visualize: image 1, image 2 aligned, difference
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Display images with percentile clipping
v1 = np.percentile(img1, [1, 99])
v2 = np.percentile(img2_aligned, [1, 99])
# Difference uses symmetric scale around 0
diff_range = np.percentile(np.abs(diff), 99)

axes[0].imshow(img1, cmap='gray', origin='lower', vmin=v1[0], vmax=v1[1])
axes[0].set_title("Image 1")

axes[1].imshow(img2_aligned, cmap='gray', origin='lower', vmin=v2[0], vmax=v2[1])
axes[1].set_title("Image 2 (aligned)")

axes[2].imshow(diff, cmap='RdBu', origin='lower',
               vmin=-diff_range, vmax=diff_range)
axes[2].set_title("Difference (img2 - img1)")

plt.tight_layout()
plt.savefig('outputs/difference.png', dpi=100)
plt.show()

# Save the difference image
fits.writeto('data/difference.fits', diff.astype(np.float32), overwrite=True)
print("\nSaved difference image to data/difference.fits")