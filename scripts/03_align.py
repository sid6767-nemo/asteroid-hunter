import numpy as np
import astroalign as aa
import matplotlib.pyplot as plt
from astropy.io import fits

# Load both images
with fits.open('data/image1.fits') as hdul:
    img1 = hdul[0].data.astype(np.float64)

with fits.open('data/image2.fits') as hdul:
    img2 = hdul[0].data.astype(np.float64)

print("Image 1 shape:", img1.shape)
print("Image 2 shape:", img2.shape)

# Replace any NaN values with 0 (alignment doesn't like NaN)
img1 = np.nan_to_num(img1)
img2 = np.nan_to_num(img2)

# Align image 2 onto image 1's pixel grid
print("\nAligning... this can take 10-30 seconds")
aligned_img2, footprint = aa.register(img2, img1)
print("Alignment done!")

# Visualize: image 1, image 2 (original), image 2 (aligned)
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

vmin, vmax = np.percentile(img1, [1, 99])

axes[0].imshow(img1, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
axes[0].set_title("Image 1 (reference)")

axes[1].imshow(img2, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
axes[1].set_title("Image 2 (original)")

axes[2].imshow(aligned_img2, cmap='gray', origin='lower', vmin=vmin, vmax=vmax)
axes[2].set_title("Image 2 (aligned to Image 1)")

plt.tight_layout()
plt.savefig('outputs/alignment.png', dpi=100)
plt.show()

# Save the aligned image for later use
fits.writeto('data/image2_aligned.fits', aligned_img2.astype(np.float32), overwrite=True)
print("\nSaved aligned image to data/image2_aligned.fits")