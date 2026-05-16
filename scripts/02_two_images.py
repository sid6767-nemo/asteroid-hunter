from astropy.io import fits
import matplotlib.pyplot as plt
import numpy as np

# Load both images
with fits.open('data/image1.fits') as hdul:
    img1 = hdul[0].data
    header1 = hdul[0].header

with fits.open('data/image2.fits') as hdul:
    img2 = hdul[0].data
    header2 = hdul[0].header

# Print info about each
print("=== Image 1 ===")
print("Shape:", img1.shape)
print("Min:", np.nanmin(img1), "Max:", np.nanmax(img1))
print("Filter:", header1.get('FILTER', 'unknown'))

print("\n=== Image 2 ===")
print("Shape:", img2.shape)
print("Min:", np.nanmin(img2), "Max:", np.nanmax(img2))
print("Filter:", header2.get('FILTER', 'unknown'))

# Display them side by side
fig, axes = plt.subplots(1, 2, figsize=(12, 6))

# Percentile clipping so faint stuff is visible
vmin1, vmax1 = np.nanpercentile(img1, [1, 99])
vmin2, vmax2 = np.nanpercentile(img2, [1, 99])

axes[0].imshow(img1, cmap='gray', origin='lower', vmin=vmin1, vmax=vmax1)
axes[0].set_title(f"Image 1 - filter: {header1.get('FILTER', '?')}")

axes[1].imshow(img2, cmap='gray', origin='lower', vmin=vmin2, vmax=vmax2)
axes[1].set_title(f"Image 2 - filter: {header2.get('FILTER', '?')}")

plt.tight_layout()
plt.savefig('outputs/two_images.png', dpi=100)
plt.show()