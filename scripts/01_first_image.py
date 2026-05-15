from astropy.io import fits
from astropy.utils.data import get_pkg_data_filename
import matplotlib.pyplot as plt

# Get a sample astronomy image that ships with astropy
filename = get_pkg_data_filename('tutorials/FITS-images/HorseHead.fits')

# Open the FITS file
with fits.open(filename) as hdul:
    image_data = hdul[0].data

# Show some basic info
print("Image shape:", image_data.shape)
print("Min pixel value:", image_data.min())
print("Max pixel value:", image_data.max())

# Display the image
plt.imshow(image_data, cmap='gray', origin='lower')
plt.colorbar()
plt.title('Horsehead Nebula - first FITS image')
plt.savefig('outputs/horsehead.png')
plt.show()