# exp_motion_rate.py
# Day 15 experiment: turn pixel motion into real sky motion (arcsec/day)
from astropy.io import fits

PIXEL_SCALE = 0.25  # arcsec per pixel, from CDELT1 in the header

# read the exposure time (MJD-OBS) out of each image's header
h1 = fits.open('data/warp1.fits')[0].header
h2 = fits.open('data/warp2.fits')[0].header
gap_days = abs(h2['MJD-OBS'] - h1['MJD-OBS'])
print(f"Time between the two images: {gap_days:.2f} days")
print()

# the three final candidates and how far each moved, in pixels
candidates = {'#1': 39.0, '#2': 34.0, '#3': 32.4}

def classify(rate):
    """First-pass classification by apparent motion (arcsec/day). Approximate."""
    if rate < 1:
        return "too slow → proper-motion star or artifact, NOT a solar-system object"
    elif rate < 50:
        return "slow → distant-object / TNO-like range"
    elif rate < 500:
        return "main-belt asteroid range"
    else:
        return "fast → NEO-like"
    
for name, pixels in candidates.items():
    arcsec_moved = pixels * PIXEL_SCALE
    rate = arcsec_moved / gap_days
    print(f"{name}: {rate:.4f} arcsec/day  →  {classify(rate)}")