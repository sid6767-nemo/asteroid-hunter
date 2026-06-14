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

def flag_for_review(rate):
    """Worth a human's attention if it moves at solar-system speeds.
    Near-stationary things (likely stars/artifacts) are kept but marked
    low-priority, not thrown away."""
    return rate >= 1.0
   
flagged = 0
for name, pixels in candidates.items():
    arcsec_moved = pixels * PIXEL_SCALE
    rate = arcsec_moved / gap_days
    if flag_for_review(rate):
        tag = "** FLAG FOR REVIEW **"
        flagged += 1
    else:
        tag = "(low priority - likely star/artifact, kept anyway)"
    print(f"{name}: {rate:.4f} arcsec/day  ->  {classify(rate)}  {tag}")

print(f"\n{flagged} of {len(candidates)} candidate(s) flagged for review.")

print("\n--- flag logic check on sample rates ---")
for test_rate in [0.02, 5, 120, 800]:
    print(f"{test_rate} arcsec/day -> flagged? {flag_for_review(test_rate)}  ({classify(test_rate)})")