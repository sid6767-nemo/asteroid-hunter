"""
exp_artifact_filter.py — experiment: reject candidates near saturated/bright regions.
Reuses the pipeline, then checks each candidate against the original images.
"""

import sys
import numpy as np
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
import astroalign as aa

MIN_MOVE = 8
MAX_MOVE = 40

def load_fits(path):
    with fits.open(path) as hdul:
        data = hdul[0].data.astype(np.float64)
    return np.nan_to_num(data)

def detect_sources(image):
    mean, median, std = sigma_clipped_stats(image, sigma=3.0)
    finder = DAOStarFinder(fwhm=3.0, threshold=5.0 * std)
    s = finder(image - median)
    if s is None:
        return np.empty((0, 2))
    return np.array([s['xcentroid'], s['ycentroid']]).T

def nearest(point, cands):
    if len(cands) == 0:
        return None, np.inf
    d = np.hypot(cands[:,0]-point[0], cands[:,1]-point[1])
    i = np.argmin(d)
    return i, d[i]

def find_candidates(pos, neg):
    out = []
    for i, p in enumerate(pos):
        j, dpn = nearest(p, neg)
        if j is None or not (MIN_MOVE <= dpn <= MAX_MOVE):
            continue
        k, _ = nearest(neg[j], pos)
        if k == i:
            out.append({'pos': p, 'neg': neg[j], 'dist': dpn})
    return out

def local_brightness(image, x, y, box=10):
    """Return the max pixel value in a box around (x,y) in the image."""
    h, w = image.shape
    x0, x1 = max(0, int(x)-box), min(w, int(x)+box)
    y0, y1 = max(0, int(y)-box), min(h, int(y)+box)
    region = image[y0:y1, x0:x1]
    if region.size == 0:
        return 0
    return np.max(region)

def main():
    img1 = load_fits(sys.argv[1])
    img2 = load_fits(sys.argv[2])
    aligned2 = aa.register(img2, img1)[0]

    diff = (aligned2 - np.median(aligned2[aligned2!=0])) - (img1 - np.median(img1[img1!=0]))

    pos = detect_sources(diff)
    neg = detect_sources(-diff)
    cands = find_candidates(pos, neg)
    print(f"Candidates before artifact filter: {len(cands)}")

    # Saturation threshold: how bright is "too bright" (near a saturated star)
    # Use a high percentile of the original image as the cutoff
    sat_level = np.percentile(img1[img1 > 0], 99.5)
    print(f"Saturation cutoff (99.5th percentile of img1): {sat_level:.0f}")

    kept = []
    for n, c in enumerate(cands):
        x, y = c['pos']
        b1 = local_brightness(img1, x, y)
        b2 = local_brightness(aligned2, x, y)
        near_sat = (b1 > sat_level) or (b2 > sat_level)
        status = "REJECT (near saturation)" if near_sat else "KEEP"
        print(f"  #{n+1} at ({x:.0f},{y:.0f}) moved {c['dist']:.1f}px "
              f"| local max img1={b1:.0f} img2={b2:.0f} | {status}")
        if not near_sat:
            kept.append(c)

    print(f"\nCandidates after artifact filter: {len(kept)}")

if __name__ == "__main__":
    main()