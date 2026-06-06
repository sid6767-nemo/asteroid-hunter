"""
exp_mask_filter.py — experiment: reject candidates on masked (no-data) regions.
Masked pixels in the warps are 0. A real detection needs real data at BOTH
the before-position (negative, in img1) and after-position (positive, in img2).
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

def fraction_masked(image, x, y, box=5):
    """Fraction of pixels that are exactly 0 (masked) in a box around (x,y)."""
    h, w = image.shape
    x0, x1 = max(0, int(x)-box), min(w, int(x)+box)
    y0, y1 = max(0, int(y)-box), min(h, int(y)+box)
    region = image[y0:y1, x0:x1]
    if region.size == 0:
        return 1.0
    return np.mean(region == 0)

def main():
    img1 = load_fits(sys.argv[1])
    img2 = load_fits(sys.argv[2])
    aligned2 = aa.register(img2, img1)[0]
    diff = (aligned2 - np.median(aligned2[aligned2!=0])) - (img1 - np.median(img1[img1!=0]))

    pos = detect_sources(diff)
    neg = detect_sources(-diff)
    cands = find_candidates(pos, neg)
    print(f"Candidates before mask filter: {len(cands)}")

    MASK_LIMIT = 0.5   # reject if more than 50% of the box is masked

    kept = []
    for n, c in enumerate(cands):
        # positive location checked in image 2 (aligned), negative in image 1
        m_pos = fraction_masked(aligned2, c['pos'][0], c['pos'][1])
        m_neg = fraction_masked(img1, c['neg'][0], c['neg'][1])
        masked = (m_pos > MASK_LIMIT) or (m_neg > MASK_LIMIT)
        status = "REJECT (masked region)" if masked else "KEEP"
        print(f"  #{n+1} pos({c['pos'][0]:.0f},{c['pos'][1]:.0f}) neg({c['neg'][0]:.0f},{c['neg'][1]:.0f}) "
              f"| masked frac pos={m_pos:.2f} neg={m_neg:.2f} | {status}")
        if not masked:
            kept.append(c)

    print(f"\nCandidates after mask filter: {len(kept)}")

if __name__ == "__main__":
    main()