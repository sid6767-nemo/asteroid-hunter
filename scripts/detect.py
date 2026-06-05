"""
detect.py — Asteroid Hunter pipeline
Takes two FITS images of the same field, aligns them, subtracts,
detects sources, and finds candidate moving objects.

Usage:
    python scripts/detect.py data/warp1.fits data/warp2.fits
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
import astroalign as aa

# --- Settings ---
MIN_MOVE = 8      # min pixel motion (raised from 3 to beat alignment jitter)
MAX_MOVE = 40     # max pixel motion
DETECT_FWHM = 3.0
DETECT_SIGMA = 5.0


def load_fits(path):
    with fits.open(path) as hdul:
        data = hdul[0].data.astype(np.float64)
        header = hdul[0].header
    return np.nan_to_num(data), header


def align(img2, img1):
    """Align img2 onto img1's pixel grid."""
    aligned, footprint = aa.register(img2, img1)
    return aligned


def normalize(img):
    """Subtract median, return centered image."""
    valid = img[img != 0]
    return img - np.median(valid)


def detect_sources(image):
    """Detect point sources brighter than threshold."""
    mean, median, std = sigma_clipped_stats(image, sigma=3.0)
    finder = DAOStarFinder(fwhm=DETECT_FWHM, threshold=DETECT_SIGMA * std)
    sources = finder(image - median)
    if sources is None:
        return np.empty((0, 2))
    return np.array([sources['xcentroid'], sources['ycentroid']]).T


def nearest(point, candidates):
    if len(candidates) == 0:
        return None, np.inf
    dists = np.hypot(candidates[:, 0] - point[0], candidates[:, 1] - point[1])
    idx = np.argmin(dists)
    return idx, dists[idx]

def local_brightness(image, x, y, box=10):
    """Max pixel value in a box around (x,y) — high means near a bright/saturated star."""
    h, w = image.shape
    x0, x1 = max(0, int(x)-box), min(w, int(x)+box)
    y0, y1 = max(0, int(y)-box), min(h, int(y)+box)
    region = image[y0:y1, x0:x1]
    if region.size == 0:
        return 0
    return np.max(region)

def reject_near_saturation(candidates, img1, img2):
    """Remove candidates sitting near saturated/bright regions."""
    sat_level = np.percentile(img1[img1 > 0], 99.5)
    kept = []
    for c in candidates:
        x, y = c['pos']
        b1 = local_brightness(img1, x, y)
        b2 = local_brightness(img2, x, y)
        if b1 <= sat_level and b2 <= sat_level:
            kept.append(c)
    return kept, sat_level    
def find_candidates(pos_xy, neg_xy):
    """Mutual-nearest-neighbour matching within motion window."""
    confirmed = []
    for i, p in enumerate(pos_xy):
        j, d_pn = nearest(p, neg_xy)
        if j is None or not (MIN_MOVE <= d_pn <= MAX_MOVE):
            continue
        k, _ = nearest(neg_xy[j], pos_xy)
        if k == i:
            confirmed.append({'pos': p, 'neg': neg_xy[j], 'dist': d_pn})
    return confirmed


def main():
    if len(sys.argv) != 3:
        print("Usage: python scripts/detect.py <image1.fits> <image2.fits>")
        sys.exit(1)

    path1, path2 = sys.argv[1], sys.argv[2]
    print(f"Loading {path1} and {path2}...")
    img1, hdr1 = load_fits(path1)
    img2, hdr2 = load_fits(path2)

    print("Aligning...")
    aligned2 = align(img2, img1)

    print("Subtracting...")
    diff = normalize(aligned2) - normalize(img1)

    print("Detecting sources...")
    pos_xy = detect_sources(diff)
    neg_xy = detect_sources(-diff)
    print(f"  {len(pos_xy)} positive, {len(neg_xy)} negative")

    print("Finding candidates...")
    candidates = find_candidates(pos_xy, neg_xy)
    print("Filtering candidates near saturation...")
    before = len(candidates)
    candidates, sat_level = reject_near_saturation(candidates, img1, aligned2)
    print(f"  Saturation cutoff: {sat_level:.0f} | {before} -> {len(candidates)} after filter")
    print(f"\n{len(candidates)} candidate moving object(s):")
    for n, c in enumerate(candidates):
        print(f"  #{n+1}: moved {c['dist']:.1f} px, "
              f"({c['neg'][0]:.0f},{c['neg'][1]:.0f}) -> "
              f"({c['pos'][0]:.0f},{c['pos'][1]:.0f})")

    # Visualize
    fig, ax = plt.subplots(figsize=(10, 10))
    dr = np.percentile(np.abs(diff[diff != 0]), 99)
    ax.imshow(diff, cmap='gray', origin='lower', vmin=-dr, vmax=dr)
    for c in candidates:
        ax.annotate('', xy=c['pos'], xytext=c['neg'],
                    arrowprops=dict(arrowstyle='->', color='lime', lw=2))
        ax.scatter(*c['neg'], s=120, edgecolor='red', facecolor='none', lw=1.8)
        ax.scatter(*c['pos'], s=120, edgecolor='cyan', facecolor='none', lw=1.8)
    ax.set_title(f"Candidate moving objects: {len(candidates)}")
    plt.tight_layout()
    plt.savefig('outputs/pipeline_result.png', dpi=100)
    plt.show()
    print("\nSaved to outputs/pipeline_result.png")


if __name__ == "__main__":
    main()