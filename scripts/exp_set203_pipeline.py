# exp_set203_pipeline.py - 4-image asteroid detection on IASC set203
import glob
import numpy as np
from astropy.io import fits
import astroalign as aa
from astropy.stats import sigma_clipped_stats
from photutils.detection import DAOStarFinder
from scipy.spatial import KDTree
import matplotlib.pyplot as plt

PIXEL_SCALE = 0.256  # arcsec/px

# --- load all 4 frames ---
files = sorted(glob.glob('data/set203/*.fits'))
print(f"Loading {len(files)} frames...")

images = []
mjds = []
for f in files:
    with fits.open(f) as hdul:
        data = hdul[0].data.astype(float)
        mjd = hdul[0].header['MJD-OBS']
        images.append(data)
        mjds.append(mjd)
        print(f"  {f.split('/')[-1][:20]}...  MJD={mjd:.6f}")

gaps = [(mjds[i+1] - mjds[i]) * 24 * 60 for i in range(3)]
print(f"\nTime gaps: {gaps[0]:.1f} min, {gaps[1]:.1f} min, {gaps[2]:.1f} min")

# --- align frames 2,3,4 to frame 1 ---
print("\nAligning frames to frame 1...")
aligned = [images[0]]
for i in range(1, 4):
    try:
        registered, _ = aa.register(images[i], images[0])
        aligned.append(registered)
        print(f"  Frame {i+1} aligned OK")
    except Exception as e:
        print(f"  Frame {i+1} FAILED: {e}")

print(f"\nAligned {len(aligned)} of 4 frames successfully.")

# --- detect sources in each aligned frame ---
mask0 = (aligned[0] == 0)
_, ref_median, ref_std = sigma_clipped_stats(aligned[0], sigma=3.0, mask=mask0)
THRESHOLD = 4.2 * ref_std
print(f"\nReference std={ref_std:.1f}, global threshold={THRESHOLD:.1f} counts")
# --- find saturated/very bright regions (bright stars) to avoid ---
# use a fixed multiple of the noise above the median, ignoring zero-padding
valid_pixels = aligned[0][aligned[0] > 0]
img_median = np.median(valid_pixels)
SATURATION_LEVEL = img_median + 200 * ref_std  # very bright = bright star
bright_mask = aligned[0] > SATURATION_LEVEL
bright_ys, bright_xs = np.where(bright_mask)
bright_pixels = np.column_stack([bright_xs, bright_ys])
if len(bright_pixels) > 0:
    bright_tree = KDTree(bright_pixels)
    print(f"Found {len(bright_pixels)} saturated pixels (level>{SATURATION_LEVEL:.0f})")
else:
    bright_tree = None
    print("No saturated pixels found")

BRIGHT_REJECT_RADIUS = 200  # reject candidates within this many px of a bright star
print("\nDetecting sources in each frame...")
all_sources = []
for i, img in enumerate(aligned):
    mask = (img == 0)
    _, median, _ = sigma_clipped_stats(img, sigma=3.0, mask=mask)
    dao = DAOStarFinder(fwhm=5.0, threshold=THRESHOLD)
    sources = dao(img - median, mask=mask)
    if sources is None:
        sources = []
        print(f"  Frame {i+1}: 0 sources")
    else:
        print(f"  Frame {i+1}: {len(sources)} sources")
    all_sources.append(sources)

# --- setup positions and trees for all frames ---
MIN_MOVE = 3
MAX_MOVE = 100
CONFIRM_RADIUS = 5

def classify(rate):
    if rate < 1:
        return "too slow -> likely star or artifact"
    elif rate < 50:
        return "slow -> TNO-like"
    elif rate < 500:
        return "main-belt asteroid range"
    else:
        return "fast -> NEO-like"

pos1 = np.array([[s['x_centroid'], s['y_centroid']] for s in all_sources[0]])
pos2 = np.array([[s['x_centroid'], s['y_centroid']] for s in all_sources[1]])
pos3 = np.array([[s['x_centroid'], s['y_centroid']] for s in all_sources[2]])
pos4 = np.array([[s['x_centroid'], s['y_centroid']] for s in all_sources[3]])

tree1 = KDTree(pos1)
tree2 = KDTree(pos2)
tree3 = KDTree(pos3)
tree4 = KDTree(pos4)

# --- find candidates between frame 1 and frame 2 (mutual nearest-neighbor) ---
print("\nFinding candidates (frame 1 -> frame 2)...")

dists_fwd, idxs_fwd = tree2.query(pos1)
dists_bwd, idxs_bwd = tree1.query(pos2)

candidates = []
for i, (dist, j) in enumerate(zip(dists_fwd, idxs_fwd)):
    if MIN_MOVE < dist < MAX_MOVE:
        if idxs_bwd[j] == i:
            candidates.append({'pos1': pos1[i], 'pos2': pos2[j], 'dist': dist})

print(f"  {len(candidates)} candidate(s) between frame 1 and 2")

# --- also find candidates from frame 2->3 and frame 3->4 ---
print("\nFinding additional candidates from frame 2->3 and frame 3->4...")

def find_mutual_candidates(posA, posB, gap_idx):
    treeA = KDTree(posA)
    treeB = KDTree(posB)
    dists_fwd, idxs_fwd = treeB.query(posA)
    dists_bwd, idxs_bwd = treeA.query(posB)
    cands = []
    for i, (dist, j) in enumerate(zip(dists_fwd, idxs_fwd)):
        if MIN_MOVE < dist < MAX_MOVE:
            if idxs_bwd[j] == i:
                cands.append({'pos1': posA[i], 'pos2': posB[j],
                              'dist': dist, 'gap_idx': gap_idx})
    return cands

cands_23 = find_mutual_candidates(pos2, pos3, 1)
cands_34 = find_mutual_candidates(pos3, pos4, 2)
print(f"  Frame 2->3: {len(cands_23)} candidates")
print(f"  Frame 3->4: {len(cands_34)} candidates")

# --- unified line-check across all starting pairs ---
print("\nRunning unified line-check across all frame pairs...")
confirmed = []

def line_check_and_add(cands, start_gap):
    for c in cands:
        A = c['pos1']
        B = c['pos2']
        step = B - A
        frame_positions = [A - step*start_gap + step*i for i in range(4)]
        hits = 0
        frame_dists = []
        for fi, (tree, pred) in enumerate(zip([tree1, tree2, tree3, tree4], frame_positions)):
            d, _ = tree.query(pred)
            frame_dists.append(d)
            if d < CONFIRM_RADIUS:
                hits += 1
        if hits >= 3:
            rate = c['dist'] * PIXEL_SCALE / (gaps[start_gap] / (24*60))
            frame1_pos = frame_positions[0]   # where this track would be in frame 1

            # reject if near a saturated bright star
            if bright_tree is not None:
                bright_dist, _ = bright_tree.query([A[0], A[1]])
                if bright_dist < BRIGHT_REJECT_RADIUS:
                    continue  # too close to a bright star, skip
            is_dup = any(
                np.hypot(frame1_pos[0]-e['frame1_pos'][0],
                         frame1_pos[1]-e['frame1_pos'][1]) < 30
                for e in confirmed
            )
            if not is_dup:
                confirmed.append({**c, 'rate': rate, 'hits': hits,
                                 'frame_dists': frame_dists,
                                 'frame1_pos': frame1_pos})
                near_obj1 = np.hypot(A[0]-664, A[1]-2077) < 50
                flag = " *** OBJ0001/GE56 ***" if near_obj1 else ""
                print(f"  CONFIRMED #{len(confirmed)}: {c['dist']:.1f}px = "
                      f"{rate:.0f} arcsec/day "
                      f"({A[0]:.0f},{A[1]:.0f})->({B[0]:.0f},{B[1]:.0f}) "
                      f"| {hits}/4 frames{flag}")

line_check_and_add(candidates, 0)
line_check_and_add(cands_23, 1)
line_check_and_add(cands_34, 2)

print(f"\n{len(confirmed)} total confirmed candidate(s)")
print("\n=== FINAL CONFIRMED CANDIDATES ===")
for n, c in enumerate(confirmed):
    rate = c['dist'] * PIXEL_SCALE / (gaps[0] / (24*60))
    A, B = c['pos1'], c['pos2']
    print(f"\nCandidate #{n+1}")
    print(f"  Speed:     {rate:.0f} arcsec/day  ->  {classify(rate)}")
    print(f"  Track:     ({A[0]:.0f},{A[1]:.0f}) -> ({B[0]:.0f},{B[1]:.0f})")
    print(f"  Frame dists: {[f'{d:.1f}' for d in c['frame_dists']]}")
    print(f"  Confirmed: {c['hits']}/4 frames")


# --- visualize all 4 frames with confirmed candidates ---
print("\nGenerating visualization...")
fig, axes = plt.subplots(1, 4, figsize=(20, 6))
fig.suptitle(f'Confirmed moving objects: {len(confirmed)}', fontsize=14)

colors = ['cyan', 'orange', 'red']
for frame_idx, (ax, img) in enumerate(zip(axes, aligned)):
    mean, med, std = sigma_clipped_stats(img, sigma=3.0)
    ax.imshow(img, cmap='gray', vmin=med-2*std, vmax=med+4*std, origin='upper')
    ax.set_title(f'Frame {frame_idx+1}')
    for n, c in enumerate(confirmed):
        A, B = c['pos1'], c['pos2']
        step = B - A
        color = colors[n % len(colors)]
        if frame_idx == 0:
            pos = A
        elif frame_idx == 1:
            pos = B
        elif frame_idx == 2:
            pos = B + step
        else:
            pos = B + 2 * step
        circle = plt.Circle((pos[0], pos[1]), 15,
                            color=color, fill=False, linewidth=2)
        ax.add_patch(circle)
        ax.text(pos[0]+18, pos[1]+18, f'#{n+1}',
               color=color, fontsize=8, fontweight='bold')

plt.tight_layout()
plt.savefig('outputs/set203_tracks.png', dpi=150, bbox_inches='tight')
print("Saved to outputs/set203_tracks.png")
plt.show()

# --- zoom in on candidate #3 ---
if len(confirmed) >= 3:
    fig2, axes2 = plt.subplots(1, 4, figsize=(20, 5))
    fig2.suptitle('Candidate #3 zoom (NEO-like) 4/4 frames', fontsize=13)
    c = confirmed[2]
    A = c['pos1']
    B = c['pos2']
    step = B - A
    positions = [A, B, B+step, B+2*step]
    ZOOM = 80
    for frame_idx, (ax, img, pos) in enumerate(zip(axes2, aligned, positions)):
        mean, med, std = sigma_clipped_stats(img, sigma=3.0)
        x, y = int(pos[0]), int(pos[1])
        x1, x2 = max(0, x-ZOOM), min(img.shape[1], x+ZOOM)
        y1, y2 = max(0, y-ZOOM), min(img.shape[0], y+ZOOM)
        cutout = img[y1:y2, x1:x2]
        ax.imshow(cutout, cmap='gray', vmin=med-2*std, vmax=med+4*std, origin='upper')
        ax.set_title(f'Frame {frame_idx+1}  ({x},{y})')
        circle = plt.Circle((pos[0]-x1, pos[1]-y1), 10,
                            color='red', fill=False, linewidth=2)
        ax.add_patch(circle)
    plt.tight_layout()
    plt.savefig('outputs/set203_candidate3_zoom.png', dpi=150, bbox_inches='tight')
    print("Saved to outputs/set203_candidate3_zoom.png")
    plt.show()

# --- convert candidate positions to RA/Dec ---
from astropy.wcs import WCS

print("\nConverting candidate positions to sky coordinates...")

with fits.open(files[0]) as hdul:
    hdr = hdul[0].header

w = WCS(naxis=2)
w.wcs.crpix = [hdr['CRPIX1'], hdr['CRPIX2']]
w.wcs.cdelt = [hdr['CDELT1'], hdr['CDELT2']]
w.wcs.crval = [hdr['CRVAL1'], hdr['CRVAL2']]
w.wcs.ctype = [hdr['CTYPE1'], hdr['CTYPE2']]
crota = np.radians(hdr['CROTA2'])
w.wcs.pc = [[np.cos(crota), -np.sin(crota)],
            [np.sin(crota), np.cos(crota)]]

for n, c in enumerate(confirmed):
    A = c['pos1']
    ra, dec = w.all_pix2world(A[0], A[1], 0)
    print(f"  Candidate #{n+1}: pixel ({A[0]:.0f},{A[1]:.0f}) -> RA={float(ra):.5f}, Dec={float(dec):.5f}")

# --- query SkyBoT ---
from astroquery.imcce import Skybot
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u

print("\nQuerying SkyBoT (via astroquery) for known objects...")
epoch = Time(mjds[0], format='mjd')

for n, c in enumerate(confirmed):
    A = c['pos1']
    ra, dec = w.all_pix2world(A[0], A[1], 0)
    ra, dec = float(ra), float(dec)
    field = SkyCoord(ra=ra, dec=dec, unit='deg')
    print(f"\n--- Candidate #{n+1} (RA={ra:.5f}, Dec={dec:.5f}) ---")
    try:
        results = Skybot.cone_search(field, 0.1 * u.deg, epoch)
        print(results)
    except Exception as e:
        print(f"  No known objects found, or query failed: {e}")

# --- find ALL known asteroids in field ---
print("\nQuerying SkyBoT for ALL known objects in the field...")
center_ra, center_dec = w.all_pix2world(aligned[0].shape[1]/2, aligned[0].shape[0]/2, 0)
center_ra, center_dec = float(center_ra), float(center_dec)
field_center = SkyCoord(ra=center_ra, dec=center_dec, unit='deg')

def safe_float(x):
    try:
        return float(x.value)
    except AttributeError:
        return float(x)

try:
    field_results = Skybot.cone_search(field_center, 0.2 * u.deg, epoch)
    print(f"\nFound {len(field_results)} known object(s) in field:\n")
    for row in field_results:
        try:
            ra_obj = safe_float(row['RA'])
            dec_obj = safe_float(row['DEC'])
            name = row['Name']
            v_mag = safe_float(row['V'])
            px, py = w.all_world2pix(ra_obj, dec_obj, 0)
            px, py = safe_float(px), safe_float(py)
            in_frame = (0 <= px <= aligned[0].shape[1]) and (0 <= py <= aligned[0].shape[0])
            flag = "IN FRAME" if in_frame else "outside frame"
            print(f"  {name}: RA={ra_obj:.5f}, Dec={dec_obj:.5f}, "
                  f"V={v_mag:.2f}, pixel=({px:.0f},{py:.0f})  [{flag}]")
        except Exception as row_e:
            print(f"  (skipped a row: {row_e})")
except Exception as e:
    print(f"  Query failed: {e}")

# --- diagnose missed objects ---
print("\n\nDiagnosing missed known objects...")
missed = {
    "2002 GE56": (171.27461, 3.26445),
    "2021 RY128": (171.37278, 3.30067),
    "2015 RM287": (171.37743, 3.22112),
    "2002 QK157": (171.41576, 3.25610),
}
tree_f1 = KDTree(pos1)
for name, (ra_m, dec_m) in missed.items():
    px, py = w.all_world2pix(ra_m, dec_m, 0)
    px, py = float(px), float(py)
    dist, idx = tree_f1.query([px, py])
    print(f"\n{name}: predicted pixel ({px:.0f},{py:.0f})")
    if dist < 15:
        print(f"  -> DETECTED in frame 1, {dist:.1f}px from a real source.")
    else:
        print(f"  -> NOT detected in frame 1 (nearest source is {dist:.1f}px away).")