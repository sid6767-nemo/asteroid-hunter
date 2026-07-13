# asteroid-hunter

Finds moving asteroids in telescope images automatically — no clicking through frames by hand.

Give it a set of images of the same patch of sky taken minutes apart, and it lines them up, finds every point of light, picks out the ones that move in a straight line (asteroids move, stars don't), throws out false alarms, and checks each one against a database of known asteroids.

## Results

Run on the IASC **set203** practice dataset (4 frames, ~21 minutes apart), the pipeline confirmed **3 real named asteroids**, each matched to the official catalog within ~3.8 arcseconds:

| Object | Speed | Brightness (V) |
|--------|-------|----------------|
| 2015 RM287 | 975 ″/day | 21.75 |
| 2004 RH62 | 745 ″/day | 19.90 |
| 2002 GE56 | 818 ″/day | 19.16 |

![detected asteroids](outputs/set203_tracks.png)

## How it works

1. **Align** — line up the frames so the stars sit in the same place (using `astroalign`).
2. **Detect** — find every point of light in each frame (`photutils` segmentation + deblending).
3. **Track** — stars stay put, asteroids move. The pipeline looks for points that shift in a straight line across the frames.
4. **Filter** — throw out false alarms: junk near bright stars, and anything that isn't a clean straight-line mover with steady brightness across all 4 frames.
5. **Cross-match** — check each surviving candidate against the SkyBoT database to see if it's a known asteroid.

## How to run it

```bash
# clone the repo
git clone https://github.com/sid6767-nemo/asteroid-hunter.git
cd asteroid-hunter

# set up a virtual environment
python3 -m venv .venv
source .venv/bin/activate      # on Windows: .venv\Scripts\activate

# install dependencies
pip install -r requirements.txt

# run the pipeline on the included sample data
python scripts/exp_set203_pipeline.py
```

The sample dataset (`data/set203/`) is included, so it works straight out of the box. Results are saved to the `outputs/` folder.

## Known limitations

This is a work in progress. Right now:

- The two faintest known asteroids in the field (V≈22.5 and 23.6) are below the detection threshold — too dim to see in single frames without stacking.
- One catchable object (2007 DT63, V≈20.1) isn't being picked up yet — it's on the to-do list.

**Goal:** detect all 7 known asteroids in the set203 field cleanly and automatically, with no hardcoded positions.

## Data credit

The sample images come from the **International Astronomical Search Collaboration (IASC)** and the **Pan-STARRS** survey. They're used here for educational and practice purposes.

## Faint recovery: hearing below the detection limit

My pipeline detects asteroids by finding them in single frames first, then
linking the ones that move. That means it can never see anything too faint
for one exposure. The "hunt by ear" score works differently — it pools light
from all 4 frames along a guessed track (shift-and-stack) before deciding.
In theory that reaches fainter objects. I tested whether it actually does.

**The experiment** (`scripts/exp_faint_recovery.py`): ask SkyBoT for every
known asteroid in the set203 field, convert each one's predicted sky position
to a pixel track, and run the blind roam score there — plus the same search
on 12 patches of empty sky, as a "fake floor" for what pure noise can fluke.

**Result:**

| asteroid    | V    | pipeline | roam score | recovered velocity vs truth       |
|-------------|------|----------|-----------|-----------------------------------|
| 2002 GE56   | 19.2 | found    | 25.3      | 820"/d @ 90° vs 788 @ 90° — match |
| 2004 RH62   | 19.9 | found    | 17.6      | 745"/d @ 80° vs 725 @ 82° — match |
| 2015 RM287  | 21.8 | found    | 1.9       | 995"/d @ 80° vs 947 @ 80° — match |
| 2002 QK157  | 22.5 | *missed* | 0.46      | 845"/d @ 68° vs 826 @ 71° — match |
| 2021 RY128  | 23.6 | *missed* | 0.23      | 125° off — no match (noise)       |

Empty-sky floor (same search effort): 0.30.

**The recovery.** At the predicted position of 2002 QK157 — V22.5, below my
pipeline's single-frame limit — the blind search returned its strongest
response at 845"/day, 68°, matching the asteroid's true motion to within 3°
and 2% in speed. The score amplitude alone (1.5× the empty-sky floor) would
not be conclusive; the velocity agreement is what makes it a recovery. Blank
sky can fluke a score, but it flukes at a random velocity out of ~8,640
searched — not at the exact motion of a cataloged asteroid.

**The honest half.** 2021 RY128 (V23.6) also crept slightly above the floor,
but its recovered velocity was 125° off the truth — noise, and the method
correctly rejects it. And the scores track brightness exactly as physics
predicts (V19.2 → 25 down to V23.6 → noise), which is how I know the score
is measuring something real.

What I learned: my earlier "false positives" question had two answers at
once. One faint signal was a real asteroid my pipeline couldn't see; another
was noise. The velocity-match test is what separates them.