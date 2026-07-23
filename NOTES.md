# Asteroid GPU Project Notes

## Days 1–6
First few days were just the basics — setting up folders and files for this project on my laptop. The actual GPU work comes later in the project and I'll use Google Colab's free T4 access rather than buying hardware.

I made the repo public on GitHub from day one — figured if I'm going to do this, I want it visible.

By the time the first week was over I had aligned a few images with astroalign. I learned how the stars are used as reference for the alignment. I don't really understand what a FITS file is, although I know it stores an image with a header.

## Day 7
Worked on subtraction and learned about the colour differences between the g and r filters.

## Day 8
Improved subtraction using warps. Failed at first because of a lot of differences between the warps, but after selecting good warps it worked better.

## Days 9–10
Added detection using photutils DAOStarFinder. After running it I got 205 positives (parts where warp2 is brighter) and 11 negatives (parts where warp1 is brighter) — most were just stars. After matching them I got 35 candidates, but each negative (mostly artifacts — saturated stars, edge effects) matched to many nearby positive stars.

## Day 11
Fixed it by joining mutually nearest positives and negatives together, between 3 and 40 pixels apart.

## Day 12
Combined all my separate scripts into one pipeline — detect.py. Now runs with one command on any two FITS files. Raised MIN_MOVE from 3 to 8 px to kill the alignment errors.

## Day 13
Added a filter that rejects candidates sitting near bright stars. These spots threw up false detections because stars don't subtract cleanly.

## Day 14
Added a filter for no-data regions, because the warps have areas where pixels are just 0 and a candidate sitting on one isn't real.

I also understood that my pipeline only detects sideways (transverse) motion — an object moving straight toward or away from the telescope wouldn't change position between the two images, so my method can't catch that kind of motion. A real limitation to keep in mind.

## Day 15 : converting pixel motion to arcsec/day
I converted how far each candidate moved from pixels into how fast it
moved across the sky, in arcseconds per day. That let me actually
measure whether something was moving like an asteroid. All 3 of my
candidates came out too slow, so they're most likely stars or
artifacts, not asteroids.

## Day 16 : flag for review
I added a flag-for-review step. It flags any candidate moving fast
enough to maybe be a real solar-system object, and keeps the slow
ones as low-priority instead of deleting them. On my data all 3 were
too slow, so 0 of 3 got flagged, which is the correct result since
there were no real moving objects in this field.

## Day 17: 4-image line-check

Switched to a practice dataset containing four FITS images taken 21 minutes apart,
rather than two. The pipeline found 96 candidates between the first two frames. To
filter those down, I predicted where each candidate should appear in frames 3 and 4
based on its speed and direction. Anything that didn't show up near the predicted
position in at least 2 frames was rejected. 3 candidates survived (739 arcsec/day
NEO-like, 153 arcsec/day main-belt, 70 arcsec/day main-belt).

## Day 18: visualization and SkyBot cross-match

Built a four-frame visualization showing each candidate's confirmed positions
across all four images, plus a zoomed-in view of candidate #3. Converted each
candidate's position into real sky coordinates and checked them against
SkyBot, a database of known asteroids, to see if each was new or already
found. Candidate #3 matched a real asteroid within 4 arcseconds, strong proof
the pipeline works. Candidates #1 and #2 had no match, probably meaning they
were artifacts rather than real objects.

## Day 19: de-blending

Candidates #1 and #2 were false asteroid detections, so using photutils, 
I tried de-blending the code to remove any noise that the nearby stars might 
have caused.

## Day 20-21: figuring out the issue

GE56(near candidate #1) was never a deblending problem. The code was detecting it fine,
the old bright-star filter was deleting it for being near a star. That old filter rejected 
everything within 200px of any saturated pixel → killed 23% of the whole image. So I only 
rejected near the few giant stars (the ones with spikes), zone size scaled to star size. 
Small stars left alone. I got GE56 but 38 extra candidates, I tweaked certain values and
got it to 3  actual candidates. The goal get all 7 asteroid to be detected.

The objects the pipeline recovers are real main-belt asteroids orbiting between Mars and Jupiter. 
Example: 2002 GE56, ~6 km across, V=19.16 — recovered cleanly despite sitting next 
to a saturated star.

## Day 22 — Flask upload page
I built a Flask upload page that runs the real pipeline on uploaded frames and 
shows the result image and candidates.
I wanted a more user-friendly option in case I turn it into a website.
I tested it on a fresh practice set the pipeline
had never seen, as a proof-check that it works on new data. I ran into a problem though: 
some faint real asteroids and false positives were hard to tell apart. That left me stuck.

## Day 23 — Solving false positives with stack-along-track
I figured out a stack-along-track concentration test to separate faint asteroids from 
false positives. I stacked cutouts centered on the tracked position across all the frames. 
A real object piles up into a concentrated point, 
while a fake smears out. I validated it on a few more datasets and got
almost all of them right, four out of five correct.
The one limitation left was false positives from stars sitting
near the detector's grid lines, which drift slightly between frames.
I also couldn't fully check for missed asteroids yet,
since my Astrometrica overlay wasn't working. That's it for now.

## Day 24 — Orbit determination, Stage 1
I started the orbit goal in stages. Stage 1 makes the pipeline
output each asteroid's real sky position (RA/Dec) at each frame's
exact time, and also in the standard MPC format that orbit tools read.
Validated on set203: the positions match the catalog within ~4 arcsec.
I learned that a single night's short arc can't determine a real
orbit, so I use known orbits from JPL Horizons for the visualization.

## Day 25 — 3D orbit viewer in the web app
I built a 3D orbit viewer with Three.js showing the asteroid's real
orbit around the Sun with Earth, animated over time with speed and
zoom controls and a rich star and nebula background.
I wired it into the Flask app: the pipeline now identifies detected
asteroids via SkyBoT, and each named one gets a View orbit button
that fetches its real orbit from JPL Horizons.
There is also a View all orbits together button that shows every
detected asteroid's orbit on one shared map.

## Day 26 - Sonification for accessibility
I built a way to find asteroids by ear, for blind and people with poor vision.
You move a cursor over the sky and the sound gets louder and higher as
you near a real asteroid, and peaks when you land on it.
I tried a asteroid detection version first but it drowned faint asteroids
in star noise, so I built the sound-field from the real detections,
which is reliable.
It is embedded below the blink viewer in the web app.
Two things to fix next: the readings are uneven from different sides,
and I want it to speak the asteroid's name and coordinates out loud.

## Days 27–29 — SkyBoT matching fixed (twice), standalone demos, README rewrite
The SkyBoT name matching had been silently reading the wrong columns, so matches were coming out wrong or missing. Fixed it to actually read the RA/Dec columns properly. Also centered the sonar beacon on the frame-1 position so the sound field lines up with the backdrop image instead of drifting off it.

Split off standalone versions of the detect-by-ear tool and orbit viewer that don't need the full pipeline running, so I can demo either one on its own. Cleaned up the repo so generated outputs and textures aren't tracked in git anymore.

The SkyBoT fix from Day 27 wasn't fully right either — had to go back in and fix the name matching a second time. Also added per-frame positions for hotspots so I can see where a candidate actually sat in every frame, not just the summary track.

Rewrote the README properly: web app, orbit viewer, hazard assessment, sonification, and an honest limitations section instead of just pipeline notes. It was still describing the project as a command-line pipeline when it's a web app now.

## Days 30–32 — hunt-by-ear score engine, QK157 recovery
Built the actual shift-and-stack score engine behind hunt-by-ear and ported it to JS so it runs live in the browser, wired into new `/hunt` routes. Before this the sonification was just a demo; now it's scoring a real stacked signal as you move.

Ran an experiment asking whether the ear score could recover known asteroids too faint for the pipeline to see in a single frame. Pulled every known asteroid in the set203 field from SkyBoT, turned each one's sky position into a pixel track, and scored it two ways: once knowing its real velocity, once blind (guessing, like a user would). Compared against empty sky and against the 3 pipeline-confirmed asteroids as a reference.

2002 QK157 (V=22.5, well below what the pipeline can see in one frame) came back with a high score and a velocity match to the catalog within 3° direction and 2% speed — a real recovery, not a lucky noise spike.

## Days 33–34 — 4/7 in the README, hunt-by-ear hardened
Updated the README results to 4/7 known asteroids recovered (3 pipeline + QK157 via shift-and-stack), and put hunt-by-ear on the front page instead of burying it as a side feature. Added real web app screenshots with captions so the README shows what it actually looks like instead of just describing it.

Turned the QK157 score of 0.45 into an actual threshold: below 0.35 is noise floor (refused), 0.35–0.45 is flagged low-confidence, above 0.45 is a candidate. Added bright-star rejection zones and cursor-anchored zoom so the tool isn't blindsided by saturated stars, and made hunt sessions per-browser instead of shared. 

Rebuilt the upload flow so importing frames goes straight into sonification instead of stopping at the results page first. Also added per-frame candidate confirmation — hit space to mark a candidate in each frame — feeding into a new review page that checks marks against the catalog and reads back a spoken summary.

## Days 35–36 — click-safe confirming mode, auto-blink review
The per-frame confirmation flow from Day 34 broke easily — a stray mouse move could wipe out a seeded cursor position, Tab could jump focus out of the flow, and the finish button had no guard against double-submits. Fixed all three, added predicted-position cursor seeding so the cursor starts near where the object should be in each frame, and the review blink viewer now shows my per-frame clicks as green circles so I can see what I actually marked.

Built an auto-blink comparator into the review page so frames cycle on their own instead of me stepping through by hand, with every candidate drawn in its own color. Matched candidates get a direct link to their orbit viewer; unmatched single-night candidates get an honest note that there's no orbit to show, instead of silently omitting anything. Fixed a display bug where the review canvas stretched across the full-width wrapper instead of hugging the actual image — circles were landing in the wrong spot relative to what was on screen.