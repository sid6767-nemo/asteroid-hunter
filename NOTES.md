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
