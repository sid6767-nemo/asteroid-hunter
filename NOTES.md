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