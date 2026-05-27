# Asteroid GPU Project Notes

First few days, was just the basics setting up folders and files for this project on my laptop. The actual GPU work comes later in the project and I'll use Google Colab's free T4 access rather than buying hardware.

I made the repo public on GitHub from day one — figured if I'm going to do this, I want it visible.

By the time the first week was over i had aligned a few images to astroalign. I learned how the stars are used as reference for the alignment. I don't really understand what a FITS file is although i know it stores an image with a header. 

On day 7 I worked on on subtraction and learned about the colour differences between g and r filter. Day 8, I improved subtraction using warps. Failed at first because of a lot of differences in warps, but after selecting good warps it worked better.

For day 9 and 10 I added detection using photutils DAOStarFinder. Affter running it I got 205 positive (basically parts where warp2 is brighter) and 11 negatives (parts where warp1 is brighter), most were stars. After working it out i got 35 candidates but each negative (mostly artifacts — saturated stars, edge effects) matched to many nearby positive stars. On day 11 I fixed my code by joining the mutually nearest positives and negatives together between 3 and 40 pixels.
