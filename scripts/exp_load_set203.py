# exp_load_set203.py - inspect the IASC 4-image set
import glob
from astropy.io import fits

files = sorted(glob.glob('data/set203/*.fits'))
print(f"Found {len(files)} files:\n")

for f in files:
    with fits.open(f) as hdul:
        hdr = hdul[0].header
        data = hdul[0].data
        print(f)
        print(f"   shape:    {data.shape}")
        print(f"   DATE-OBS: {hdr.get('DATE-OBS')}")
        print(f"   MJD-OBS:  {hdr.get('MJD-OBS')}")
        print(f"   EXPTIME:  {hdr.get('EXPTIME')}")
        print()

print("=== Full header of first file ===")
with fits.open(files[0]) as hdul:
    print(repr(hdul[0].header))