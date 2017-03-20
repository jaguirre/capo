import numpy as np
import healpy as hp
from astropy import units as u
from astropy import constants as c

# Longest baseline
b_max = 1000*u.m
# Highest frequency
nu_max = 200.*u.MHz
# Lowest frequency
nu_min = 100.*u.MHz
# Highest resolution
res_max = np.degrees(1./(b_max * nu_max/c.c).to(u.dimensionless_unscaled).value)*60.
# Oversampling of synthesizbeam
pixel_oversamp = 2.

# Calculate FoV
FoV = np.degrees(((c.c/nu_min)/(14.*u.m)).to(u.dimensionless_unscaled).value)
FoV_oversamp = 2.

# Calculate number of healpix pixels
nside = 1024
npix_hp = hp.nside2npix(nside)
res = hp.nside2resol(nside,arcmin=True)

area_arcmin2 = np.power(res_max/pixel_oversamp,2)

# Calculate number of pixels in 12 hour equatorial map
sq_arcmin_eq = (FoV*60.*FoV_oversamp) * (180.*60.)
npix_eq = sq_arcmin_eq/area_arcmin2

# Calcualte number of pixels in 1 hour snapshot
sq_arcmin_snap = np.power(FoV*60.*FoV_oversamp,2)
npix_snap = sq_arcmin_snap/area_arcmin2

# Number of channels
nchan = 200.

# Conversion to size
bytes_per_float = 8
mb_per_vox = bytes_per_float / 1e6

size_eq = npix_eq * nchan * mb_per_vox
size_snap = npix_snap * nchan * mb_per_vox
