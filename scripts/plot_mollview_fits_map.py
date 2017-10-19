import numpy as np
import math
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as py
from matplotlib.colors import LogNorm
import healpy as hp

noise = hp.read_map('../output/nside-1/shot_noise_zmin_0.010_zmax_0.015_Nside_0001.fits')

hp.mollview(noise,coord='G')

py.savefig('../output/nside-1/shot_noise_zmin_0.010_zmax_0.015_Nside_0001.pdf')

py.close()

map = hp.read_map('../output/nside-1/map_zmin_0.010_zmax_0.015_Nside_0001.fits')

hp.mollview(map,coord='G')

py.savefig('../output/nside-1/map_zmin_0.010_zmax_0.015_Nside_0001.pdf')

exit()
