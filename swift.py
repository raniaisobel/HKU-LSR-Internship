import numpy as np
import matplotlib.pyplot as plt
from scipy import optimize
from scipy.fft import fft, fftfreq
import matplotlib.patches as mpatches
from astropy.io import fits
import astropy.units as u
from astropy.coordinates import SkyCoord
import io

# PNe Hash Database
dtype = [('idPNMain', 'i4'),('PNG', 'U10'),('Name', 'U10'),('PNstat', 'U1'),('RAJ2000', 'U12'),('DECJ2000', 'U12'),('DRAJ2000', 'f8'),('DDECJ2000', 'f8'),('Glon', 'f8'),('Glat', 'f8'),('MajDiam', 'f8'),('mainClass', 'U1'),('subClass', 'U1')]
data2 = np.genfromtxt('all_true_all_PNe.csv', delimiter=',',dtype=dtype,names=True)

go2 = data2['Glon']
ga2 = data2['Glat']
glon2 = go2 * u.deg
glat2 = ga2 * u.deg

g_coords2 = SkyCoord(l=glon2,b=glat2, frame='galactic')

# Swift Database
dtypeee = [('blank', 'f8'),('source_number', 'i8'), ('name', 'U25'), ('ra', 'U15'), ('dec', 'U15'), ('lii', 'f8'), ('bii', 'f8'), ('error_radius', 'f8'), ('count_rate_fb', 'f8'), ('count_rate_fb_pos_err', 'f8'), ('count_rate_fb_neg_err', 'f8'), ('pow_flux_a', 'f8'), ('pow_flux_a_pos_err', 'f8'), ('pow_flux_a_neg_err', 'f8'), ('apec_flux_a', 'f8'), ('apec_flux_a_pos_err', 'f8'), ('apec_flux_a_neg_err', 'f8')]
data5 = np.genfromtxt('swift.txt', delimiter='|', skip_header=3, dtype=dtypeee, invalid_raise=False)

go5 = data5['lii']
ga5 = data5['bii']

glon5 = go5 * u.deg
glat5 = ga5 * u.deg

g_coords5 = SkyCoord(l=glon5,b=glat5, frame='galactic')

# Matrix
print(len(g_coords2)) # - 3345
print(len(g_coords5)) # = 206,335
# so its (would be) a matrix of 3345 x 206,335

"""
ds = np.zeros((len(g_coords2), len(g_coords5)))
ds = g_coords2[:, None].separation(g_coords5[None, :]).deg
print(ds)
"""

# Testing chunk
