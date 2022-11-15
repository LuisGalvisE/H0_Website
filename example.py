from H0calculate import *
from astropy import constants as const
import matplotlib.pyplot as plt
from time import time

start = time ()

cval = const.c.to('km/s').value

skymap = 'GW170817_skymap.fits.gz'

# sky position of galaxy
ra_deg_NGC4993 = 360/24*(13+9/60+47/3660)
dec_deg_NGC4993 = -(23 + 23/60 + 4/3600)
ra_rad_NGC4993 = ra_deg_NGC4993*np.pi/180
dec_rad_NGC4993 = dec_deg_NGC4993*np.pi/180

#redshift distribution (normal distribution: mu, sigma)
sig_cz_NGC4993 = 166
mu_cz_NGC4993 = 3017


H0_calculate = H0calculate (skymap, mu_cz_NGC4993/cval, sig_cz_NGC4993/cval, ra_rad_NGC4993, dec_rad_NGC4993, H0max=200)

H0_array, pH0 = H0_calculate.probH0 ()

plt.plot ( H0_array, pH0)

plt.xlim (H0_array [0], H0_array [-1])

plt.xlabel (r'$H_{0}$', size=15)
plt.ylabel (r'$p(H_{0})$', size=15)

plt.tight_layout ()

plt.savefig ('pH0_GW170817_test')

print ('Time:', time()-start, 'sec')


