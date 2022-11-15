'''
Script for fast calculation of H0 posterior for Bright GW event 

Tathagata ghosh
'''

import numpy as np
import healpy as hp
from scipy.stats import norm, truncnorm
from scipy.integrate import simpson
from scipy.interpolate import interp1d
from astropy.cosmology import FlatLambdaCDM
from astropy import constants as const, units as u
from ligo.skymap.io.fits import read_sky_map


class H0calculate :

	def __init__ (self, skymap, counterpart_z, counterpart_sigmaz, counterpart_ra, counterpart_dec, Om0=0.3, H0min=20, H0max=140, H0bins=100, zcut=None) :

		# redshift-luminosity distance relation
		cosmoref = FlatLambdaCDM(H0=70, Om0=Om0, Tcmb0=2.725)
		if zcut is None :
			zcut = 10
		zref = np.arange (0,zcut+0.01,0.01)
		dlref = cosmoref.luminosity_distance (zref).value
		self.dlH02z = interp1d (dlref*70, zref)
		diff_comoving_vol = interp1d(zref, cosmoref.differential_comoving_volume(zref).to(u.Gpc**3 / u.sr).value)
		
		# skymap
		(self.prob, self.distmu, self.distsigma, self.distnorm), metadata = read_sky_map (skymap, distances=True, moc=False, nest=True)

		# pixel corresponding to sky position of identified galaxy		
		npix = len(self.prob)
		nside = hp.npix2nside (npix)
		self.counterpart_pix = hp.ang2pix (nside, np.pi/2 - counterpart_dec, counterpart_ra, nest=True)
		
		# minimum and maximum distance of GW event
		mu = self.distmu[(self.distmu < np.inf) & (self.distmu > 0)]
		dlGWmin = 0.5*min(mu)
		dlGWmax = 2*max(mu)
		self.dlGW_array = np.linspace (dlGWmin, dlGWmax,1000)
		
		# redshift samples
		a = (0.0 - counterpart_z) / counterpart_sigmaz
		zsmear = truncnorm.rvs (a,5, loc=counterpart_z, scale=counterpart_sigmaz, size=10000)
		self.counterpart_z_sample = zsmear[np.where(zsmear<zcut)[0]].flatten()
		
		# H0 array
		self.H0bins = H0bins
		self.H0_array = np.linspace ( H0min, H0max, self.H0bins)
		self.Om0 = Om0
		
		# redshift prior
		self.dVc_by_dz = diff_comoving_vol (self.counterpart_z_sample)


	def __RedshiftEvolutionPowerLaw (self, z, power) :
		
		return (1+z)**(power-1)
	
	def probH0 (self) :
	
		distmu_los = self.distmu [self.counterpart_pix]
		distsigma_los = self.distsigma [self.counterpart_pix]
		likelihood_x_dl_skymap = norm (distmu_los, distsigma_los)
		
		# redshift prior
		pz = self.__RedshiftEvolutionPowerLaw (self.counterpart_z_sample, 0)*self.dVc_by_dz
		
		pH0 = np.zeros (self.H0bins)
		
		for hh, H0 in enumerate (self.H0_array) :
		
			cosmo = FlatLambdaCDM(H0=H0, Om0=self.Om0, Tcmb0=2.725)
			
			zmin = self.dlH02z(self.dlGW_array [0]*H0) *0.5
			zmax = self.dlH02z(self.dlGW_array [-1]*H0) *2
			redshift_bins = 10000
			zGW_array_temp = np.linspace (zmin,zmax,redshift_bins)
			
			dl_array_temp = cosmo.luminosity_distance (zGW_array_temp).value

			likelihood_x_z_H0= likelihood_x_dl_skymap.pdf(dl_array_temp)
			likelihood_x_z_H0 /= simpson (likelihood_x_z_H0, zGW_array_temp)

			px_z_H0_interp = interp1d(zGW_array_temp,likelihood_x_z_H0,kind='linear',bounds_error=False,fill_value=0)
			px_z_H0 = px_z_H0_interp(self.counterpart_z_sample)
		
			pH0 [hh] = np.sum(px_z_H0*pz)/H0**3
		
			
		pH0 /= simpson (pH0, self.H0_array)
		
		return self.H0_array, pH0



