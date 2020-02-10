import pandas as pd
import numpy as np
from astropy.io import fits
from scipy.stats.distributions import chi2
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70,Om0=.3)

filename = 'uds_20to23kmag_zabove1.fits'

with fits.open(filename) as hdu:
    data = hdu[1].data

    mask = data['stellar_mass_50'] > 10
    data = data[mask]

    mask = data['redshift_50'] > 2
    data = data[mask]

    mask = (chi2.sf(data['chisq_phot'],data['n_bands'])) > 3e-7
    data = data[mask]

    ###########doesn't seem to filter anything out - what's up with that?
    age = cosmo.age(data['redshift_50']).value
    factor = 0.2/age
    mask = 10**data['sSFR_50'] < factor
    data = data[mask]

    print(data.shape)
    newhdu = fits.BinTableHDU(data=data)
    newhdu.writeto('sample.fits')
