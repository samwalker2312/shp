import pandas as pd
import numpy as np
from astropy.io import fits
from scipy.stats.distributions import chi2
from astropy.cosmology import FlatLambdaCDM
import math

cosmo = FlatLambdaCDM(H0=70,Om0=.3)

filename = 'output_cats/mass_select.fits'

with fits.open(filename) as hdu:
    data = hdu[1].data

    mask = data['stellar_mass_50'] > 10
    data = data[mask]

    mask = data['redshift_50'] > 2
    data = data[mask]

    mask = (chi2.sf(data['chisq_phot'],(data['n_bands']-10))) > (0.5*math.erfc(5*(2**-0.5)))
    data = data[mask]

    age = cosmo.age(data['redshift_50']).value*(10**9)
    factor = 0.2/age
    mask = 10**data['sSFR_50'] < factor
    data = data[mask]

    print(data.shape)
    newhdu = fits.BinTableHDU(data=data)
    newhdu.writeto('output_cats/sample.fits')
