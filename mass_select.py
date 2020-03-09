import pandas as pd
import numpy as np
from astropy.io import fits
from scipy.stats.distributions import chi2
from astropy.cosmology import FlatLambdaCDM

cosmo = FlatLambdaCDM(H0=70,Om0=.3)

filename = 'output_cats/uds_20to23kmag.fits'

with fits.open(filename) as hdu:
    data = hdu[1].data

    mask = data['stellar_mass_50'] > 10
    data = data[mask]

    mask = data['redshift_50'] > 2
    data = data[mask]

    print(data.shape)
    index = (data['#ID'].astype(str))
    index = (index.view(np.ndarray)).astype(int)
    np.savetxt('masscut_indices.txt', index)