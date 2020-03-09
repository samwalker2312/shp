from astropy.io import fits
import numpy as np

redshifts = np.genfromtxt('redshift_percentiles_full').astype(float)
ssfrs = np.genfromtxt('ssfr_percentiles_full').astype(float)
zcol = fits.Column(name='redshift_2.5', format = 'D', array = redshifts)
ssfrcol = fits.Column(name = 'ssfr_97.5', format = 'D', array = ssfrs)
hdu_new = fits.BinTableHDU.from_columns([zcol, ssfrcol])

with fits.open('output_cats/sample.fits') as hdu_old:
    new_columns = hdu_old[1].columns
    new_columns += hdu_new.columns
    sample_hdu = fits.BinTableHDU.from_columns(new_columns)
    sample_hdu.writeto('output_cats/sample_fulldata.fits')
