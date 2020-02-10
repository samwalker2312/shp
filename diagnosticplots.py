import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import pandas as pd
filename = 'uds_20to23kmag_zabove1.fits'

def magnitudecalc(f):
    mag = 23.9 - 2.5*np.log(f)
    return mag

with fits.open(filename) as hdu:
    fitted_data = hdu[1].data

mask = fitted_data['stellar_mass_50'] > 0
fitted_data = fitted_data[mask]
IDs = np.array(fitted_data['#ID']).astype('float')
print(IDs)

data = pd.read_csv('K_selected_uds_photoz_masses_photom_upload.cat', header=0, delim_whitespace=True, index_col=0)
k_mag = data['K_iso'].apply(magnitudecalc)
fluxes = data.loc[:,['U_2as', 'B_2as', 'V_2as', 'R_2as', 'i_2as', 'znew_2as', 'Y_2as', 'J_2as', 'H_2as', 'K_2as', 'ch1_flux', 'ch2_flux']]
data['filter'] = np.sum(fluxes.values == -99., axis=1)
data = data[(data['zmed'] > 1) & (20 < k_mag) & (k_mag < 23) & (data['H_2as'] > -99) & (data['filter'] < 3)]

data = data.loc[IDs, :]
print(data.shape)
print(fitted_data.shape)

zmed = np.array(data['zmed'])
plt.plot(fitted_data['redshift_50'], zmed, 'o')
plt.xlabel('Fitted median redshift')
plt.ylabel("Derek's redshift")
plt.show()

plt.plot(fitted_data['stellar_mass_50'], data['log10(Mass_iso)'], 'o')
plt.xlabel('Fitted stellar mass')
plt.ylabel("Derek's mass")
plt.show()

plt.plot(np.log10(fitted_data['SFR_50']), data['log10(SFR_iso)'],'o')
plt.xlabel('Fitted log(SFR)')
plt.ylabel("Derek's log(SFR)")
plt.show()

plt.plot(fitted_data['sSFR_50'], data['log10(sSFR_iso)'],'o')
plt.xlabel('Fitted log(sSFR)')
plt.ylabel("Derek's log(sSFR)")
plt.show()

plt.plot(fitted_data['stellar_mass_50'],np.log10(fitted_data['SFR_50']), 'o')
plt.ylabel('Fitted log(SFR)')
plt.xlabel('Fitted log(stellar mass)')
plt.show()

plt.plot(fitted_data['dust:Av_50'],fitted_data['sSFR_50'],'o')
plt.xlabel('V-band dust attenuation')
plt.ylabel('sSFR')
plt.show()
