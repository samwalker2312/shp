import numpy as np
from astropy.io import fits
import pandas as pd

##################################################
#for photometric data

def loaduds(ID):
    dataset = pd.read_csv('K_selected_uds_photoz_masses_photom_upload.cat', header=0, delim_whitespace=True, index_col=0)
    #IDcol = np.genfromtxt('K_selected_uds_photoz_masses_photom_upload.cat', skip_header=2, usecols=0)

    #extracts ra and dec columns
    #coords = dataset.loc[ID,'RA':'Dec']

    #extracts values with z > 1 (presented as an example)
    #data = dataset[dataset['zmed'] > 1.]
    #data = data[data.index.isin(ID)]

    ID = int(ID)

    isofactor = dataset.loc[ID,'isofactor']

    fluxes_nonisophot = isofactor*(dataset.loc[ID, ['U_2as', 'B_2as', 'V_2as', 'R_2as', 'i_2as', 'z_2as','znew_2as', 'Y_2as', 'J_2as', 'H_2as', 'K_2as']].to_numpy())
    fluxerrs_nonisophot = isofactor*(dataset.loc[ID, ['U_2as_err', 'B_2as_err', 'V_2as_err', 'R_2as_err', 'i_2as_err', 'z_2as_err', 'znew_2as_err', 'Y_2as_err', 'J_2as_err', 'H_2as_err', 'K_2as_err']].to_numpy())
    fluxes_irac = dataset.loc[ID, ['ch1_flux', 'ch2_flux']].to_numpy()
    fluxerrs_irac = dataset.loc[ID, ['ch1_err', 'ch2_err']].to_numpy()
    fluxes = np.concatenate((fluxes_nonisophot, fluxes_irac))
    fluxerrs = np.concatenate((fluxerrs_nonisophot, fluxerrs_irac))
    fluxes[fluxes<-98] = 0
    fluxerrs[fluxerrs<-98] = 9.9*10**99.
    output = np.c_[fluxes, fluxerrs]
    for i in range(len(output)):
	if i < len(fluxes_nonisophot):
		max_snr = 20.
	else:
		max_snr = 10.
	if output[i,0]/output[i,1]>max_snr:
		output[i,1] = output[i,0]/max_snr

    print(output)
    return output

def load_goodss(ID):
    """ Load UltraVISTA photometry from catalogue. """

    # load up the relevant columns from the catalogue.
    cat = np.loadtxt("hlsp_candels_hst_wfc3_goodss-tot-multiband_f160w_v1-1photom_cat.txt",
                     usecols=(10, 13, 16, 19, 22, 25, 28, 31, 34, 37, 40, 43, 46, 49, 52, 55,
                              11, 14, 17, 20, 23, 26, 29, 32, 35, 38, 41, 44, 47, 50, 53, 56))

    # Find the correct row for the object we want.
    row = int(ID) - 1
    # Extract the object we want from the catalogue.
    fluxes = cat[row, :16]
    fluxerrs = cat[row, 16:]

    hdul = fits.open('vandels_catalogues/VANDELS_CDFS_HST_PHOT_v1.0.fits')
    data = hdul[1].data
    IDfield = data.field('ID')
    IDlist = IDfield.tolist()
    if ID in IDlist:
        mask = (IDfield == ID)
        fluxes[11] = data.field('Ks_HAWKI_flux')[mask][0]
        fluxerrs[11] = data.field('Ks_HAWKI_err')[mask][0]
    hdul.close()

    # Turn these into a 2D array.
    photometry = np.c_[fluxes, fluxerrs]

    # blow up the errors associated with any missing fluxes.
    for i in range(len(photometry)):
        if (photometry[i, 0] == 0.) or (photometry[i, 1] <= 0):
            photometry[i,:] = [0., 9.9*10**99.]

    print(photometry)

    # Enforce a maximum SNR of 20, or 10 in the IRAC channels.
    for i in range(len(photometry)):
        if i < 12:
            max_snr = 20.

        else:
            max_snr = 10.

        if photometry[i, 0]/photometry[i, 1] > max_snr:
            photometry[i, 1] = photometry[i, 0]/max_snr

    return photometry

def load_goodss_old(ID):
    """ Load UltraVISTA photometry from catalogue. """

    # load up the relevant columns from the catalogue.
    cat = np.loadtxt("hlsp_candels_hst_wfc3_goodss-tot-multiband_f160w_v1-1photom_cat.txt",
                     usecols=(10, 13, 16, 19, 25, 28, 31, 34, 37, 43, 46, 49, 52, 55,
                              11, 14, 17, 20, 26, 29, 32, 35, 38, 44, 47, 50, 53, 56))

    # Find the correct row for the object we want.
    row = int(ID) - 1
    # Extract the object we want from the catalogue.
    fluxes = cat[row, :14]
    fluxerrs = cat[row, 14:]

    # Turn these into a 2D array.
    photometry = np.c_[fluxes, fluxerrs]

    # blow up the errors associated with any missing fluxes.
    for i in range(len(photometry)):
        if (photometry[i, 0] == 0.) or (photometry[i, 1] <= 0):
            photometry[i,:] = [0., 9.9*10**99.]

    # Enforce a maximum SNR of 20, or 10 in the IRAC channels.
    for i in range(len(photometry)):
        if i < 10:
            max_snr = 20.

        else:
            max_snr = 10.

        if photometry[i, 0]/photometry[i, 1] > max_snr:
            photometry[i, 1] = photometry[i, 0]/max_snr

    return photometry


####################################################################
#for spectroscopic data

def bin(spectrum, binn):
    """ Bins up two or three column spectral data by a specified factor. """

    binn = int(binn)
    nbins = len(spectrum)/binn
    binspec = np.zeros((nbins, spectrum.shape[1]))

    for i in range(binspec.shape[0]):
        spec_slice = spectrum[i*binn:(i+1)*binn, :]
        binspec[i, 0] = np.mean(spec_slice[:, 0])
        binspec[i, 1] = np.mean(spec_slice[:, 1])

        if spectrum.shape[1] == 3:
            binspec[i,2] = (1./float(binn)
                            *np.sqrt(np.sum(spec_slice[:, 2]**2)))

    return binspec


def load_vandels_spec(ID):
    """ Loads VANDELS spectroscopic data from file. """

    hdulist = fits.open("VANDELS_CDFS_" + ID + ".fits")

    spectrum = np.c_[hdulist[1].data["WAVE"][0],
                     hdulist[1].data["FLUX"][0],
                     hdulist[1].data["ERR"][0]]

    mask = (spectrum[:,0] < 9250.) & (spectrum[:,0] > 5250.)

    return bin(spectrum[mask], 2)
