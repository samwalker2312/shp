from loaddata import loaduds
#from loaddata import load_vandels_spec
import bagpipes as pipes
import numpy as np
#from astropy.cosmology import FlatLambdaCDM
import os
import paramiko
from astropy.io import fits

filt_list = np.loadtxt("filters_uds/uds_filt_list.txt", dtype="str")

dblplaw = {}
dblplaw["tau"] = (0.1, 15.)
dblplaw['tau_prior'] = 'uniform'

dblplaw["alpha"] = (0.01, 1000.)
dblplaw["alpha_prior"] = "log_10"

dblplaw["beta"] = (0.01, 1000.)
dblplaw["beta_prior"] = "log_10"

dblplaw["massformed"] = (0., 13.)
dblplaw["massformed_prior"] = "uniform"

dblplaw["metallicity"] = (0.1, 2.5)
dblplaw["metallicity_prior"] = "log_10"

dust = {}
dust["type"] = "Salim"
dust["Av"] = (0.,8.)
dust['Av_prior'] = 'uniform'

dust['delta'] = (-.3,.3)
dust['delta_prior'] = 'Gaussian'
dust['delta_prior_mu'] = 0.
dust['delta_prior_sigma'] = .1

dust['B'] = (0,5)
dust['B_prior'] = 'uniform'

dust["eta"] = 2.

nebular = {}
nebular["logU"] = -3.

fit_instructions = {}
fit_instructions["redshift"] = (0., 10.)
fit_instructions["redshift_prior"] = 'uniform'
fit_instructions["dblplaw"] = dblplaw
fit_instructions["dust"] = dust
fit_instructions["nebular"] = nebular
fit_instructions['t_bc'] = 0.01

#cosmo = FlatLambdaCDM(H0=70, Om0 = .3, Tcmb0=2.725)
def analysisfunc(fit):
	ID = str(fit.galaxy.ID)
	fit.posterior.get_advanced_quantities()

	if "redshift" in fit.fitted_model.params:
		redshift = np.median(fit.posterior.samples["redshift"])
	else:
		redshift = fit.fitted_model.model_components["redshift"]

	# Plot the posterior photometry and full spectrum.
	wavs = (fit.posterior.model_galaxy.wavelengths)*(1+redshift)/10

	posterior = fit.posterior.samples["spectrum_full"]
	spectrum = np.median(posterior, axis=0)
	data = np.c_[wavs, spectrum]

	fname = 'uds_spectra/' + ID + '_spectrum.dat'
	np.savetxt(fname, data)

filename = 'output_cats/subsample.fits'

with fits.open(filename) as hdu:
    data = hdu[1].data

mask = data['redshift_50'] > 3.5
data = data[mask]
IDs = (data['#ID'].astype(str))
print(IDs)

IDs = (IDs.view(np.ndarray)).astype(int)
fit_cat = pipes.fit_catalogue(IDs, fit_instructions, loaduds,\
 spectrum_exists=False, cat_filt_list=filt_list, run="uds_20to23kmag", make_plots = True, analysis_function=analysisfunc, full_catalogue=True)
fit_cat.fit(verbose=False, n_live = 1000)
