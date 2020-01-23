from loaddata import load_goodss
from loaddata import load_vandels_spec
import bagpipes as pipes
import numpy as np

filt_list = np.loadtxt("filters_guo/goodss_filt_list.txt", dtype="str")

dblplaw = {}
dblplaw["tau"] = (0., 15.)
dblplaw['tau_prior'] = 'uniform'

dblplaw["alpha"] = (0.1, 1000.)
dblplaw["alpha_prior"] = "log_10"

dblplaw["beta"] = (0.1, 1000.)
dblplaw["beta_prior"] = "log_10"

dblplaw["massformed"] = (1., 13.)
dblplaw["massformed_prior"] = "log_10"

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
fit_instructions["t_bc"] = 0.01
fit_instructions["nebular"] = nebular

def analysisfunc(fit):
	ID = str(fit.galaxy.ID)
	fit.posterior.get_advanced_quantities()

	if "redshift" in fit.fitted_model.params:
		redshift = np.median(fit.posterior.samples["redshift"])
	else:
		redshift = fit.fitted_model.model_components["redshift"]

	# Plot the posterior photometry and full spectrum.
	wavs = (fit.posterior.model_galaxy.wavelengths*(1.+redshift))*10

	posterior = fit.posterior.samples["spectrum_full"]
	spectrum = np.median(posterior, axis=0)
	data = np.c_[wavs, spectrum]
	fname = 'G13_spectra/' + ID + '_spectrum.dat'
	np.savetxt(fname, data)



IDs = np.array([2408,1424,2782,18180,8785,22211,22085,3973,9209])
fit_cat = pipes.fit_catalogue(IDs, fit_instructions, load_goodss,\
 spectrum_exists=False, cat_filt_list=filt_list, run="guo_2020atopten", make_plots = True, analysis_function=analysisfunc, full_catalogue=True)
fit_cat.fit(verbose=False)
