import numpy as np
import pandas as pd
from loaddata import loaduds
import bagpipes as pipes

def magnitudecalc(f):
    mag = 23.9 - 2.5*np.log10(f)
    return mag

filt_list = np.loadtxt('filters_uds/uds_filt_list.txt', dtype='str')

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

########### calculate IDs for fitting here
IDs = np.genfromtxt('masscut_indices.txt')
print(len(IDs))
fit_cat = pipes.fit_catalogue(IDs, fit_instructions, loaduds,spectrum_exists=False, cat_filt_list=filt_list, run="uds_20to23kmag", make_plots =False, full_catalogue = True)
#fit_cat.fit(verbose=False, n_live=1000, mpi_serial = True)
