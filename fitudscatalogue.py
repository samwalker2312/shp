import numpy as np
import pandas as pd
from loaddata import loaduds
import bagpipes as pipes

def magnitudecalc(f):
    zpt = 3631
    mag = -2.5*np.log(np.divide(f*(10**-6),zpt))
    return mag

filt_list = np.loadtxt('filters_uds/uds_filt_list.txt', dtype='str')

dblplaw = {}
dblplaw["tau"] = (0., 15.)
dblplaw['tau_prior'] = 'uniform'

dblplaw["alpha"] = (0.1, 1000.)
dblplaw["alpha_prior"] = "log_10"

dblplaw["beta"] = (0.1, 1000.)
dblplaw["beta_prior"] = "log_10"

dblplaw["massformed"] = (0., 13.)
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

########### calculate IDs for fitting here
data = pd.read_csv('K_selected_uds_photoz_masses_photom_upload.cat', header=0, delim_whitespace=True, index_col=0)
k_mag = data['K_iso'].apply(magnitudecalc)
#print(k_mag)
IDs = data.index[(data['zmed'] > 1) & (data['SED_FLAG'] == 0)].tolist()
print(len(IDs))
print(IDs[:1])

#perhaps truncate list using slices to first 10 or so for initial tests
fit_cat = pipes.fit_catalogue(IDs[:1], fit_instructions, loaduds,\
 spectrum_exists=False, cat_filt_list=filt_list, run="uds_testrun", make_plots = True)
fit_cat.fit(verbose=False, n_live=400)
