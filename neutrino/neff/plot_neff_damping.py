import matplotlib.pyplot as plt
import numpy as np
import cPickle as pickle
from hpylib.util.remote_data import *
from hpylib.cosmology.read_planck_data import *

planck_pspec_file = 'data/planck_data/2015/powspec/COM_PowerSpect_CMB_R2.02.fits'
planck_pspec = read_planck_pspec_fits(sync_from_remote('midway',planck_pspec_file))

dl_neff_file = '~/data_midway/projects/neutrino/neff/param_sample/Dl_base_TT_lowP_lensing_const_omegab_zeq_thetas.pkl'
dl_load = pickle.load(open(sync_from_remote('midway',dl_neff_file,update=True),'rb'))

num_sample = dl_load['num_sample']
dl_neff    = dl_load['dl_lensed']

fig, ax = plt.subplots()
ax.set_position([0.1,0.1,0.85,0.85])
