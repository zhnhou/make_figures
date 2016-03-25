from hpylib.util.remote_data import *
import cPickle as pickle

dl_base_file = '~/data_midway/projects/cmb_pspec/baseline/dl_lensed_base_TT_lowP_lensing.pkl'
dl_lensed = pickle.load(open(sync_from_remote('midway',dl_base_file,update=False),'rb'))
