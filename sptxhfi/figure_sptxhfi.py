import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
from sys import exit

from hpylib.util.remote_data import *

def restore_save(savfile):

    n = readsav(savfile)
    key = n.keys()
    if (len(key) != 1):
        exit(".sav file is not a combined end2end file")

    num_bands = int(n[key[0]]['num_bands'][0])
    bands = n[key[0]]['bands'][0]
    dbs_data = n[key[0]]['dbs_data_combined'][0]
    dbs_sims = n[key[0]]['dbs_sims_combined'][0] # (nsims, nspecs, nbands)

    winminell = int(n[key[0]]['winminell'][0])
    winmaxell = int(n[key[0]]['winmaxell'][0])

    winfunc_data = n[key[0]]['winfunc_data_combined'][0]
    winfunc_sims = n[key[0]]['winfunc_sims_combined'][0]

    cov_sv    = n[key[0]]['cov_sv_combined'][0]
    cov_noise = n[key[0]]['cov_noise_combined'][0]

    d = {'num_bands':num_bands, 'bands':bands, 
         'dbs_data':dbs_data, 'dbs_sims':dbs_sims,
         'winminell':winminell, 'winmaxell':winmaxell,
         'winfunc_data':winfunc_data, 'winfunc_sims':winfunc_sims,
         'cov_sv':cov_sv, 'cov_noise':cov_noise}

    return d


def plot_spt150hfi143_bandpower(pdf_file=None):
    
    spt150xspt150_file = '~/data_midway/projects/sptxhfi/pspec/bandpower_spt_sn_spt_sn/end_combined_spt150sn_spt150sn.sav'
    spt150xhfi143_file = '~/data_midway/projects/sptxhfi/pspec/bandpower_spt_sn_hfi_sn/end_combined_spt150sn_hfi143sn.sav'
    #spt150xhfi217_file = '~/data_midway/projects/sptxhfi/pspec/bandpower_spt_sn_hfi_sn/end_combined_spt150sn_hfi217sn.sav'
    hfi143xhfi143_file = '~/data_midway/projects/sptxhfi/pspec/bandpower_hfi_sn_hfi_sn/end_combined_hfi143sn_hfi143sn.sav'

    s150     = restore_save( sync_from_remote('midway', spt150xspt150_file) )
    s150h143 = restore_save( sync_from_remote('midway', spt150xhfi143_file) )
    #s150h217 = restore_save( sync_from_remote('midway', spt150xhfi217_file) )
    h143     = restore_save( sync_from_remote('midway', hfi143xhfi143_file) )

    dls_150_file     = '~/data_midway/projects/sptxhfi/simulations/input/dls_input_spt_150.txt'
    dls_220_file     = '~/data_midway/projects/sptxhfi/simulations/input/dls_input_spt_220.txt'
    dls_150x220_file = '~/data_midway/projects/sptxhfi/simulations/input/dls_ave_150x220.txt'

    dls_theory_150     = np.loadtxt(sync_from_remote('midway', dls_150_file), usecols=[1])
    dls_theory_220     = np.loadtxt(sync_from_remote('midway', dls_220_file), usecols=[1])
    #dls_theory_150x220 = np.loadtxt(sync_from_remote('midway', dls_150x220_file), usecols=[1])

    dbs_ave_150     = np.mean(s150['dbs_sims'][:,1,:], axis=0)
    dbs_ave_143     = np.mean(h143['dbs_sims'][:,1,:], axis=0)
    dbs_ave_150x143 = np.mean(s150h143['dbs_sims'][:,1,:], axis=0)
    #dbs_ave_150x217 = np.mean(s150h217['dbs_sims'][:,1,:], axis=0)

    dbs_err_150     = np.sqrt(np.diag(s150['cov_sv'][1,:,1,:]))
    dbs_err_143     = np.sqrt(np.diag(h143['cov_sv'][1,:,1,:]))
    dbs_err_150x143 = np.sqrt(np.diag(s150h143['cov_sv'][1,:,1,:]))
    #dbs_err_150x217 = np.sqrt(np.diag(s150h217['cov_sv'][1,:,1,:]))


    rescale = 1.0100

    dbs_data_150     = s150['dbs_data'][1,:] * rescale**2
    dbs_data_150x143 = (s150h143['dbs_data'][1,:] - (dbs_ave_150x143 - dbs_ave_150)) * rescale
    dbs_data_143     = h143['dbs_data'][1,:]     - (dbs_ave_143 - dbs_ave_150)
    #dbs_data_150x217 = (s150h217['dbs_data'][1,:] - (dbs_ave_150x217 - dbs_ave_150)) * rescale


    fig, ax = plt.subplots()
    ax.set_position([0.1,0.1,0.85,0.85])

    ax.errorbar(s150['bands'], dbs_data_150, yerr=dbs_err_150, fmt='o', markersize='0', elinewidth=1., capsize=1., capthick=1., label=r'$\mathrm{SPT^{150}_{half1}\;\times\;SPT^{150}_{half2}}$')
    ax.errorbar(s150h143['bands']-10, dbs_data_150x143, yerr=dbs_err_150x143, fmt='o', markersize='0', elinewidth=1., capsize=1., capthick=1., label=r'$\mathrm{SPT^{150}_{full}\;\times\;HFI^{143}_{full}}$')
    ax.errorbar(h143['bands']+10, dbs_data_143, yerr=dbs_err_143, fmt='o', markersize='0', elinewidth=1., capsize=1., capthick=1., label=r'$\mathrm{HFI^{143}_{half1}\;\times\;HFI^{143}_{half2}}$')
    
    ax.legend()

    if pdf_file is None:
        pdf_file = 'bandpower.pdf'

    plt.xlim([625,3000])
    plt.ylim([20,3000])
    plt.yscale('log')
    plt.xlabel(r'$\ell$', fontsize=16)
    plt.ylabel(r'$\mathcal{D}_{\ell}\ [\mathrm{\mu K^2}]$', fontsize=16)
    plt.savefig(pdf_file, format='pdf')


def plot_spt150hfi217_bandpower(pdf_file=None):
    
    spt150xspt150_file = '~/data_midway/projects/sptxhfi/pspec/bandpower_spt_sn_spt_sn/end_combined_spt150sn_spt150sn.sav'
    spt150xhfi217_file = '~/data_midway/projects/sptxhfi/pspec/bandpower_spt_sn_hfi_sn/end_combined_spt150sn_hfi217sn.sav'
    hfi217xhfi217_file = '~/data_midway/projects/sptxhfi/pspec/bandpower_hfi_sn_hfi_sn/end_combined_hfi217sn_hfi217sn.sav'

    s150     = restore_save( sync_from_remote('midway', spt150xspt150_file) )
    s150h217 = restore_save( sync_from_remote('midway', spt150xhfi217_file) )
    h217     = restore_save( sync_from_remote('midway', hfi217xhfi217_file) )

    dls_150_file     = '~/data_midway/projects/sptxhfi/simulations/input/dls_input_spt_150.txt'
    dls_220_file     = '~/data_midway/projects/sptxhfi/simulations/input/dls_input_spt_220.txt'
    dls_150x220_file = '~/data_midway/projects/sptxhfi/simulations/input/dls_ave_150x220.txt'

    dls_theory_150     = np.loadtxt(sync_from_remote('midway', dls_150_file), usecols=[1])
    dls_theory_220     = np.loadtxt(sync_from_remote('midway', dls_220_file), usecols=[1])
    dls_theory_150x220 = np.loadtxt(sync_from_remote('midway', dls_150x220_file), usecols=[1])

    dbs_ave_150     = np.mean(s150['dbs_sims'][:,1,:], axis=0)
    dbs_ave_217     = np.mean(h217['dbs_sims'][:,1,:], axis=0)
    dbs_ave_150x217 = np.mean(s150h217['dbs_sims'][:,1,:], axis=0)

    dbs_err_150     = np.sqrt(np.diag(s150['cov_sv'][1,:,1,:]))
    dbs_err_217     = np.sqrt(np.diag(h217['cov_sv'][1,:,1,:]))
    dbs_err_150x217 = np.sqrt(np.diag(s150h217['cov_sv'][1,:,1,:]))


    rescale = 1.0100

    dbs_data_150     = s150['dbs_data'][1,:] * rescale**2
    dbs_data_217     = h217['dbs_data'][1,:]     - (dbs_ave_217 - dbs_ave_150)
    dbs_data_150x217 = (s150h217['dbs_data'][1,:] - (dbs_ave_150x217 - dbs_ave_150)) * rescale


    fig, ax = plt.subplots()
    ax.set_position([0.1,0.1,0.85,0.85])

    ax.errorbar(s150['bands'], dbs_data_150, yerr=dbs_err_150, fmt='o', markersize='0', elinewidth=1., capsize=1., capthick=1., label=r'$\mathrm{SPT^{150}_{half1}\;\times\;SPT^{150}_{half2}}$')
    ax.errorbar(s150h217['bands']-10, dbs_data_150x217, yerr=dbs_err_150x217, fmt='o', markersize='0', elinewidth=1., capsize=1., capthick=1., label=r'$\mathrm{SPT^{150}_{full}\;\times\;HFI^{217}_{full}}$')
    ax.errorbar(h217['bands']+10, dbs_data_217, yerr=dbs_err_217, fmt='o', markersize='0', elinewidth=1., capsize=1., capthick=1., label=r'$\mathrm{HFI^{217}_{half1}\;\times\;HFI^{217}_{half2}}$')
    
    ax.legend()

    if pdf_file is None:
        pdf_file = 'bandpower.pdf'

    plt.xlim([625,3000])
    plt.ylim([20,3000])
    plt.yscale('log')
    plt.xlabel(r'$\ell$', fontsize=16)
    plt.ylabel(r'$\mathcal{D}_{\ell}\ [\mathrm{\mu K^2}]$', fontsize=16)
    plt.savefig(pdf_file, format='pdf')
