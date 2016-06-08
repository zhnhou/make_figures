import numpy as np
import matplotlib.pyplot as plt
from scipy.io.idl import readsav
from mpl_toolkits.axes_grid1 import make_axes_locatable

from hpylib.util.remote_data import *
from hpylib.mapping.sptsz_map import *
from hpylib.util.sptsz_end2end import *

class create_map_figure(object):
    def __init__(self, fits_file):
        self.fits_file = fits_file
        m = read_sptsz_fits(fits_file)

        self.ra0dec0 = [m['ra0'], m['dec0']]
        self.nside   = [m['nsidex'], m['nsidey']]
        self.map_data = np.flipud(m['map_data'])
        self.reso_arcmin = m['reso_arcmin']

        self.setup_coord()
        
    def setup_coord(self):
        reso_deg = self.reso_arcmin / 60.00

        self.xra = (self.ra0dec0[0] - 0.5*self.nside[0]*reso_deg, self.ra0dec0[0] + 0.5*self.nside[0]*reso_deg)
        self.yra = (self.ra0dec0[1] - 0.5*self.nside[1]*reso_deg, self.ra0dec0[1] + 0.5*self.nside[1]*reso_deg)

    def cut_map(self, xra, yra, replace=True):
        reso_deg = self.reso_arcmin / 60.00

        xarr = np.arange(0,self.nside[0])*reso_deg + self.ra0dec0[0] - 0.5*self.nside[0]*reso_deg
        yarr = np.arange(0,self.nside[1])*reso_deg + self.ra0dec0[1] - 0.5*self.nside[1]*reso_deg

        ipx = np.where((xarr >= min(xra)) & (xarr <= max(xra)))
        ipy = np.where((yarr >= min(yra)) & (yarr <= max(yra)))

        map2d_cut = self.map_data[ipy[0].min():(ipy[0].max()+1),ipx[0].min():(ipx[0].max()+1)]

        if (replace):
            self.map_data = map2d_cut
            self.xra = xra
            self.yra = yra

        return map2d_cut

class create_residual_figure(object):

    def __init__(self, end1, end2, dim1=1, dim2=1, rescale1=1.00, rescale2=1.00, res_beam_cov=None):

        self.end1 = end1
        self.end2 = end2

        self.dim1 = dim1
        self.dim2 = dim2

        self.res_beam_cov = res_beam_cov

        if (not np.array_equal(end1['bands'], end2['bands']) ):
            print "The band definition of two end2end files are different"
            exit()

    def process_end(self):
        
        end1_sims_mean = self.end1['dbs_sims'][:,self.dim1,:].mean(axis=0)
        end2_sims_mean = self.end2['dbs_sims'][:,self.dim2,:].mean(axis=0)

        winfunc_corr = end1_sims_mean - end2_sims_mean

        res_end1_end2 = self.end1['dbs_data'][self.dim1,:]*self.rescale1 - self.end2['dbs_data'][self.dim2,:]*self.rescale2 - winfunc_corr

        res_sims = self.end1['dbs_sims'][:,self.dim1,:] - self.end2['dbs_sims'][:,self.dim2,:]
        cov = np.cov(res_sims.transpose())

        if (not self.res_beam_cov):
            cov += self.res_beam_cov

        d = {'res_data':res_end1_end2, 'res_cov':cov}

        self.res_info = d

    def make_residual_figure(self):

        fig, ax = plt.subplot()

        diag = self.res_info['res_cov']

        ax.set_position([0.1,0.1,0.85,0.85])
#ax.errorbar(self.end1['bands'], self.res_info['res_data'], yerr=, fmt='o', markersize='0', elinewidth=1., capsize=1., capthick=1., label=None)


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
