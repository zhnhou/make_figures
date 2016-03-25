import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import mpl_toolkits.axisartist as AA

import numpy as np
import cPickle as pickle
from hpylib.util.remote_data import *
from hpylib.cosmology.read_planck_data import *
from hpylib.cosmology.read_spt_data import *

pdf_file = 'neff_phase_zoom.pdf'

planck_pspec_file = '~/data_midway/planck_data/2015/powspec/COM_PowerSpect_CMB_R2.02.fits'
plk = read_planck_pspec_fits(sync_from_remote('midway',planck_pspec_file))

spt = read_spt_bandpower('/Users/zhenhou/data_local/data_spt/bandpowers_lps12')

dl_neff_file = '~/data_midway/projects/neutrino/neff/param_sample/Dl_base_TT_lowP_lensing_const_omegab_zeq_thetas_thetad.pkl'
dl_load = pickle.load(open(sync_from_remote('midway',dl_neff_file,update=False),'rb'))

num_sample = dl_load['num_sample']
dl_neff    = dl_load['dl_lensed']
lmax       = dl_load['lmax']

fig, ax = plt.subplots()
ax.set_position([0.1,0.1,0.85,0.65])

plk_ell = plk['highl_TT_bands']
plk_l4cl = plk['highl_TT_dbs']*plk_ell**2
plk_err = plk['highl_TT_err']*plk_ell**2

spt_ell = spt['band']
spt_l4cl = spt['dbs_fgsub']*spt_ell**2
spt_err = spt['err']*spt_ell**2

ell = np.arange(lmax+1)

lines = []
colors = ['red','darkorange','deepskyblue','green','purple']
for i in np.arange(1,6):
    nm = dl_neff[400,0,0] / dl_neff[400,0,i]
    line = plt.plot(ell, dl_neff[:,0,i]*ell**2 * nm, color=colors[i-1], zorder=1, label=r'$\mathrm{N_{eff}}='+str(i)+'$')

err1 = ax.errorbar(plk_ell, plk_l4cl, yerr=plk_err, fmt='o', markersize='0', elinewidth=1.5, capsize=1.5, capthick=1.5, 
color='black', zorder=5)
#ax.errorbar(spt_ell, spt_l4cl*1.01, yerr=spt_err, fmt='o', markersize='0', elinewidth=1.5, capsize=1.5, capthick=1.5, 
#color='blue', zorder=5)

#ax.errorbar(670,6e8, yerr=0.6e8, fmt='o', markersize='0', elinewidth=1.5, capsize=1.5, capthick=1.5, 
#color='black', zorder=5)
#plt.text(750, 5.35e8, 'Planck 2015', horizontalalignment='left', verticalalignment='bottom', fontsize=18)

#ax.errorbar(670,4e8, yerr=0.6e8, fmt='o', markersize='0', elinewidth=1.5, capsize=1.5, capthick=1.5, 
#color='blue', zorder=5)
#plt.text(750, 3.35e8, 'SPT-SZ', horizontalalignment='left', verticalalignment='bottom', fontsize=18)

ax.spines['top'].set_color('white')
ax.spines['right'].set_color('white')

ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.set_xlim([1250,1850])
ax.set_ylim([1e9,1.9e9])

plt.yticks([0.9e9,1.5e9],['10','15'])

plt.legend(frameon=False, prop={'size':18})

ax.set_xlabel(r'$\ell$', size=20)
ax.set_ylabel(r'$\ell^2\mathcal{D}_\ell\ \mathrm{[10^8\mu K^2]}$', size=20)
ax.set_title(r'fixing $\Omega_b h^2$, $z_{\rm EQ}$, $\theta_s$, $\theta_d$', fontsize=25, loc='left')

plt.savefig(pdf_file, format='pdf')
