import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

import numpy as np
import cPickle as pickle
from hpylib.util.remote_data import *
from hpylib.cosmology.read_planck_data import *

pdf_file = 'clpp_mnu.pdf'

dl_mnu_file = '~/data_midway/projects/neutrino/neff/param_sample/Dl_mnu_const_omegab_zeq_thetas.pkl'
dl_load = pickle.load(open(sync_from_remote('midway',dl_mnu_file,update=False),'rb'))

num_sample = dl_load['num_sample']
lmax       = dl_load['lmax']
clkk = dl_load['cl_lenspoten'] * 1.0e7

plk_lens = np.loadtxt('/Users/zhenhou/data_local/planck_data/2015/lensing/L4CLphiphi.txt')
xerr = (plk_lens[:,1]-plk_lens[:,0])*0.5
band = (plk_lens[:,1]+plk_lens[:,0])*0.5
clpp_plk = plk_lens[:,2]
yerr = plk_lens[:,3]

mnu=np.arange(0,3.0125,0.0125)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=mnu[0], vmax=mnu[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
scalarMap._A = []

fig, ax = plt.subplots()
ax.set_position([0.1,0.1,0.88,0.60])

ell = np.arange(lmax+1)

for i in np.arange(len(mnu)):
    colorVal = scalarMap.to_rgba(mnu[i])
    plt.plot(ell[2:], clkk[2:,i], color=colorVal, zorder=1, linewidth=0.01)

ax.errorbar(band, clpp_plk, xerr=xerr, yerr=yerr, fmt='o', markersize='0', elinewidth=1.5, capsize=1.5, capthick=1.5, 
color='black', zorder=5, label='Planck 2015')

cbar = plt.colorbar(scalarMap, ticks=[0,1.0,2.0,3.0])
cbar.set_label(r'$\Sigma m_\nu\ [\,\rm{eV}\,]$', fontsize=18)

ax.set_xscale('log')
ax.set_yscale('log')

ax.spines['top'].set_color('white')
ax.spines['right'].set_color('white')

ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.set_xlim([2,2000])
ax.set_ylim([1e-2,2.0])

#plt.yticks([0.5e9,1e9,1.5e9],['5','10','15'])

plt.legend(frameon=False, prop={'size':18}, loc=3)

ax.set_xlabel(r'$L$', size=20)
ax.set_ylabel(r'$[L(L+1)]^2C_L^{\phi\phi}/(2\pi)\ [\times 10^7]$', size=20)
ax.set_title(r'fixing $\Omega_b h^2$, $\Omega_c h^2$, $\theta_s$', fontsize=25, loc='left')

plt.savefig(pdf_file, format='pdf')
