import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx

import numpy as np
import cPickle as pickle
from hpylib.util.remote_data import *
from hpylib.cosmology.read_planck_data import *

pdf_file = 'pk_mnu.pdf'

dl_mnu_file = '~/data_midway/projects/neutrino/neff/param_sample/Dl_mnu_const_omegab_zeq_thetas.pkl'
tmp = pickle.load(open(sync_from_remote('midway',dl_mnu_file,update=False),'rb'))

num_sample = tmp['num_sample']
pk = tmp['pk']
kh = tmp['kh']

mnu=np.arange(0,3.0125,0.0125)
jet = cm = plt.get_cmap('jet') 
cNorm  = colors.Normalize(vmin=mnu[0], vmax=mnu[-1])
scalarMap = cmx.ScalarMappable(norm=cNorm, cmap=jet)
scalarMap._A = []

fig, ax = plt.subplots()
ax.set_position([0.1,0.15,0.88,0.65])

for i in np.arange(len(mnu)):
    colorVal = scalarMap.to_rgba(mnu[i])
    plt.plot(kh[:,i], pk[:,i], color=colorVal, zorder=1, linewidth=0.01)

cbar = plt.colorbar(scalarMap, ticks=[0,1.0,2.0,3.0])
cbar.set_label(r'$\Sigma m_\nu\ [\,\rm{eV}\,]$', fontsize=18)

ax.set_xscale('log')
ax.set_yscale('log')

ax.spines['top'].set_color('white')
ax.spines['right'].set_color('white')

ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.set_xlim([1e-4,1])
ax.set_ylim([10,6e4])

#plt.yticks([0.5e9,1e9,1.5e9],['5','10','15'])
#plt.legend(frameon=False, prop={'size':18}, loc=3)

ax.set_xlabel(r'$k\ [\,\mathrm{h\;Mpc^{-1}}\,]$', size=19)
ax.set_ylabel(r'$P(k)\ [\,\mathrm{h^{-3}\;Mpc^3}\,]$', size=20)
ax.set_title(r'fixing $\Omega_b h^2$, $\Omega_c h^2$, $\theta_s$', fontsize=25, loc='left')

plt.savefig(pdf_file, format='pdf')
