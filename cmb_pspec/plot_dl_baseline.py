import matplotlib.pyplot as plt
from hpylib.util.remote_data import *
import cPickle as pickle

pdf_file = 'dl_baseline.pdf'

dl_base_file = '~/data_midway/projects/cmb_pspec/baseline/dl_lensed_base_TT_lowP_lensing.pkl'
dl_lensed = pickle.load(open(sync_from_remote('midway', dl_base_file, update=False), 'rb'))

dl_lensed_TT = dl_lensed[:,0]
lmax = np.shape(dl_lensed_TT)[0]-1

ell = np.arange(lmax+1)

fig, ax = plt.subplots()
ax.set_position([0.1,0.1,0.85,0.85])

plt.plot(ell, dl_lensed_TT, color='black')

ax.spines['top'].set_color('white')
ax.spines['right'].set_color('white')

ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.set_xlim([2,2500])
ax.set_ylim([60,6000])

#ax.set_xscale('log')
ax.set_yscale('log')

#plt.yticks([0.5e9,1e9,1.5e9],['5','10','15'])

ax.set_xlabel(r'$\ell$', size=20)
ax.set_ylabel(r'$\ell^2\mathcal{D}_\ell\ \mathrm{[\mu K^2]}$', size=20)

plt.savefig(pdf_file, format='pdf')
plt.clf()
plt.cla()

"""
change omegabh2
"""

pdf_file = 'dl_change_ombh2.pdf'

dl_ombh2_file = '~/data_midway/projects/cmb_pspec/baseline/dl_lensed_change_ombh2.pkl'
dl_lensed = pickle.load(open(sync_from_remote('midway', dl_ombh2_file, update=False), 'rb'))

dl_lensed_TT = dl_lensed[:,0,:]
lmax = np.shape(dl_lensed_TT)[0]-1
num_sample = np.shape(dl_lensed_TT)[1]

ell = np.arange(lmax+1)

fig, ax = plt.subplots()
ax.set_position([0.1,0.1,0.85,0.85])

colors = ['red','darkorange','deepskyblue','green','purple']
labels = [r'$\omega_b + 20\sigma_{\omega_b}$',
          r'$\omega_b + 10\sigma_{\omega_b}$',
          r'$\omega_b = 0.02226$',
          r'$\omega_b - 10\sigma_{\omega_b}$',
          r'$\omega_b + 20\sigma_{\omega_b}$']
for i in np.arange(num_sample):
    plt.plot(ell, dl_lensed_TT[:,num_sample-i-1], color=colors[i], label=labels[i])

plt.legend(frameon=False, prop={'size':18})

ax.spines['top'].set_color('white')
ax.spines['right'].set_color('white')

ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.set_xlim([2,2500])
ax.set_ylim([60,8000])

#ax.set_xscale('log')
ax.set_yscale('log')

#plt.yticks([0.5e9,1e9,1.5e9],['5','10','15'])

ax.set_xlabel(r'$\ell$', size=20)
ax.set_ylabel(r'$\mathcal{D}_\ell\ \mathrm{[\mu K^2]}$', size=20)

plt.savefig(pdf_file, format='pdf')
plt.clf()
plt.cla()


"""
change omegach2
"""

pdf_file = 'dl_change_omch2.pdf'

dl_ombh2_file = '~/data_midway/projects/cmb_pspec/baseline/dl_lensed_change_omch2.pkl'
dl_lensed = pickle.load(open(sync_from_remote('midway', dl_ombh2_file, update=False), 'rb'))

dl_lensed_TT = dl_lensed[:,0,:]
lmax = np.shape(dl_lensed_TT)[0]-1
num_sample = np.shape(dl_lensed_TT)[1]

ell = np.arange(lmax+1)

fig, ax = plt.subplots()
ax.set_position([0.1,0.1,0.85,0.85])

colors = ['red','darkorange','deepskyblue','green','purple']
labels = [r'$\omega_c + 10\sigma_{\omega_c}$',
          r'$\omega_c + 5\sigma_{\omega_c}$',
          r'$\omega_c = 0.1186$',
          r'$\omega_c - 5\sigma_{\omega_c}$',
          r'$\omega_c + 10\sigma_{\omega_c}$']
for i in np.arange(num_sample):
    plt.plot(ell, dl_lensed_TT[:,num_sample-i-1], color=colors[i], label=labels[i])

plt.legend(frameon=False, prop={'size':18})

ax.spines['top'].set_color('white')
ax.spines['right'].set_color('white')

ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.set_xlim([2,2500])
ax.set_ylim([60,8000])

#ax.set_xscale('log')
ax.set_yscale('log')

#plt.yticks([0.5e9,1e9,1.5e9],['5','10','15'])

ax.set_xlabel(r'$\ell$', size=20)
ax.set_ylabel(r'$\mathcal{D}_\ell\ \mathrm{[\mu K^2]}$', size=20)

plt.savefig(pdf_file, format='pdf')
plt.clf()
plt.cla()


"""
change omegak
"""

pdf_file = 'dl_change_omk.pdf'

dl_ombh2_file = '~/data_midway/projects/cmb_pspec/baseline/dl_lensed_change_omk.pkl'
dl_lensed = pickle.load(open(sync_from_remote('midway', dl_ombh2_file, update=True), 'rb'))

dl_lensed_TT = dl_lensed[:,0,:]
lmax = np.shape(dl_lensed_TT)[0]-1
num_sample = np.shape(dl_lensed_TT)[1]

ell = np.arange(lmax+1)

fig, ax = plt.subplots()
ax.set_position([0.1,0.1,0.85,0.85])

colors = ['red','darkorange','deepskyblue','green','purple']
labels = [r'$\Omega_k + 2\sigma_{\Omega_k}$',
          r'$\Omega_k + \sigma_{\Omega_k}$',
          r'$\Omega_k = 0.00$',
          r'$\Omega_k - \sigma_{\Omega_k}$',
          r'$\Omega_k + 2\sigma_{\Omega_k}$']
for i in np.arange(num_sample):
    plt.plot(ell, dl_lensed_TT[:,num_sample-i-1], color=colors[i], label=labels[i])

plt.legend(frameon=False, prop={'size':18})

ax.spines['top'].set_color('white')
ax.spines['right'].set_color('white')

ax.get_xaxis().tick_bottom()
ax.get_yaxis().tick_left()

ax.set_xlim([2,2500])
ax.set_ylim([60,8000])

#ax.set_xscale('log')
ax.set_yscale('log')

#plt.yticks([0.5e9,1e9,1.5e9],['5','10','15'])

ax.set_xlabel(r'$\ell$', size=20)
ax.set_ylabel(r'$\mathcal{D}_\ell\ \mathrm{[\mu K^2]}$', size=20)

plt.savefig(pdf_file, format='pdf')
plt.clf()
plt.cla()
