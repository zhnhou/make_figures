import os
import matplotlib.pyplot as plt
from figure_sptxhfi import *

home = os.getenv('HOME')
host = os.getenv('ENV_HOSTNAME')

fits_file = home+'/data_'+host+'/spt_data/coadds/coadd_ra23h30dec-55_2008_150_lps14_0702.fits'
mf_spt = create_map_figure(fits_file)

xra = [mf_spt.ra0dec0[0]-5.00, mf_spt.ra0dec0[0]+5.00] # +/- 8 deg
yra = [mf_spt.ra0dec0[1]-5.00, mf_spt.ra0dec0[1]+5.00] # +/- 8 deg

print xra
print yra

map2d_spt = np.flipud(mf_spt.cut_map(xra,yra)) * 1e6 * 0.825

yticks = [-60, -55, -50]
xticks = [348, 352, 356]

vmax = 120
vmin = -120

fig, ax = plt.subplots()

im = plt.imshow(map2d_spt, cmap=plt.get_cmap('bone'), extent=(xra[0],xra[1],yra[0],yra[1]), vmax=vmax, vmin=vmin,
                interpolation='bicubic', aspect='equal')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.2)

cbar = fig.colorbar(im, ticks=(-120,-60,0,60,120), orientation='vertical', cax=cax,
                    drawedges=False)

cbar.solids.set_edgecolor("face")
#cbar.ax.set_ylabel('$\Delta T\;[\mu\mathrm{K}]$', fontsize=14)
cbar.ax.set_ylabel('  ', fontsize=16)
cbar.ax.axes.set_yticklabels(["$-120$","$-60$","$0$","$60$","$120$"], fontsize=16)
cbar.outline.set_edgecolor('black')
cbar.outline.set_linewidth(1.2)

#plt.setp(ax.xaxis.get_ticklines(), 'markersize', 3, 'markeredgewidth', 2)
#plt.setp(ax.yaxis.get_ticklines(), 'markersize', 3, 'markeredgewidth', 2)

ax.set_xticks(xticks)
ax.set_yticks(yticks)

#ax.set_xlabel(r"$\mathrm{RA}$", fontsize=16)
ax.set_xlabel(" ", fontsize=16)
ax.set_ylabel('$\mathrm{Dec}$', fontsize=16)

#ax.axes.set_xticklabels(["$348^{\circ}$","$352^{\circ}$","$356^{\circ}$"], fontsize=16)
ax.axes.set_xticklabels([" "," "," "], fontsize=16)
ax.axes.set_yticklabels(["$-60^{\circ}$", "$-55^{\circ}$","$-50^{\circ}$"], fontsize=16)

fig_file = 'map_spt150.pdf'
plt.savefig(fig_file, format='pdf', transparent=True)
plt.clf()
