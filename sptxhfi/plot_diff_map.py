import os
import matplotlib.pyplot as plt
from figure_sptxhfi import *

home = os.getenv('HOME')
host = os.getenv('ENV_HOSTNAME')

fits_file = home+'/data_'+host+'/spt_data/coadds/coadd_ra23h30dec-55_2008_150_lps14_0702.fits'
mf_spt = create_map_figure(fits_file)

xra = [mf_spt.ra0dec0[0]-5.00, mf_spt.ra0dec0[0]+5.00] # +/- 8 deg
yra = [mf_spt.ra0dec0[1]-5.00, mf_spt.ra0dec0[1]+5.00] # +/- 8 deg

fits_file = home+'/data_'+host+'/projects/spt_x_planck/planck_2015/scan_filter/coadd_scanmap_fits/full/coadd_scanmap_planck_143_R2.00_full_ra23h30dec-55_2008.fits'
mf_hfi = create_map_figure(fits_file)

map2d_spt = np.flipud(mf_spt.cut_map(xra,yra)) * 1e6 * 0.825
map2d_hfi = np.flipud(mf_hfi.cut_map(xra,yra)) * 1e6

yticks = [-60, -55, -50]
xticks = [348, 352, 356]

vmax = 50
vmin = -50

fig, ax = plt.subplots()

im = plt.imshow(map2d_spt-map2d_hfi, cmap=plt.get_cmap('bone'), extent=(xra[0],xra[1],yra[0],yra[1]), vmax=vmax, vmin=vmin,
                interpolation='bicubic', aspect='equal')

divider = make_axes_locatable(ax)
cax = divider.append_axes("right", size="5%", pad=0.2)

cbar = fig.colorbar(im, ticks=(-50,-25,0,25,50), orientation='vertical', cax=cax,
                    drawedges=False)

cbar.solids.set_edgecolor("face")
cbar.ax.set_ylabel('$\Delta T\;[\mu\mathrm{K}]$', fontsize=16)
cbar.ax.axes.set_yticklabels(["$-50$","$-25$","$0$","$25$","$50$"], fontsize=16)
cbar.outline.set_edgecolor('black')
cbar.outline.set_linewidth(1.2)

#plt.setp(ax.xaxis.get_ticklines(), 'markersize', 3, 'markeredgewidth', 2)
#plt.setp(ax.yaxis.get_ticklines(), 'markersize', 3, 'markeredgewidth', 2)

ax.set_xticks(xticks)
ax.set_yticks(yticks)

ax.set_xlabel("$\mathrm{RA}$", fontsize=16)
ax.set_ylabel(' ', fontsize=16)

ax.axes.set_xticklabels(["$348^{\circ}$","$352^{\circ}$","$356^{\circ}$"], fontsize=16)
ax.axes.set_yticklabels([" ", " "," "," "," "], fontsize=16)

fig_file = 'map_diff.pdf'
plt.savefig(fig_file, format='pdf', transparent=True)
plt.clf()
