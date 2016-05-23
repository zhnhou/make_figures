import os
from figure_sptxhfi import *

home = os.getenv('HOME')
host = os.getenv('ENV_HOSTNAME')

fits_file = home+'/data_'+host+'/spt_data/coadds/coadd_ra5h30dec-55_2008_150_lps14_0702.fits'
mf_spt = create_map_figure(fits_file)

xra = [mf_spt.ra0dec0[0]-4.00, mf_spt.ra0dec0[0]+4.00] # +/- 8 deg
yra = [mf_spt.ra0dec0[1]-4.00, mf_spt.ra0dec0[1]+4.00] # +/- 8 deg

map2d_spt = mf_spt.cut_map(xra,yra) * 1e6 * 0.875

xticks = [-59, -57, -55, -53, -51]
yticks = [79, 81, 83, 85]
mf_spt.make_figure(map_image=map2d_spt, vmin=-120, vmax=120)


fits_file = home+'/data_'+host+'/projects/spt_x_planck/planck_2015/scan_filter/coadd_scanmap_fits/full/coadd_scanmap_planck_143_R2.00_full_ra5h30dec-55_2008.fits'
mf_hfi = create_map_figure(fits_file)

map2d_hfi = mf_hfi.cut_map(xra,yra) * 1e6
mf_hfi.make_figure(map_image=map2d_hfi, vmin=-300, vmax=300)
