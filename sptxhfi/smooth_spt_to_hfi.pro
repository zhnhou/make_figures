pro smooth_spt_to_hfi

    spt_map = '/home/zhenhou/data_midway/spt_data/coadds/coadd_ra23h30dec-55_2008_150_lps14_0702.fits'
    mask    = '/home/zhenhou/data/spt_data/masks/lps14/mask_ra23h30dec-55_2008.sav'
    hfi_map = '/home/zhenhou/data_midway/projects/spt_x_planck/planck_2015/scan_filter/coadd_scanmap_fits/full/coadd_scanmap_planck_143_R2.00_full_ra23h30dec-55_2008.fits'

    hfi_beam_file = '/project/kicp/zhenhou/planck_data/2015/beams/hfi_beam_143x143_nominal_R2.00.txt'
    spt_beam_file = '/project/kicp/zhenhou/spt_data/beams/blgrid_2008_150_rewrite.txt'

    readcol, hfi_beam_file, il_hfi, bl_hfi, format='(I,D)'
    readcol, spt_beam_file, il_spt, bl_spt, format='(I,D)'

    wpix = healpixwindow(2048)
    bl_effective = bl_hfi * wpix[*,0] / bl_spt

    restore, mask
    m = read_spt_fits(spt_map)
    map_spt = m.map.map * apod * 0.825

    
    m = read_spt_fits(hfi_map)
    map_hfi = m.map.map * apod
    
    ell = findgen(4001)
    map_spt_smoothed = convolve_flatsky(map_spt, 1.00d0, ell, bl_effective)

    tvim, (map_spt_smoothed - map_hfi)*1e6, range=[-100,100]

    save, map_spt_smoothed, map_hfi, filename='/home/zhenhou/data_midway/spt_data/coadds/coadd_ra23h30dec-55_2008_150_smoothed_to_hfi143.sav'

    stop
end
