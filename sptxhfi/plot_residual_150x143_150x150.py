import numpy as np
from figure_sptxhfi import *
from hpylib.util.sptsz_end2end import *
from scipy.io.idl import readsav

run_name = 'run_06p4'
end_base_file = '/project/kicp/zhenhou/projects/sptxhfi/pspec/'+run_name+'/bandpower_spt_sn_hfi_sn/end_combined_spt150s_hfi143s.sav'

end_143x143_file = '/home/zhenhou/data_midway/projects/sptxhfi/pspec/bandpower_hfi_sn_hfi_sn_run_06p4/end_combined_hfi143sn_hfi143sn.sav'

end_150x143_file = '/project/kicp/zhenhou/projects/sptxhfi/pspec/'+run_name+'/bandpower_spt_sn_hfi_sn/end_combined_spt150sn_hfi143sn.sav'

end_150x150_file = '/project/kicp/zhenhou/projects/sptxhfi/pspec/bandpower_spt_sn_spt_sn/end_combined_spt150sn_spt150sn.sav'

beam_cov_residual_file = '/home/zhenhou/Projects/projects/sptxhfi/pspec/residual/beam_cov_residual.sav'

end_base = restore_end_save(end_base_file, ellmin=650, ellmax=2500)
end_143x143 = restore_end_save(end_143x143_file, ellmin=650, ellmax=2500)
end_150x143 = restore_end_save(end_150x143_file, ellmin=650, ellmax=2500)
end_150x150 = restore_end_save(end_150x150_file, ellmin=650, ellmax=2500)

end_150x143['dbs_data'] = end_base['dbs_data']

b0 = readsav(beam_cov_residual_file)
b1 = {}

b1['d21d21'] = b0['d21d21'][13:50,13:50]
b1['d31d31'] = b0['d31d31'][13:50,13:50]
b1['d21d31'] = b0['d21d31'][13:50,13:50]
b1['d31d21'] = b0['d31d21'][13:50,13:50]

residual = create_residual_figure(end_143x143, end_150x143, end_150x150, res_beam_cov=b1)

residual.process_end()
residual.make_residual_figure()



