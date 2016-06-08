import numpy as np
from figure_sptxhfi import *
from hpylib.util.sptsz_end2end import *


run_name = 'run_06p4'
end_base_file = '/project/kicp/zhenhou/projects/sptxhfi/pspec/'+run_name+'/bandpower_spt_sn_hfi_sn/end_combined_spt150s_hfi143s.sav'

end_150x143_file = '/project/kicp/zhenhou/projects/sptxhfi/pspec/'+run_name+'/bandpower_spt_sn_hfi_sn/end_combined_spt150sn_hfi143sn.sav'

end_150x150_file = '/project/kicp/zhenhou/projects/sptxhfi/pspec/bandpower_spt_sn_spt_sn/end_combined_spt150sn_spt150sn.sav'

end_base = restore_end_save(end_base_file)
end_150x143 = restore_end_save(end_150x143_file)
end_150x150 = restore_end_save(end_150x150_file)

rescale = 1.0087800

end_150x143['dbs_data'] = end_base['dbs_data']

#residual = create_residual_figure(end_150x143, end_150x150, rescale1=rescale, rescale2=rescale**2)

#residual.process_end()
