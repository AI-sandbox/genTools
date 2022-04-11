from gen_tools import *

# User Defined Paramters:
root_dir = '/home/projects/dimreduction/Mexico/masters_Asian_anc4/'
beagle_vcf_file = 'Mexico_AsianAnc_w0.2.beagle'
vit_fbk_tsv_file = 'Mexico_AsianAnc.vit'

beagle_or_vcf = 1 # 1 for beagle / 2 for vcf
vit_or_fbk_or_tsv = 1
fb_or_msp = 1 # Only used if vit_or_fbk_or_tsv <> 1
prob_thresh = 0 # Only used if fb_or_msp = 2

num_arrays = 4
num_ancestries = 4
average_parents = True
is_masked = True

save_masks = True
masks_file = 'masks_file_Asian.npz'

###########################################
masks, rs_ID_list, ind_ID_list = array_process(root_dir, beagle_vcf_file, vit_fbk_tsv_file, beagle_or_vcf, vit_or_fbk_or_tsv, 
                                               fb_or_msp, num_arrays, num_ancestries, average_parents, prob_thresh, is_masked,
                                               save_masks, masks_file)
