#!/bin/bash
mkdir cfg_attn_fp_exp
mkdir lm_attn_fp_exp
for (( m=1; m<=8; m++ )) do
{
    
	echo "

detector_ring_diameter = 786.0
crystal_size = 2.85, 2.85, 2.85, 18.1
crystal_gap_size = 0.0, 0.0, 0.0 
crystal_array_size = 35, 679 
number_of_detector_modules = 24, 1 
# TOF_information = 420, 39.06
number_of_radial_bins = 549

image_size = 239, 239, 679
voxel_size = 2.85, 2.85, 2.85

iterative_algorithm_type = 0

initial_guess = /home/rbayerlein/data/explorer/20210827/Multi-Bed_Phantom_Multi-Bed_Phantom_154523/UCD/Image/CTAC_201_mumap_kVp-140_size-256x256x828_vox-2.7344x2.7344x2.344.img


warmup_setting = 1, 1
iteration_setting = 1, 1, 1    #  1, 10, 500

regularizer_strength = 0
regularizer_model_type = 0    #0: pairwise MRF, 1: patch
regularizer_potential_function_type = 2   #0: Quadratic, 1: Hyperbola, 2: Fair, 3: Huber, 4: hw
regularizer_neighborhood_properties = 0, 1    #size, isotropic  #0: 3x3x3, 1:isotropic  # 1st(0) or 2nd(1), aniso(0) or iso(1)
regularizer_buildin_parameter_list = 1e-10

input_raw_data_file = /home/rbayerlein/ssd/IDQVJTYNOC/lm_reorder_f0_prompts.$[$m].lm

input_raw_data_format_type = 0

reconstruction_output_setting =./lm_attn_fp_exp, ./lm_prompts_f$[$m]_attn_fp_exp.raw

	 " >  ./cfg_attn_fp_exp/lmrecon_attn_fp_exp_f$m.cfg

}
done






