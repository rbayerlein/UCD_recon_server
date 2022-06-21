#!/bin/bash
study_dir=/home/rbayerlein/data/explorer/20210827/Multi-Bed_Phantom_Multi-Bed_Phantom_154523
server_temp_dir=/home/rbayerlein/ssd/YXZEFSTBRV
cfg_folder=${server_temp_dir}/cfg_attn_fp_exp
lm_folder=${server_temp_dir}/lm_attn_fp_exp

if [[ ! -d $cfg_folder ]]; then
	mkdir $cfg_folder
	chmod -R 775 $cfg_folder
fi

if [[ ! -d $lm_folder ]]; then
	mkdir $lm_folder
	chmod -R 775 $lm_folder
fi

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

initial_guess = ${study_dir}/UCD/Image/CTAC_201_mumap_kVp-140_size-239x239x239_vox-2.85x2.85x2.85.img


warmup_setting = 1, 1
iteration_setting = 1, 1, 1    #  1, 10, 500

regularizer_strength = 0 # mandatory line, otherwise core dumps
regularizer_model_type = 0    #0: pairwise MRF, 1: patch
regularizer_potential_function_type = 2   #0: Quadratic, 1: Hyperbola, 2: Fair, 3: Huber, 4: hw
regularizer_neighborhood_properties = 0, 1    #size, isotropic  #0: 3x3x3, 1:isotropic  # 1st(0) or 2nd(1), aniso(0) or iso(1)
regularizer_buildin_parameter_list = 1e-10

input_raw_data_file = ${server_temp_dir}/lm_reorder_f0_prompts.$[$m].lm

input_raw_data_format_type = 0

reconstruction_output_setting =${lm_folder}, ./lm_prompts_f$[$m]_attn_fp_exp.raw

	 " > ${cfg_folder}/lmrecon_attn_fp_exp_f$m.cfg

}
done






