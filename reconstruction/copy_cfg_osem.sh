

for (( m=0; m<=10; m++ )) do
{
    
	echo "

detector_ring_diameter = 786.0
crystal_size = 2.85, 2.85, 2.85, 18.1
crystal_gap_size = 0.0, 0.0, 0.0 
crystal_array_size = 35, 679 
number_of_detector_modules = 24, 1 

TOF_information = 430, 39.0625

number_of_radial_bins = 549

image_size = 239, 239, 679
voxel_size = 2.85, 2.85, 2.85


sensitivity = ./uih3_sensimage_239s679-voxelsize2_85mm-sum-brd0_4-single.smap


iPSF_model = ./lmrecon_Siddon_psf_trans_axial_239x239x679/tker_57121x57121.pmat, ./lmrecon_Siddon_psf_trans_axial_239x239x679/aker_679x679_rd339.pmat

iterative_algorithm_type = 0     #  KEM   2

# initial_guess = ./lmrecon_tof_475x475x1355_200MBq_20min_PSF.intermediate.1


input_raw_data_file = ./prompts

# /home/nas_0_10_ssd/data/xzzhang/UIH_data/RawData_EB/lm_reorder_f$[$m]_prompts


input_raw_data_format_type = 0

reconstruction_output_setting = ./, ./lmrecon_tof_OSEM_f$m


	 

	 " >>  ./cfg_187fs/uih_lmrecon_osem_psf_f$m.cfg
       
}
done







