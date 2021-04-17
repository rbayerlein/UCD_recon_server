function make_explorer_lmacc_config(reconOutfolder,numFrames)

fname_cfg = ''; 

for m = 0:(numFrames-1)


fname_cfg = ''; 
fname_cfg = [reconOutfolder,'lmacc_scanner_parameter_f',num2str(m),'.cfg'];

fid1 = fopen(fname_cfg,'w'); 



str = 'detector_ring_diameter = 786.0\n\n';
fprintf(fid1,str);
str = ''; 


str = 'crystal_size = 2.85, 2.85, 2.85, 18.1\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'crystal_gap_size = 0.0, 0.0, 0.0\n\n'; 
fprintf(fid1,str); 
str = ''; 

 
str = 'crystal_array_size = 35, 679\n\n'; 
fprintf(fid1,str); 
str = ''; 

 
str = 'number_of_detector_modules = 24, 1\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'TOF_information = 430, 39.0625\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'number_of_radial_bins = 549\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'image_size = 239, 239, 679\n\n';
fprintf(fid1,str); 
str = ''; 


str = 'voxel_size = 2.85, 2.85, 2.85\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'sensitivity = /run/media/meduser/data/lmrecon_explorer/uih3_sensimage_239s679-voxelsize2_85mm-sum-brd0_4-single.smap\n\n'; 
fprintf(fid1,str); 
str = ''; 


%str = 'iPSF_model = /run/media/meduser/data/lmrecon_explorer/aker_679x679_rd339.pmat\n\n'; 


str = 'iPSF_model = /run/media/meduser/data/lmrecon_explorer/tker_57121x57121.pmat, /run/media/meduser/data/lmrecon_explorer/aker_679x679_rd339.pmat\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'iterative_algorithm_type = 0\n\n';
fprintf(fid1,str); 
str = ''; 

str = 'warmup_setting = -1, 1\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_strength = 0\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_model_type = 0\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_potential_function_type = 2   #0: Quadratic, 1: Hyperbola, 2: Fair, 3: Huber, 4: hw\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_neighborhood_properties = 0, 1    #size, isotropic  #0: 3x3x3, 1:isotropic  # 1st(0) or 2nd(1), aniso(0) or iso(1)\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_buildin_parameter_list = 1e-10\n\n'; 
fprintf(fid1,str); 
str = '';


str = 'iteration_setting = 20, 1, 3    #  subsets, save step, iteration number\n\n'; 
fprintf(fid1,str); 
str = '';


%     #  KEM   2

%# initial_guess = ./lmrecon_tof_475x475x1355_200MBq_20min_PSF.intermediate.1


str = ['input_raw_data_file = ',reconOutfolder,'lm_reorder_f',num2str(m),'_prompts\n\n']; 
fprintf(fid1,str); 
str = ''; 

%input_raw_data_file = ./prompts

%# /home/nas_0_10_ssd/data/xzzhang/UIH_data/RawData_EB/lm_reorder_f$[$m]_prompts


str = 'input_raw_data_format_type = 0\n\n'; 
fprintf(fid1,str); 
str = ''; 

str = ['reconstruction_output_setting = ',reconOutfolder,', ./lmrecon_tof_OSEM_f',num2str(m),'\n\n']; 
fprintf(fid1,str); 


fclose(fid1); 



end
