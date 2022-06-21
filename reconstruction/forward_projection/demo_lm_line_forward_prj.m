%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Forward Projection and Attenuation Factor for uEXPLORER
% Xuezhu Zhang
% Qi Lab
% 2016-2017
% 



p=genpath('../PETsystem');
addpath(p);
scanner = buildPET('explorer2000mm_unitedimaging');


numpix_attn = 255;
numslices_attn = 672;

fprintf('fread attmap data:\n');

dir_atnimage = './UIH_Data/UMAP/'; 
image_atn_map3d = fread(fopen(strcat(dir_atnimage, 'mumap_bed_0_flip_x.map'), 'rb'), inf, 'float');
image_atn_map3d = double(reshape(image_atn_map3d, num_voxels, num_voxels, num_slices));

att_image_size = [numpix_attn numpix_attn numslices_attn];
attvoxelsize_x = 2.44;
attvoxelsize_z = 2.85;
att_voxel_size = [attvoxelsize_x  attvoxelsize_x  attvoxelsize_z]


% listmode events or LORs
% dir_lmdata = './'
% lmdata = fread(fopen(strcat(dir_lmdata, 'lmdata_nontof_idxlor_5xint16.raw'), 'rb'), inf, 'int16');
% lmdata = reshape(lmdata, 5, numrad*numang*numring_wgap*numring_wgap);


lmdata_attn = scanner.doListModeForwardProjectionNonTOF(image_atn_map3d, att_image_size, att_voxel_size, int16(lmdata)); 
data_fname = strcat('lmdata_line_forward_projection_float.raw');
fid = fopen(data_fname, 'w');        
fwrite(fid, lmdata_attn_exp, 'single');        
fclose('all');  


lmdata_attn_exp = exp(-lmdata_attn);

data_fname = strcat('lmdata_attenuation_exp_float.raw');
fid = fopen(data_fname, 'w');        
fwrite(fid, lmdata_attn_exp, 'single');        
fclose('all');  







