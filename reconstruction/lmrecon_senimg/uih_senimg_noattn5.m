clear all

p=genpath('./PETsystem');
addpath(p);

scanner=buildPET('explorer2000mm_unitedimaging');
obj = scanner;

% scanner.system_parms



% num_bins_radial = 549;
num_bins_axial_reduced = 84;
% num_bins_angualr_reduced = 35;
blockringdiff = 4;   




num_blockrings = 8;
num_gap = 7;
num_axialxtals = 672;
num_transxtals = 840;
num_axialxtals_wgap = num_axialxtals+num_gap;


fname_plaeff_wgap = './plane_efficiency_679x679_float.NC'
plaeff_wgap = fread(fopen(fname_plaeff_wgap, 'rb'), inf, 'float'); 
plaeff_wgap = reshape(plaeff_wgap, num_axialxtals_wgap, num_axialxtals_wgap);

fname_cryseff_wgap = './crystal_efficiency_679_840_float.NC'
cryseff_wgap = fread(fopen(fname_cryseff_wgap, 'rb'), inf, 'float'); 
cryseff_wgap = reshape(cryseff_wgap, num_axialxtals_wgap, num_transxtals);



fprintf('calculate the sensitivity image:\n');

numpix_attn = 255;
numslices_attn = 672;

%fprintf('fread attmap data:\n');
%fname_mumap = '../UMAP/mumap_bed_0.map'
%att_image = fread(fopen(fname_mumap, 'rb'), inf, 'float'); 
%att_image = reshape(att_image, numpix_attn, numpix_attn, numslices_attn);
%att_image = att_image(end:-1:1, :, :);

att_image_size = [numpix_attn numpix_attn numslices_attn];
attvoxelsize_x = 2.44;
attvoxelsize_z = 2.85;
att_voxel_size = [attvoxelsize_x  attvoxelsize_x  attvoxelsize_z]; 
att_image = zeros(numpix_attn, numpix_attn, numslices_attn); 


num_voxels = 239;   
num_slices = 679;  

image_size = [num_voxels, num_voxels, num_slices];
imagevoxelsize = 2.85;   
voxel_size = [imagevoxelsize, imagevoxelsize, imagevoxelsize];


bp1 = zeros(num_voxels, num_voxels, num_slices);


for blockringdiff_no = 4:4

    bp0 = scanner.cal_senimg_uih(plaeff_wgap, cryseff_wgap, att_image, att_image_size, att_voxel_size, image_size, voxel_size, blockringdiff_no, num_blockrings, num_bins_axial_reduced); 

    size(bp0)
    
    bp1 = bp1 + bp0;  
    %bp1(:) = bp1(:) + bp0;

end


fwrite(fopen(strcat('uih3_sensimage_', num2str(num_voxels), 's', num2str(num_slices), '-voxelsize2_85mm-sum-brd0_4-single_noattn4.smap'), 'w'), bp1, 'single');
fclose('all');






