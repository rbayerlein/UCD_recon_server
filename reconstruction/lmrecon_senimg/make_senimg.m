function make_senimg(senimg_name, att_image, att_voxel_size, cryseff_wgap, plaeff_wgap)

num_voxels = 239;   
num_slices = 679;  
imagevoxelsize = 2.85; 

image_size = [num_voxels, num_voxels, num_slices];
voxel_size = [imagevoxelsize, imagevoxelsize, imagevoxelsize];

att_image_size = size(att_image); 

padd=genpath('./PETsystem');
addpath(padd);

scanner=buildPET('explorer2000mm_unitedimaging');
obj = scanner;

% num_bins_radial = 549;
num_bins_axial_reduced = 84;
% num_bins_angualr_reduced = 35;
blockringdiff = 4;   


num_blockrings = 8;
num_gap = 7;
num_axialxtals = 672;
num_transxtals = 840;
num_axialxtals_wgap = num_axialxtals+num_gap;


%fname_plaeff_wgap = './plane_efficiency_679x679_float.NC'
%plaeff_wgap = fread(fopen(fname_plaeff_wgap, 'rb'), inf, 'float'); 
%plaeff_wgap = reshape(plaeff_wgap, num_axialxtals_wgap, num_axialxtals_wgap);

%fname_cryseff_wgap = './crystal_efficiency_679_840_float.NC'
%cryseff_wgap = fread(fopen(fname_cryseff_wgap, 'rb'), inf, 'float'); 
%cryseff_wgap = reshape(cryseff_wgap, num_axialxtals_wgap, num_transxtals);



fprintf('calculate the sensitivity image:\n');

numpix_attn = 256; %255;
numslices_attn = 828; %672;

numpix_attn = size(att_image, 1); 
numslices_attn = size(att_image, 3); 

%fprintf('fread attmap data:\n');
%fname_mumap = './ct_ac.img'
%att_image = fread(fopen(fname_mumap, 'rb'), inf, 'float'); 
%att_image = reshape(att_image, numpix_attn, numpix_attn, numslices_attn);
%att_image = att_image(end:-1:1, :, :);

%bp1 = zeros(num_voxels, num_voxels, num_slices);

parpool(5)

parfor blockringdiff_no = 0:blockringdiff
    if blockringdiff_no == 0
    	bp0 = scanner.cal_senimg_uih(plaeff_wgap, cryseff_wgap, att_image, att_image_size, att_voxel_size, image_size, voxel_size, 0, num_blockrings, num_bins_axial_reduced);  
    end
    if blockringdiff_no == 1
    	bp1 = scanner.cal_senimg_uih(plaeff_wgap, cryseff_wgap, att_image, att_image_size, att_voxel_size, image_size, voxel_size, 1, num_blockrings, num_bins_axial_reduced);
    end
    if blockringdiff_no == 2
    	bp2 = scanner.cal_senimg_uih(plaeff_wgap, cryseff_wgap, att_image, att_image_size, att_voxel_size, image_size, voxel_size, 2, num_blockrings, num_bins_axial_reduced);
    end
    if blockringdiff_no == 3
    	bp3 = scanner.cal_senimg_uih(plaeff_wgap, cryseff_wgap, att_image, att_image_size, att_voxel_size, image_size, voxel_size, 3, num_blockrings, num_bins_axial_reduced);
    end
    if blockringdiff_no == 4
    	bp4 = scanner.cal_senimg_uih(plaeff_wgap, cryseff_wgap, att_image, att_image_size, att_voxel_size, image_size, voxel_size, 4, num_blockrings, num_bins_axial_reduced);
    end
end
bp = bp0 + bp1 + bp2 + bp3 + bp4;

fwrite(fopen(senimg_name, 'w'), bp, 'single');
fclose('all');


rmpath('padd'); 



