function run_attn_cor

img_path = uigetdir('Choose image folder');

pp = genpath('../../image/'); 
addpath(pp);  

[img, dcmimg_info] = dcm2raw_fast(img_path); 

rmpath(pp); 

img = double(img); 
 
 
img_size = size(img);  

%img = img + double(dcmimg_info.RescaleIntercept);  

img(img < 0) = 0; 

vox_size = [double(dcmimg_info.PixelSpacing(1)) double(dcmimg_info.PixelSpacing(1)) double(dcmimg_info.SliceThickness)]

sum_img1 = sum(img(:))

while img_size(3) > 828.5
	img = img(:,:,1:2:end) + img(:,:,2:2:end); 
	img = img ./ 2; 
	vox_size(3) = vox_size(3) * 2;
	img_size = size(img)
end



 

while img_size(1) > 256.5

	img = img(1:2:end,:,:) + img(2:2:end,:,:); 
	img = img(:,1:2:end,:) + img(:,2:2:end,:); 
	img = img ./ 4; 
	img_size = size(img); 

	vox_size(1) = vox_size(1)*2; 
	vox_size(2) = vox_size(2)*2;
end

 


% rotate the image to match UCD reconstruction
for zz = 1:img_size(3)
	img(:,:,zz) = rot90(img(:,:,zz), 2); 
end

%img = img - 1000; 

sum_img2 = sum(img(:))


CT_kVp = double(dcmimg_info.KVP)


img_size
vox_size

b = 4.71e-3; 
a = 5.1e-6; 

if CT_kVp == 80
	b = 6.26e-3; 
	a = 3.64e-6; 
elseif CT_kVp  == 120
	b = 4.71e-3; 
	a = 5.1e-6; 
elseif CT_kVp == 140
	b = 4.08e-3; 
	a = 5.64e-6; 
else
	errordlg('Invalid CT kVp'); 
	return; 
end



mu_map = img;  
mu_map(img<1050) = img(img < 1050) .* (9.6e-6); 
mu_map(img>=1050) = b + (a .* img(img >= 1050)); 


lmdata_path = '/run/media/meduser/data/MARTINEZ_RAYMOND_7429193_132851/PET/RawData/1.2.156.112605.18587648329783.191206212852.9.7332.115298/test_attn2/lm_reorder_f1_prompts.lm'; 
fproj_path = '/run/media/meduser/data/MARTINEZ_RAYMOND_7429193_132851/PET/RawData/1.2.156.112605.18587648329783.191206212852.9.7332.115298/test_attn2/lm_reorder_f1_prompts.atn_fac';


lm_mumap_forwardproj(lmdata_path, fproj_path, mu_map, img_size, vox_size)
	
	
	

	
	
