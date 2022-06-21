function show_mips_dcm
img_path = uigetdir('Choose image folder'); 

[img, dcmimg_info] = dcm2raw_fast(img_path, 3); 
 
 
img_size = size(img);  

%img_mask = ones(img_size(1), img_size(2)); 

%for xx = 1:img_size(1)
%	for yy = 1:img_size(2)
%		r = sqrt((xx - 120)^2 + (yy - 120)^2); 
%		if r > 110
%			img_mask(xx,yy) = 0; 
%		end
%	end
%end


%img_mask = repmat(img_mask, 1, 1, img_size(3)); 

%img = img .* img_mask; 

%for zz = 1:img_size(3)
%	img(:,:,zz) = rot90(img(:,:,zz), 3); 
%end



mip_cor = squeeze(max(img, [], 1));
mip_sag = squeeze(max(img, [], 2));  

figure
imagesc(mip_cor, [0 max(mip_cor(:))]); 


figure
imagesc(mip_sag, [0 max(mip_sag(:))]); 



figure
imagesc( img(:,:,round(img_size(3)/2)) ); 
