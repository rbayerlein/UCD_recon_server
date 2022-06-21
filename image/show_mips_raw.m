function mips = show_mips_raw(img_size)

[img_fname, img_path] = uigetfile('*', 'Choose raw image file'); 

img_fname = fullfile(img_path, img_fname); 

fid_img = fopen(img_fname,'r'); 
img = fread(fid_img, inf, 'float'); 
img = reshape(img, img_size); 


img_mask = ones(img_size(1), img_size(2)); 

for xx = 1:img_size(1)
	for yy = 1:img_size(2)
		r = sqrt((xx - 120)^2 + (yy - 120)^2); 
		if r > 110
			img_mask(xx,yy) = 0; 
		end
	end
end


img_mask = repmat(img_mask, 1, 1, img_size(3)); 

img = img .* img_mask; 



mip_cor = squeeze(max(img, [], 1));
mip_sag = squeeze(max(img, [], 2));  

figure
imagesc(mip_cor, [0 max(mip_cor(:))]); 


figure
imagesc(mip_sag, [0 max(mip_sag(:))]); 



figure
imagesc( img(:,:,round(img_size(3)/2)) ); 





