function choose_z_motion(dcm_fname_dir)


% read dicom image



 
ring_d = 786;
ax_fov = 1944; 
num_crys = 840; % explorer
num_crys_rings_wgap = 679; 
num_crys_rings_nogap = 672; 




[img_raw, dcm_info] = dcm2raw_fast(dcm_fname_dir); 




% make mips


mip_sag = squeeze(max(img_raw, [], 1)); 
mip_cor = squeeze(max(img_raw, [], 2)); 


% show mips for checking, can be removed
figure
imagesc(mip_cor, [0 0.5*max(mip_cor(:))])
colormap(flipud(gray(1024))); 
axis image; 


figure
imagesc(mip_sag, [0 0.5*max(mip_sag(:))])
colormap(flipud(gray(1024)));
axis image;



%  make one image with both mips

mip_both = cat(1, mip_cor, mip_sag); 
mip_both_sz = size(mip_both); 

f1 = figure
set(f1, 'Units', 'Normalized', 'OuterPosition', [0.1, 0.1, 0.9, 0.9]);
imagesc(mip_both, [0 0.5*max(mip_sag(:))]); 
colormap(flipud(gray(1024)));
axis image; 


r1 = imrect(); 
p = wait(r1) 

r2 = imrect(); 


%r1 = drawrectangle('Position', [floor(mip_both_sz(2)/2), 1, 50, floor(mip_both_sz(1)/2)], 'DrawingArea', [1 1 mip_both_sz(2), floor(mip_both_sz(1)/2)]); 



%r2 = drawrectangle('Position', [floor(mip_both_sz(2)/2), floor(mip_both_sz(1)/2)+2, 50, floor(mip_both_sz(1)/2)-1], 'DrawingArea', [1 1 mip_both_sz(2), floor(mip_both_sz(1)/2)]); 


pause



figure
imagesc(mip_cor); 
axis image; 
[row1, column1] = ginput(1); 


figure
imagesc(mip_sag); 
axis image; 
[row2, column2] = ginput(1); 


slice_select = round((column1 + column2)/2); 

figure
imagesc(mean(img_store(:,:,(slice_select-2):(slice_select+2)), 3)); 
axis image; 
[row3, column3] = ginput(1); 


x_c = (row1 + column3) / 2; 
y_c = (row2 + row3) / 2; 
z_c = (column1 + column2) / 2; 

%voi_center = [round(x_c), round(y_c), round(z_c)]; 
voi_center = [x_c, y_c, z_c]

voi_radius = 30; 
voi_length = 50; 












% steps
% 1. choose center point of voi
% 2. choose size / contour of voi
% 3. 










