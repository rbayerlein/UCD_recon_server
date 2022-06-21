function [img_raw, dcm_info] = dcm2raw_fast(dcm_fname_dir, rot_angle)

% Does fast conversion of stack of dicom images to 3D raw image volume
%
% input: 
% 	dcm_fname_dir is path to dicom image directory (e.g. '/mnt/data/rbayerlein/dcm_image_test') 
%
% output:
% 	img_raw is 3D image (matlab variable)
% 	dcm_info is dicom info struct 

% Eric Berg, January 2020
% tested for uEXPLORER only 


% get list of files in directory
lst = dir(dcm_fname_dir); 

img_fname = '';
i = 1;
while ~contains(img_fname, '.dcm') && i < length(lst)
  img_fname = lst(i).name;
  i = i + 1; 
end 

% get the dicom info - image size, rescale slope
img_fname = fullfile(dcm_fname_dir, img_fname); 
dcm_info = dicominfo(img_fname)
trans_fov = double(dcm_info.PixelSpacing(1)) * double(dcm_info.Rows);
pix_size = [double(dcm_info.PixelSpacing(1)), double(dcm_info.PixelSpacing(1)), double(dcm_info.SliceThickness)];
try 
	img_size = [double(dcm_info.Rows) double(dcm_info.Columns) double(dcm_info.NumberOfSlices)];
catch ME
	img_size = [double(dcm_info.Rows) double(dcm_info.Columns) double(dcm_info.ImagesInAcquisition)];
	dcm_info.NumberOfSlices= dcm_info.ImagesInAcquisition; 
end
rescale_slope = double(dcm_info.RescaleSlope); 
rescale_intercept = double(dcm_info.RescaleIntercept); 


% check z position of first, second, and last and second last files to determine file order
% to do
% right now, uEXPLORER images are in correct order in the directory

% get list of filenames in order
lst_new = {}; 
c = 0; 
for k = 1:length(lst)
	img_fname = lst(k).name;
	if contains(img_fname, '.dcm')
		str_temp = img_fname(1:(end-4)); 
		fnum = str2num(str_temp); 
		lst_new{fnum} = img_fname;
		c = c+1; 
	end
end

if round(double(dcm_info.NumberOfSlices)) ~= c
	errordlg('Error when loading dicom images'); 
	return
end


% do binary file read instead of matlab dicomread (SLOW)

img_raw = zeros(img_size);

for i = 1:length(lst_new) 
  img_fname = lst_new{i};
  img_fname = fullfile(dcm_fname_dir, img_fname); 

  numbytesoffset = img_size(1) * img_size(2) * 2; 
  
  fid = fopen(img_fname, 'r');
  fseek(fid, -1*numbytesoffset, 'eof'); 
  img_temp = fread(fid, inf, 'uint16'); 
  fclose(fid); 

  img_temp = (double(img_temp) .* rescale_slope) + rescale_intercept; 
  img_temp = reshape(img_temp, img_size(1:2)); 
  img_raw(:,:, i) = img_temp; 
  img_temp = []; 
end

for zz = 1:size(img_raw,3) 
  img_raw(:,:,zz) = rot90(img_raw(:,:,zz), rot_angle); 
end



mip_sag = squeeze(max(img_raw, [], 1)); 
mip_cor = squeeze(max(img_raw, [], 2)); 


% show mips for checking, can be removed
%figure
%imagesc(mip_cor, [0 0.5*max(mip_cor(:))])
%colormap(flipud(gray(1024))); 
%axis image; 


%figure
%imagesc(mip_sag, [0 0.5*max(mip_sag(:))])
%colormap(flipud(gray(1024)));
%axis image;

