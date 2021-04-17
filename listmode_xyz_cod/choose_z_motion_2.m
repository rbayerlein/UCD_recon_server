function choose_z_motion



imshow('pears.png')


r1 = drawrectangle; 


pos = customWait(r1)



r1pos = get(r1, 'Position')





pause(1)



close all

pause(0.5)

imshow('pears.png')


r2pos_init = r1pos; 
r2pos_init(3) = r2pos_init(3) + 50;
r2pos_init(4) = r2pos_init(4) - 50; 

r2 = drawrectangle('Position', r2pos_init); 
pos = customWait(r2)



r2pos = get(r2, 'Position')


pause


% read dicom image

img_temp = ''; 
i = 1;

while ~contains(img_temp, '.dcm') && i < length(lst)
  img_temp = lst(i).name;
  i = i + 1; 
end 

img_temp = fullfile(udir, img_temp)
dcm_info = dicominfo(img_temp)

 
ring_d = 786;
ax_fov = 1944; 
trans_fov = double(dcm_info.PixelSpacing(1)) * double(dcm_info.Rows);
num_crys = 840; % explorer
num_crys_rings_wgap = 679; 
num_crys_rings_nogap = 672; 
pix_size = [double(dcm_info.PixelSpacing(1)), double(dcm_info.PixelSpacing(1)), double(dcm_info.SliceThickness)];


img_temp = ''; 
i = 1;

img_store = zeros(dcm_info.Rows, dcm_info.Rows, dcm_info.NumberOfSlices); 
zpos_store = zeros(1, dcm_info.NumberOfSlices); 
counter = 1; 

for i = 1:length(lst) 
  img_temp = lst(i).name;
  i = i + 1;
  if contains(img_temp, '.dcm')
    %img_temp = fullfile(udir, img_temp); 
    %dcm_info_temp = dicominfo(img_temp);
    %img_temp = dicomread(dcm_info_temp);
    %img_temp = zeros(600, 600); 
    %img_temp = double(img_temp) .* double(dcm_info_temp.RescaleSlope); 
    %img_store(:,:,counter) =  img_temp; 
    %zpos_store(1, counter) = double(dcm_info_temp.SliceLocation); 
    %counter = counter + 1; 
  end 
end 


%[q, iSort] = sort(zpos_store); 

%zpos_store = zpos_store(iSort); 

%img_store = img_store(:,:,iSort); 

%zpos_store = zpos_store - min(zpos_store); % range 0 to ~1940 mm
%min(zpos_store)
%max(zpos_store)





mip_cor = max(img_store, [], 1); 
mip_sag = max(img_store, [], 2); 


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






function pos = customWait(hROI)

% Listen for mouse clicks on the ROI
l = addlistener(hROI,'ROIClicked',@clickCallback);

% Block program execution
uiwait;

% Remove listener
delete(l);

% Return the current position
pos = hROI.Position;




function clickCallback(~,evt)

if strcmp(evt.SelectionType,'double')
    uiresume;
end




