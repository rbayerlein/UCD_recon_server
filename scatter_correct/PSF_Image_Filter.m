% filter for reconstructed images when PSF is OFF
% very large image values outside the scanner ring can occur and simset
% will not be able to execute a simulation

function PSF_Image_Filter(fname_in, img_size)

%% get raw image and parameters
fid_in = fopen(fname_in, 'rb');
img_raw = fread(fid_in, inf, 'float');
img_raw = reshape(img_raw, img_size);

img_center = ceil(img_size(1)/2);
filter_radius = 100; % in voxels
filter_factor = 4;% filter any value that is larger than X times the max value in the center region

%% display raw image

% slice = reshape(img_raw(:,:,ceil(img_size(3)/2)), img_size(2), img_size(1));
% imshow(slice, []);
% colorbar;
% colormap('bone')

%% filter extreme values

% get max value along axial center of image

max = 0;
for ax = 1 : img_size(3)
    for x = 1 : img_size(1)
        for y = 1 : img_size(2)
%             fprintf('x=%d, y=%d => val = %d\n',x, y, sqrt(power(x-img_center,2) + power(y-img_center,2)));
            if sqrt(power(x-img_center,2) + power(y-img_center,2)) < filter_radius
                if img_raw(x,y,ax) > max
                   max = img_raw(x,y,ax);
                end
            end
        end
    end
end

fprintf('max value across axial range within radius %d: %d\n', filter_radius, max);

% filter image
img_filtered = zeros(239,239,679);

for ax = 1 : img_size(3)
    for x = 1 : img_size(1)
        for y = 1 : img_size(2)
            if img_raw(x,y,ax) > filter_factor * max
                img_filtered(x,y,ax) = 0;
            else
                img_filtered(x,y,ax) = img_raw(x,y,ax);
            end
        end
    end
end

%% display raw image

% slice_filtered = reshape(img_filtered(:,:,ceil(img_size(3)/2)), img_size(2), img_size(1));
% figure;
% imshow(slice_filtered, []);
% colorbar;
% colormap('bone')

%% save image
disp('saving')
fname_temp = [fname_in, '.temp'];
fid_out = fopen(fname_temp, 'wb');
fwrite(fid_out, img_filtered, 'float');

cmd_rm = ['rm ', fname_in];
system(cmd_rm);
pause(0.1);

cmd_mv = ['mv ', fname_temp, ' ', fname_in];
system(cmd_mv);
pause(0.1);

disp('done')

end