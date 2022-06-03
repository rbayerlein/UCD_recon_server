function data_out = Downsmapling_image(data_in,vox_size_in,img_size_out,vox_size_out)
% unit of vox_size_in and vox_size_out :mm x mm x mm
% By Zhaoheng Xie  
% this function does using replicate the boundary instead of  adding zero padding 
% check image dim to make sure output image is smaller than input image

xdim_in = size(data_in,1)*vox_size_in(1); % mm
ydim_in = size(data_in,2)*vox_size_in(2); % mm
zdim_in = size(data_in,3)*vox_size_in(3); % mm

xdim_out = img_size_out(1)*vox_size_out(1); % mm
ydim_out = img_size_out(2)*vox_size_out(2); % mm
zdim_out = img_size_out(3)*vox_size_out(3); % mm


while xdim_in < xdim_out
    disp('WARNING: input image dim (x) is smaller than output image dim!');
    disp('Padding the input images...');
    padsize=ceil((xdim_out-xdim_in)/vox_size_in(1)/2)+1; % pad padsize*2 slices at a time,using defalut 'both direction'(symmetric)
    disp(['Symmetrically  padding:' num2str(padsize) 'in both (-x and +x) direction']);
    data_in = padarray(data_in,[padsize,0,0],'replicate');
    xdim_in = size(data_in,1)*vox_size_in(1); % mm
end


while ydim_in < ydim_out
    disp('WARNING: input image dim (y) is smaller than output image dim!');
    disp('Padding the input images...');
    padsize=ceil((ydim_out-ydim_in)/vox_size_in(2)/2)+1; % pad padsize*2 slices at a time,using defalut 'both direction'(symmetric)
    disp(['Symmetrically  padding:' num2str(padsize) 'in both (-y and +y) direction']);
    data_in = padarray(data_in,[0,padsize,0],'replicate');
    ydim_in = size(data_in,2)*vox_size_in(2); % mm
end


while zdim_in < zdim_out
    disp('WARNING: input image dim (z) is smaller than output image dim!');
    disp('Padding the input images...');
    padsize=ceil((zdim_out-zdim_in)/vox_size_in(2)/2)+1; % pad padsize*2 slices at a time,using defalut 'both direction'(symmetric)
    disp(['Symmetrically padding:' num2str(padsize) 'in both (-z and +z) direction']);
    data_in = padarray(data_in,[0,0,padsize],'replicate');
    zdim_in = size(data_in,3)*vox_size_in(3); % mm
end

% create input and query (output) grids (see interp3 documentation for more info)

[X,Y,Z] = meshgrid(-xdim_in/2+vox_size_in(1)/2:vox_size_in(1):xdim_in/2-vox_size_in(1)/2, ...
    -ydim_in/2+vox_size_in(2)/2:vox_size_in(2):ydim_in/2-vox_size_in(2)/2, ...
    -zdim_in/2+vox_size_in(3)/2:vox_size_in(3):zdim_in/2-vox_size_in(3)/2);

[Xq,Yq,Zq] = meshgrid(-xdim_out/2+vox_size_out(1)/2:vox_size_out(1):xdim_out/2 -vox_size_out(1)/2, ...
    -ydim_out/2+vox_size_out(2)/2:vox_size_out(2):ydim_out/2-vox_size_out(2)/2, ...
    -zdim_out/2+vox_size_out(3)/2:vox_size_out(3):zdim_out/2-vox_size_out(3)/2);

data_out = interp3(X,Y,Z,data_in,Xq,Yq,Zq,'cubic' );

end
