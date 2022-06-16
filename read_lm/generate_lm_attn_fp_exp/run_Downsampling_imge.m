function run_Downsampling_imge(fname_mumap, img_size_out, vox_size_out)
%% input mumap file name
% fname_mumap = '/home/rbayerlein/data/explorer/20210714/NERVO_ATILLIO_7696939_144717_Raw/UCD/Image/CTAC_120_MIN_201_mumap_kVp-140_size-256x256x828_vox-2.7344x2.7344x2.344.img';

fprintf('mumap in:%s\n', fname_mumap);
% %% voxel sizes OUTPUT
% img_size_out = [239,239,679];
% vox_size_out = [2.85,2.85,2.85];
%% voxel sizes INPUT
str1 = '_size-'; 
str2 = '_vox-'; 
ind_ss1 = strfind(fname_mumap, str1); 
ind_ss2 = strfind(fname_mumap, str2); 
ind_end = strfind(fname_mumap, '.img'); 

if isempty(ind_ss1) || isempty(ind_ss2)
	disp('Invalid mu map. Return.'); 
	return; 
else
    disp('mu map is good. continue.');
end
size_str = fname_mumap((ind_ss1+length(str1)):(ind_ss2-1));
vox_str = fname_mumap((ind_ss2+length(str2)):(ind_end-1));
ind_x1 = strfind(size_str, 'x');
ind_x2 = strfind(vox_str, 'x');
img_size_in = [str2double(size_str(1:(ind_x1(1)-1))) str2double(size_str((ind_x1(1)+1):(ind_x1(2)-1))) str2double(size_str((ind_x1(2)+1):end))]
vox_size_in = [str2double(vox_str(1:(ind_x2(1)-1))) str2double(vox_str((ind_x2(1)+1):(ind_x2(2)-1))) str2double(vox_str((ind_x2(2)+1):end))]

%% open files
fid_mumap = fopen(fname_mumap, 'rb');
data_in = fread(fid_mumap, inf, 'float');
data_in = reshape(data_in, img_size_in);
%% run code
data_out = Downsmapling_image(data_in, vox_size_in, img_size_out, vox_size_out);
%% show results
% slice = reshape(data_out( 120, :,:), [239 679]);
% imshow(slice, [])
%% save results
size_pos = strfind(fname_mumap, 'size-');
file_out = [fname_mumap(1:size_pos+4), num2str(img_size_out(1)), 'x', num2str(img_size_out(2)), 'x', num2str(img_size_out(3)), '_vox-',num2str(vox_size_out(1)), 'x', num2str(vox_size_out(2)),'x', num2str(vox_size_out(3)), '.img'];
fprintf('file out:%s\n', file_out);
fid_out = fopen(file_out, 'wb');
fwrite(fid_out, data_out, 'float');
fclose(fid_out);
fclose(fid_mumap);
disp('done');

end