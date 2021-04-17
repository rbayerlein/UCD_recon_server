function plot_cod_xyz(fname_cod_xyz, t_sample, num_roi)

fid = fopen(fname_cod_xyz, 'r'); 

cod_xyz = fread(fid, inf, 'double'); 

if mod(length(cod_xyz), 4) ~= 0
	errordlg('Invalid file size, must be multiple of 4'); 
	return
end



cod_xyz = reshape(cod_xyz, 4, length(cod_xyz)/4); 


cod_xyz(1,:) = cod_xyz(1,:) ./ cod_xyz(4,:); 
cod_xyz(2,:) = cod_xyz(2,:) ./ cod_xyz(4,:); 
cod_xyz(3,:) = cod_xyz(3,:) ./ cod_xyz(4,:); 

cod_xyz = cod_xyz(1:3,:); 


for k = 1:num_roi


cod_xyz_u0 = cod_xyz(:,k:num_roi:end); 


t_plot = 0:(size(cod_xyz_u0, 2));
t_plot = t_sample .* t_plot(1:(end-1)); 

figure
hold on
plot(t_plot, cod_xyz_u0(1,:), 'r');
plot(t_plot, cod_xyz_u0(2,:), 'g');
plot(t_plot, cod_xyz_u0(3,:), 'b'); 
title(num2str(k)); 
hold off

end



