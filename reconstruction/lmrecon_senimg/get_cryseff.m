function crys_eff = get_cryseff(nc_path);

num_crys_trans = 840; 
num_crys_ax = 672; 

fid_nc = fopen(nc_path, 'rb'); 

if fid_nc < 0
	disp('Could not open .nc file, quit'); 
	return; 
end

fseek(fid_nc, 3095517*4, 'bof'); 

crys_eff = fread(fid_nc, num_crys_trans*num_crys_ax, 'float'); 
crys_eff = reshape(crys_eff, num_crys_ax, num_crys_trans); 

fclose(fid_nc); 

