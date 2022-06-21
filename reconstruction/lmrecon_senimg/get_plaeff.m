function plane_eff = get_plaeff(nc_path);

num_crys_ax = 672; 

fid_mich = fopen('../reconstruction/lmrecon_senimg/michel_lut_672x672', 'rb'); 
if fid_mich < 0
	disp('Could not open michelogram LUT, quit'); 
	return;
end

michel_lut_672x672 = fread(fid_mich, inf, 'int'); 
michel_lut_672x672 = reshape(michel_lut_672x672, num_crys_ax, num_crys_ax); 
fclose(fid_mich); 

plane_eff = ones(num_crys_ax,num_crys_ax); 

fid_nc = fopen(nc_path, 'rb'); 

if fid_nc < 0
	disp('Could not open .nc file, quit'); 
	return; 
end

fseek(fid_nc, 2079453*4, 'bof'); 

plane_eff_1 = fread(fid_nc, num_crys_ax*num_crys_ax, 'float');
fclose(fid_nc);

for k = 1:num_crys_ax
	for j = 1:num_crys_ax
		plane_eff(k,j) = plane_eff_1(michel_lut_672x672(k,j)+1); 
	end
end



 
