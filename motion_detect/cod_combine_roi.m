function cod_combine_roi


fname_cod_xyz = [handles.outfolder_path, 'cod_xyz.cg']; 
fid = fopen(fname_cod_xyz, 'r'); 
cod_xyz = fread(fid, inf, 'double'); 
fclose(fid); 

if isempty(cod_xyz)|| mod(length(cod_xyz), 4) ~= 0
	errordlg('Invalid file size, must be multiple of 4'); 
	return
end

ts = 0.5;  

cod_xyz = reshape(cod_xyz, 4, length(cod_xyz)/4); 
