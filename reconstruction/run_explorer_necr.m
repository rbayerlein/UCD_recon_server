function run_explorer_recon



num_data = 35; 



fdir_dbase = '/run/media/meduser/backup/20190603/ServicePatient-CrossAcq_ServicePatientID_162249/PET/RawData/1.2.156.112605.18587648329783.190604002242.9.6332.108223/';

fdir_obase = '/run/media/meduser/data/explorer_necr/20190603/';


lst = dir(fdir_dbase);
cc = 1; 
for k = 1:size(lst,1)
	f_temp = lst(k).name; 
	if ~isempty(strfind(f_temp,'.1.raw'))
		A{cc} = f_temp; 
		cc = cc + 1; 
	end
end

whos A

ccc = 1; 
for k = 1:length(A)
	str_test = A{k}
	ind = strfind(str_test,'.1.raw'); 
	qtest = str_test(ind - 5)
	
	f_num = str2num(qtest); 
	if strcmp(str_test(ind-6),'.') < 0.5
		f_num2 = str2num(str_test(ind - 6)); 
		f_num = f_num + 10*f_num2; 
	end
	f_num = floor(f_num / 2) + 1; 
	
	if mod(str2num(qtest), 2) == 1
		B{ccc} = str_test; 
		frame_num(ccc) = f_num; 
		ccc = ccc + 1; 
	end
end

frame_num
	
B{1}
B{2}
B{3}




pause



for N = 1:num_data


recon_out = num2str(frame_num(N)); 

num_frames = 1;
framing = '1,3600'; 

write_sino = '1'; 
write_lm = '0'; 
ext = '0'; 


%fdir_dbase = '/run/media/meduser/data/explorer_timingdata/sensitivity_sensitivity_192708/PET/RawData/1.2.156.112605.18587648329783.190525033146.9.9556.23581/';


%fdir_dbase = '/run/media/meduser/backup/20190603/ServicePatient-CrossAcq_ServicePatientID_162249/PET/RawData/1.2.156.112605.18587648329783.190604002242.9.6332.108223/';

%fdir_dbase = '/run/media/meduser/data/explorer_necr/20190603/raw_data_temp/'; 

 


recon_out = [fdir_obase,recon_out]; 



if ~exist(recon_out,'dir'); 

mkdir(recon_out);

recon_out


%raw_data_fname = [fdir_dbase,'1.2.156.112605.18587648329783.190604002255.9.1940.55689.1.raw'];
raw_data_fname = [fdir_dbase,B{N}];
 
B{N}

 
recon_info_name = '/run/media/meduser/data/read_lm/Reconstruction_Parameters_1'; 

fid=fopen(recon_info_name,'wt');
fprintf(fid,'%s\n%s\n%s\n%s\n%s\n%s\n',...
	recon_out,...
    raw_data_fname,...
    ext,...
    framing,...
    write_lm,...
    write_sino); 
                
fclose(fid);


recon_out = [recon_out,'/']; 


make_explorer_lmacc_config(recon_out,num_frames); 
pause(1)


make_run_recon_explorer(recon_out,num_frames); 
pause(1)


make_removefile_sh(recon_out,num_frames); 
pause(1)


       
fname = '/run/media/meduser/data/read_lm/process_lm_sino_explorer';
system(fname);


end
           
end     
                
                
                
                
                
                
                
