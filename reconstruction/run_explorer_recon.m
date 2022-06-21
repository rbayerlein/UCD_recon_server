function run_explorer_recon



recon_out = 'center/1'; 

num_frames = 1;
framing = '1,300'; 

write_sino = '1'; 
write_lm = '0'; 
ext = '0'; 


%fdir_dbase = '/run/media/meduser/data/explorer_timingdata/sensitivity_sensitivity_192708/PET/RawData/1.2.156.112605.18587648329783.190525033146.9.9556.23581/';


%fdir_dbase = '/run/media/meduser/backup/20190603/ServicePatient-CrossAcq_ServicePatientID_162249/PET/RawData/1.2.156.112605.18587648329783.190604002242.9.6332.108223/';

fdir_obase = '/run/media/meduser/data/explorer_sensitivity/sinos_20190603/'; 

fdir_dbase = '/run/media/meduser/backup/sensitivity_0603_sensitivity_0603_151206/PET/RawData/'; 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1

recon_out = 'center/1'; 

raw_data_fname = [fdir_dbase,'1.2.156.112605.18587648329783.190603231207.9.6956.88361/','1.2.156.112605.18587648329783.190603231406.9.1940.14753.1.raw'];


recon_out = [fdir_obase,recon_out]; 
recon_out
mkdir(recon_out);
 
 
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



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2

recon_out = 'center/2'; 

raw_data_fname = [fdir_dbase,'1.2.156.112605.18587648329783.190603231931.9.6956.26767/','1.2.156.112605.18587648329783.190603232005.9.1940.13939.1.raw'];


recon_out = [fdir_obase,recon_out]; 
recon_out
mkdir(recon_out);
 
 
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
                
                


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3

recon_out = 'center/3'; 

raw_data_fname = [fdir_dbase,'1.2.156.112605.18587648329783.190603232544.9.6956.22799/','1.2.156.112605.18587648329783.190603232642.9.1940.18029.1.raw'];


recon_out = [fdir_obase,recon_out]; 
recon_out
mkdir(recon_out);
 
 
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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4

recon_out = 'center/4'; 

raw_data_fname = [fdir_dbase,'1.2.156.112605.18587648329783.190603233321.9.6956.24853/','1.2.156.112605.18587648329783.190603233324.9.1940.16814.1.raw'];


recon_out = [fdir_obase,recon_out]; 
recon_out
mkdir(recon_out);
 
 
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5

recon_out = 'center/5'; 

raw_data_fname = [fdir_dbase,'1.2.156.112605.18587648329783.190603233841.9.6956.24809/','1.2.156.112605.18587648329783.190603233901.9.1940.14939.1.raw'];


recon_out = [fdir_obase,recon_out]; 
recon_out
mkdir(recon_out);
 
 
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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%1

recon_out = '10cm_offcenter/1'; 

raw_data_fname = [fdir_dbase,'1.2.156.112605.18587648329783.190603234413.9.6956.23351/','1.2.156.112605.18587648329783.190603234435.9.1940.16757.1.raw'];


recon_out = [fdir_obase,recon_out]; 
recon_out
mkdir(recon_out);
 
 
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%2

recon_out = '10cm_offcenter/2'; 

raw_data_fname = [fdir_dbase,'1.2.156.112605.18587648329783.190603235001.9.6956.27061/','1.2.156.112605.18587648329783.190603235032.9.1940.13894.1.raw'];


recon_out = [fdir_obase,recon_out]; 
recon_out
mkdir(recon_out);
 
 
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


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%3

recon_out = '10cm_offcenter/3'; 

raw_data_fname = [fdir_dbase,'1.2.156.112605.18587648329783.190603235551.9.6956.27893/','1.2.156.112605.18587648329783.190603235622.9.1940.13695.1.raw'];


recon_out = [fdir_obase,recon_out]; 
recon_out
mkdir(recon_out);
 
 
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





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%4

recon_out = '10cm_offcenter/4'; 

raw_data_fname = [fdir_dbase,'1.2.156.112605.18587648329783.190604000137.9.6956.25098/','1.2.156.112605.18587648329783.190604000211.9.1940.11103.1.raw'];


recon_out = [fdir_obase,recon_out]; 
recon_out
mkdir(recon_out);
 
 
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




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%5

recon_out = '10cm_offcenter/5'; 

raw_data_fname = [fdir_dbase,'1.2.156.112605.18587648329783.190604000731.9.6956.23974/','1.2.156.112605.18587648329783.190604000824.9.1940.12936.1.raw'];


recon_out = [fdir_obase,recon_out]; 
recon_out
mkdir(recon_out);
 
 
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
                
                
                
                
                
                
