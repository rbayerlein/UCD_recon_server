function make_run_recon_explorer(reconOutfolder,numFrames)

toRun_sh = ''; 

for m = 0:(numFrames-1)

toRun_sh = [reconOutfolder,'run_lmrecon_explorer_f',num2str(m),'.sh']; 


fid = fopen(toRun_sh,'w');

str = 'export OMP_NUM_THREADS=30\n\n';

fprintf(fid,str);

str = '';

str = ['/run/media/meduser/data/lmrecon_explorer/lmrecon_explorer/app/lmrecon_tof  ',reconOutfolder,'lmacc_scanner_parameter_f',num2str(m),'.cfg'];

fprintf(fid,str);

str = '';


fclose(fid);


end





