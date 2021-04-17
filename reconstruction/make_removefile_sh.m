function make_removefile_sh(reconOutfolder,numFrames)

for m = 0:(numFrames-1)

toDel_sh = [reconOutfolder,'remove_raws_f',num2str(m),'.sh']; 


fid = fopen(toDel_sh,'w');

str = ['rm ',reconOutfolder,'lm_reorder_f',num2str(m),'_prompts.lm ', reconOutfolder,'lm_reorder_f',num2str(m),'_randoms.lm ', reconOutfolder,'lm_reorder_f',num2str(m),'_sub.lm\n\n']; 
fprintf(fid,str);
str = '';


str = ['rm ',reconOutfolder,'lmrecon_tof_OSEM_f',num2str(m),'.intermediate.1 ', reconOutfolder,'lmrecon_tof_OSEM_f',num2str(m),'.intermediate.2 ', reconOutfolder,'lmrecon_tof_OSEM_f',num2str(m),'.intermediate.3 \n\n']; 
fprintf(fid,str);
str = '';

str = ['rm ',reconOutfolder,'run_lmrecon_explorer_f',num2str(m),'.sh\n\n']; 
fprintf(fid,str); 
str = ''; 

str = ['rm ', reconOutfolder,'sinogram_ssrb_f',num2str(m),'_randoms.raw\n\n']; 
fprintf(fid,str); 
str = ''; 

str = ['rm ',reconOutfolder,'lmrecon_tof_OSEM_f',num2str(m),'.likelihood\n\n']; 
fprintf(fid,str); 
str = ''; 




fclose(fid); 

end

