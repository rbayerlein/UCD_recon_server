function sum_sino_block(outfolder, frame)

encoding_type = 'int32'; 

sino_size = [91, 60, 223]; 

sino_all8_prompts = zeros(sino_size(1)*sino_size(2)*sino_size(3), 1); 
sino_all8_randoms = zeros(sino_size(1)*sino_size(2)*sino_size(3), 1);

for n = 1:8
n
	sino_name_p = ['sinogram_block_f', num2str(frame), '_prompts.', num2str(n), '.raw'];
	sino_name_r = ['sinogram_block_f', num2str(frame), '_randoms.', num2str(n), '.raw'];  
	fname_in_p = fullfile(outfolder, sino_name_p); 
	fname_in_r = fullfile(outfolder, sino_name_r);
	
	fid_p = fopen(fname_in_p, 'r'); 
	fid_r = fopen(fname_in_r, 'r');
	
	sino_p = []; 
	sino_r = []; 
	sino_p = fread(fid_p, inf, encoding_type);
	sino_r = fread(fid_r, inf, encoding_type); 
	
	size(sino_p)
	size(sino_all8_prompts)
	
	
	sino_all8_prompts = sino_all8_prompts + sino_p; 
	sino_all8_randoms = sino_all8_randoms + sino_r; 
	
	fclose(fid_p); 
	fclose(fid_r); 
	
	delete(fname_in_p); delete(fname_in_r); 
	
end

sino_name_p = ['sinogram_block_f', num2str(frame), '_prompts.raw'];
sino_name_r = ['sinogram_block_f', num2str(frame), '_randoms.raw'];

fname_out_p = fullfile(outfolder, sino_name_p); 
fname_out_r = fullfile(outfolder, sino_name_r); 


fid_out_p = fopen(fname_out_p, 'w'); 
fid_out_r = fopen(fname_out_r, 'w'); 


fwrite(fid_out_p, sino_all8_prompts, encoding_type);
fwrite(fid_out_r, sino_all8_randoms, encoding_type); 

fclose(fid_out_p); fclose(fid_out_r);  

clear sino_all8_prompts; clear sino_all8_randoms; clear sino_p; clear sino_r; 
