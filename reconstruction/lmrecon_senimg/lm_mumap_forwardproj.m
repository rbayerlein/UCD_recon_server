function lm_mumap_forwardproj(lmdata_path, fproj_path, attn_img_raw, attn_img_size, attn_vox_size)

p=genpath('./PETsystem');
addpath(p);
scanner = buildPET('explorer2000mm_unitedimaging');






fid_lm = fopen(lmdata_path, 'rb');
fseek(fid_lm, 0, 'eof'); 
pos_lm = ftell(fid_lm); 
num_events = pos_lm / 10; 
fseek(fid_lm, 0, 'bof'); 


max_events = 1e9; 


fid_out = fopen(fproj_path, 'w'); 


if num_events <= max_events
	lmdata = fread(fid_lm, inf, 'int16'); 
	lmdata = reshape(lmdata, 5, length(lmdata)/5); 
	lmdata_attn = scanner.doListModeForwardProjectionNonTOF(attn_img_raw, attn_img_size, attn_vox_size, int16(lmdata));
	%lmdata_attn = exp(-1 .* lmdata_attn); 
	fwrite(fid_out, lmdata_attn, 'float'); 
else
	event_counter = 0; 
	while event_counter < (num_events - 1)
		if (num_events - event_counter < max_events)
			lmdata = fread(fid_lm, inf, 'int16'); 
			lmdata = reshape(lmdata, 5, length(lmdata)/5);
			event_counter = event_counter + size(lmdata, 2); 
		else
			lmdata = fread(fid_lm, 10*max_events, 'int16'); 
			lmdata = reshape(lmdata, 5, length(lmdata)/5);
			event_counter = event_counter + size(lmdata, 2); 
		end
		lmdata_attn = scanner.doListModeForwardProjectionNonTOF(attn_img_raw, attn_img_size, attn_vox_size, int16(lmdata));
		
		%lmdata_attn = exp(-1 .* lmdata_attn); 
		fwrite(fid_out, lmdata_attn, 'float'); 
	end
end
	
fclose('all'); 

rmpath(p); 







