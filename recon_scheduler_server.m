function recon_scheduler_server

handles_name = '/home/rbayerlein/code/explorer-master/handles_scheduler.mat'; 

pause(0.1); 
h = load(handles_name); 
handles = h.handles;
clear h; 
pause(0.1); 

lm_counter = 1; 

wipe_directories = false; % delete temp folder and clean up lm directory

subsample_on = handles.subsample; 

if subsample_on
  reconFrames_lm = [0,1]; 
else
  reconFrames_lm = [handles.reconFrames, (handles.reconFrames(end)+1)]; % one entry more than actual number of frames
end

lm_data_ready = zeros(length(reconFrames_lm), 1);   % length one larger than actual number of frames


recon_frame_done = zeros(length(handles.reconFrames), 1); 
recon_frame_running = zeros(length(handles.reconFrames), 1); 


%% toBeDeleted
fname_log=[handles.lm_outfolder, '/recon_scheduler_server.log'];
fid_log = fopen(fname_log, 'w');
pause(0.1);

ssss = ['chmod 775 ', fname_log]; 
system(ssss); 
pause(0.1);

fprintf(fid_log, 'start of recon_scheduler_server\n');
%%

fprintf(fid_log, 'copying handles file to lm outfolder\n'); fprintf('copying handles file to lm outfolder\n');
cmd_cp = ['cp ', handles_name, ' ', handles.lm_outfolder];
system(cmd_cp);
pause(0.1);

use_scat_corr=handles.scatter_onoff.Value;
fprintf(fid_log, 'correct scatter: %d\n', double(use_scat_corr));

not_done = true; 
lm_done = false; 



while not_done
  
  % first check if the next in line lm files are all ready, by checking if all next ones are created. 
  if ~lm_done
    ch_lm = true;   % ch_lm = check list mode, i.e. check if list mode files exist
    for k = 1:8
      fname_lm_ch = [handles.lm_outfolder, 'lm_reorder_f', num2str(reconFrames_lm(lm_counter+1)), '_prompts.', num2str(k), '.lm']; 
      if ~exist(fname_lm_ch, 'file');
        ch_lm = false; 
        fprintf(fid_log, 'list mode file %d not done yet\n', k); fprintf('list mode file %d not done yet\n', k);
      end
    end
    if ch_lm

      fprintf(fid_log, 'done list mode files. \n'); fprintf('done list mode files. \n');
      m = reconFrames_lm(lm_counter);
  	  %andles.server_temp_dir{lm_counter} = outfolder_server_temp; 
      outfolder_server_temp = handles.server_temp_dir{lm_counter}; 

  	  % create outfolder server temp (the one with the weird name)
      fprintf(fid_log, 'creating outfolder_server_temp %s\n', outfolder_server_temp); fprintf('creating outfolder_server_temp %s\n', outfolder_server_temp);
      cmd_mkdir = ['cd ', handles.server_recon_data_dir, '; mkdir ',outfolder_server_temp,'; '];
  	  system(cmd_mkdir); 
  	  pause(0.1); 

      if subsample_on
        

        % combine 8 .lm files into 1 for UCD recon
        cmd_combine = ['combine_listmode_subsample ',...
        handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f ',...
        handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.1.lm ',...
        handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.2.lm ',...
        handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.3.lm ',...
        handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.4.lm ',...
        handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.5.lm ',...
        handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.6.lm ',...
        handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.7.lm ',...
        handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.8.lm '];  

        make_explorer_lmacc_config_server(handles, m, outfolder_server_temp);
        toRun_sh = [handles.lm_outfolder,'run_lmrecon_explorer_f',num2str(m)]; 
        fid = fopen(toRun_sh,'w');

        cmd_combine = [handles.recon_path_server,cmd_combine,'\n\n']; 
       fprintf(fid, cmd_combine);

        for lmk = 1:8
          str_delete_lm = ['rm ', handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.',num2str(lmk),'.lm\n\n']; 
          fprintf(fid,str_delete_lm); 
          str_delete_mul_fac = ['rm ', handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.',num2str(lmk),'.mul_fac\n\n'];
          fprintf(fid,str_delete_mul_fac); 
          str_delete_add_fac = ['rm ', handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.',num2str(lmk),'.add_fac\n\n'];
          fprintf(fid, str_delete_add_fac); 
        end

        %str = ['export OMP_NUM_THREADS=',num2str(handles.num_threads_server),'\n\n'];
      
        %fprintf(fid,str);
        %str = '';
        %str = [handles.recon_path_server, 'lmrecon_tof  ',handles.server_recon_data_dir,'/', outfolder_server_temp,'/lmacc_scanner_parameter_f',num2str(m),'.cfg\n\n'];
     
        %fprintf(fid,str);
        str = ''; 

        fclose(fid);

        ssss = ['chmod +x ', toRun_sh]; 
        system(ssss); 
        pause(0.1);


        for ks = 1:handles.num_sub_frames
          
          make_explorer_lmacc_config_server(handles, ks, outfolder_server_temp);
          toRun_sh = [handles.lm_outfolder,'run_lmrecon_explorer_f',num2str(ks)]; 
          fids = fopen(toRun_sh,'w');

          str = ['export OMP_NUM_THREADS=',num2str(handles.num_threads_server),'\n\n'];
      
          fprintf(fids,str);
          str = '';
          str = [handles.recon_path_server, 'lmrecon_tof  ',handles.server_recon_data_dir,'/', outfolder_server_temp,'/lmacc_scanner_parameter_f',num2str(ks),'.cfg\n\n'];
     
          fprintf(fids,str);
          str = ''; 

          fclose(fids); 
          
          ssss = ['chmod +x ', toRun_sh]; 
          system(ssss); 
          pause(0.1);

        end

      else  % i.e. if subsample is NOT on
  	    % combine 8 .lm files into 1 for UCD recon
  	    cmd_combine = ['combine_listmode ',...
  		  handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.lm ',...
  		  handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.1.lm ',...
  		  handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.2.lm ',...
  		  handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.3.lm ',...
  		  handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.4.lm ',...
  		  handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.5.lm ',...
  		  handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.6.lm ',...
  		  handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.7.lm ',...
  		  handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.8.lm '];  

        make_explorer_lmacc_config_server(handles, m, outfolder_server_temp); 

        toRun_sh = [handles.lm_outfolder,'run_lmrecon_explorer_f',num2str(m)]; 
        fid = fopen(toRun_sh,'w');

        cmd_combine = [handles.recon_path_server,cmd_combine,'\n\n']; 
        fprintf(fid, cmd_combine);

        for lmk = 1:8
          str_delete_lm = ['rm ', handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.',num2str(lmk),'.lm\n\n']; 
        %  fprintf(fid,str_delete_lm); 
          str_delete_mul_fac_original = ['rm ', handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.',num2str(lmk),'.mul_fac.original\n\n'];
        %  fprintf(fid,str_delete_mul_fac_original); 
          str_delete_add_fac = ['rm ', handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.',num2str(lmk),'.add_fac\n\n'];
        %  fprintf(fid, str_delete_add_fac); 
          str_delete_add_fac_original = ['rm ', handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.',num2str(lmk),'.add_fac.original\n\n'];
        %  fprintf(fid, str_delete_add_fac_original);
          str_delete_attn_fac = ['rm ', handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.',num2str(lmk),'.attn_fac\n\n'];
        %  fprintf(fid, str_delete_attn_fac); 
          str_delete_checkfile = ['rm ', handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.',num2str(lmk),'.checkfile\n\n'];
        %  fprintf(fid, str_delete_checkfile); 
        end

        str = ['export OMP_NUM_THREADS=',num2str(handles.num_threads_server),'\n\n'];
      
        fprintf(fid,str);
        str = '';
        if use_scat_corr > 0.5
          str = ['cd ', handles.install_dir_server, '/scatter_correct','; matlab -nodesktop -nodisplay -r "scat_corr_recon(', num2str(m), ')" & \n pid=$! \necho $pid > process_id\n\n']; % run scatt corr recon 
        else
          str = [handles.recon_path_server, 'lmrecon_tof  ',handles.server_recon_data_dir,'/', outfolder_server_temp,'/lmacc_scanner_parameter_f',num2str(m),'.cfg\n\n'];
        end     
        fprintf(fid,str);
        str = ''; 

        fclose(fid);

        fprintf(fid_log, 'now creating run_lmrecon_explorer_fX-script (X=frame number), which does the following:\n1) executing combine_listmode\n2) deleting individual lm files\n3) executing recon commands\n');
        fprintf('now creating run_lmrecon_explorer_fX-script (X=frame number), which does the following:\n1) executing combine_listmode\n2) deleting individual lm files\n3) executing recon commands\n');
        ssss = ['chmod +x ', toRun_sh]; 
        system(ssss); 
        pause(0.1);

      end  % if condition (subsample ON / OFF)
        
        % move the data to recon server 
      fprintf(fid_log, 'moving data to recon server directory %s\n', outfolder_server_temp); fprintf('moving data to recon server directory %s\n', outfolder_server_temp);
  	  cmd_mvdata = ['mv ', '"',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.1.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.1.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.1.add_fac" "',...
  		handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.2.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.2.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.2.add_fac" "',...
  		handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.3.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.3.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.3.add_fac" "',...
  		handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.4.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.4.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.4.add_fac" "',...
  		handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.5.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.5.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.5.add_fac" "',...
  		handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.6.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.6.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.6.add_fac" "',...
  		handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.7.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.7.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.7.add_fac" "',...
  		handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.8.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.8.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.8.add_fac" ',...
      handles.server_recon_data_dir, '/',outfolder_server_temp,'/ '];  
      cmd_cpdata = ['cp -r ', '"', toRun_sh,'" "', handles.sensitivity_path_server_backup,'" "', handles.lm_outfolder,'lmacc_scanner_parameter_f',num2str(m),'.cfg" ', handles.server_recon_data_dir, '/',outfolder_server_temp,'/ '];
      system(cmd_mvdata); 
      pause(0.1);
      system(cmd_cpdata);
      pause(0.1);

      % get path to script to generate attenuation factors
      attn_path_name = [handles.install_dir_server, '/read_lm/generate_lm_attn_fp_exp'];
      p_attn = genpath(attn_path_name); 
      addpath(p_attn); 

      % run script that adds attenuation factors to add_fac files 
      fprintf(fid_log, 'now running attn_fp to create the attenuation factors for frame number %d\n', m);fprintf('now running attn_fp to create the attenuation for frame number %d\n', m);
      attn_fp(m); % m is frame number
      pause(0.1);
      fprintf(fid_log, 'done running attn_fp.\n');fprintf('done running attn_fp.\n');

      if subsample_on
        for ks = 0:handles.num_sub_frames
          toRun_sh_sub = [handles.lm_outfolder,'run_lmrecon_explorer_f',num2str(ks)]; 
          cmd_mvdata_sub = ['cp -r ','"',toRun_sh_sub,'" ', handles.server_recon_data_dir, '/',outfolder_server_temp,'/ '];
          system(cmd_mvdata_sub);
          pause(0.01); 

          config_fname_sub = [handles.lm_outfolder, 'lmacc_scanner_parameter_f', num2str(ks), '.cfg']; 
          cmd_mvcfg_sub  = ['cp -r ','"',config_fname_sub,'" ', handles.server_recon_data_dir, '/',outfolder_server_temp,'/ '];
          system(cmd_mvcfg_sub); 
          pause(0.01);

        end
      end
  	  
  	  
  	  sensitivity_path_server_temp1 = [handles.server_recon_data_dir, '/',outfolder_server_temp, handles.sensitivity_image_name];
  	  sensitivity_path_server_temp2 = [handles.server_recon_data_dir, '/',outfolder_server_temp, remove_space_name(handles.sensitivity_image_name)];
  	  
  	  
  	  cmd_rename = ['mv -T "', sensitivity_path_server_temp1, '" "',sensitivity_path_server_temp2,'"']; 
  	  system(cmd_rename); 
  	  pause(0.1); 

  	  idn_f = strfind(cmd_mvdata, '.8.add_fac'); 
  	  cmd_deletedata = ['rm ',cmd_mvdata(6:idn_f+10)]; 
  	  %cmd_deletedata = ['rm -r "', handles.outfolder, 'lm_reorder_f', num2str(m), '_prompts.lm ',handles.outfolder,'lm_reorder_f',num2str(m),'_prompts.add_fac" ']; 
  	 % system(cmd_deletedata); 
  	  pause(0.1); 

      lm_data_ready(lm_counter) = 1; 
      lm_counter = lm_counter + 1;
    end % if ch_lm
    
  end % if ~lm_done

  if sum(lm_data_ready(:)) == (length(lm_data_ready)-1)
    lm_done = true; 
    fprintf(fid_log, 'lm done\n'); fprintf('lm done\n');
  end


  % check how many images are completed
  img_ind_start = find(recon_frame_done == 0, 1, 'first'); 
  img_ind_end = min([(img_ind_start+3*handles.max_par_recon), length(handles.reconFrames)]); % max_par_recon = 12

  for kk = img_ind_start:img_ind_end
    fprintf(fid_log, 'frame %d done: %d\n', kk, double(recon_frame_done(kk)));%toBeDeleted
    if recon_frame_done(kk) < 0.5
      ch_img = true;
      for ii = handles.osem_iter:handles.osem_iter        
        if subsample_on
          fname_img_ch = [handles.server_recon_data_dir, '/', handles.server_temp_dir{1},'/',handles.fname_recon, '_f', num2str(handles.reconFrames(kk)), '.intermediate.',num2str(ii)];
          fname_img_ch_all = [handles.server_recon_data_dir, '/', handles.server_temp_dir{1},'/',handles.fname_recon, '_f', num2str(handles.reconFrames(kk)), '.intermediate.*'];
          fname_img_ch_local = [handles.lm_outfolder, handles.fname_recon, '_f', num2str(handles.reconFrames(kk)), '.intermediate.',num2str(ii)];
        else
          fname_img_ch = [handles.server_recon_data_dir, '/', handles.server_temp_dir{kk},'/',handles.fname_recon, '_f', num2str(handles.reconFrames(kk)), '.intermediate.',num2str(ii)];
          fname_img_ch_all = [handles.server_recon_data_dir, '/', handles.server_temp_dir{kk},'/',handles.fname_recon, '_f', num2str(handles.reconFrames(kk)), '.intermediate.*'];
          fname_img_ch_local = [handles.lm_outfolder, handles.fname_recon, '_f', num2str(handles.reconFrames(kk)), '.intermediate.',num2str(ii)];
        end
        
        cmd_imgcheck = ['test -f "',fname_img_ch,'" && cp -r ',fname_img_ch_all, ' "', handles.lm_outfolder,'"']; 
        
        system(cmd_imgcheck); 
        pause(0.01); 
        
        if ~exist(fname_img_ch_local, 'file')
          ch_img = false; 
        end
      end % for loop over iterations
      
      if ch_img
        if ~subsample_on
          fprintf(fid_log, 'cleaning up the lm_outfolder %s\n', handles.lm_outfolder);
          fprintf('cleaning up the lm_outfolder %s\n', handles.lm_outfolder);

          % wipe lm_outfolder
          if wipe_directories
            cmd_rm = ['rm ', handles.lm_outfolder, 'block_sino_f', num2str(handles.reconFrames(kk)),'_prompts.* ',...
            handles.lm_outfolder, 'block_sino_f', num2str(handles.reconFrames(kk)), '_randoms.* ',...
            handles.lm_outfolder, 'deadtime_sino.* ',...
            handles.lm_outfolder, 'lm_reorder_f', num2str(handles.reconFrames(kk)), '_prompts.* ',...
            handles.lm_outfolder, 'prompts* '...
            handles.lm_outfolder, 'randoms* '...
            handles.lm_outfolder, 'singles* '];
            fprintf(fid_log, '%s\n', cmd_rm); fprintf(cmd_rm);
            system(cmd_rm);
            pause(0.1);
          end
          %copy lm and correction files from temp dir to lm_outfolder
        %  cmd_mv = ['mv ', handles.server_recon_data_dir, '/', handles.server_temp_dir{kk}, '/*prompts.add_fac ', handles.lm_outfolder ];
        %  cmd_mv2 = ['mv ', handles.server_recon_data_dir, '/', handles.server_temp_dir{kk}, '/*prompts.add_fac.original ', handles.lm_outfolder ];
        %%  cmd_mv3 = ['mv ', handles.server_recon_data_dir, '/', handles.server_temp_dir{kk}, '/*prompts.mul_fac ', handles.lm_outfolder ];
        %  cmd_mv4 = ['mv ', handles.server_recon_data_dir, '/', handles.server_temp_dir{kk}, '/*prompts.lm ', handles.lm_outfolder ];
        %  cmd_mv5 = ['mv ', handles.server_recon_data_dir, '/', handles.server_temp_dir{kk}, '/*prompts.scat_fac ', handles.lm_outfolder ];
        %  cmd_mv6 = ['mv ', handles.server_recon_data_dir, '/', handles.server_temp_dir{kk}, '/*prompts.attn_fac ', handles.lm_outfolder ];
        %  cmd_mv7 = ['mv ', handles.server_recon_data_dir, '/', handles.server_temp_dir{kk}, '/*prompts.mul_fac.original ', handles.lm_outfolder ];
        %  fprintf(fid_log, 'copying lm files to lm_outfolder:\n%s\n', cmd_mv); fprintf('copying lm files to lm_outfolder');
        %  system(cmd_mv); pause(0.1);
        %  fprintf(fid_log, 'copying lm files to lm_outfolder:\n%s\n', cmd_mv2); fprintf('copying lm files to lm_outfolder:\n%s\n', cmd_mv2);
        %  system(cmd_mv2); pause(0.1);
        %%  fprintf(fid_log, 'copying lm files to lm_outfolder:\n%s\n', cmd_mv3); fprintf('copying lm files to lm_outfolder:\n%s\n', cmd_mv3);
        %%  system(cmd_mv3); pause(0.1);
        %  fprintf(fid_log, 'copying lm files to lm_outfolder:\n%s\n', cmd_mv4); fprintf('copying lm files to lm_outfolder:\n%s\n', cmd_mv4);
        %  system(cmd_mv4); pause(0.1); 
        %  if use_scat_corr
        %    fprintf(fid_log, 'copying lm files to lm_outfolder:\n%s\n', cmd_mv5); fprintf('copying lm files to lm_outfolder:\n%s\n', cmd_mv5);
        %    system(cmd_mv5); pause(0.1);
        %  end
        %  fprintf(fid_log, 'copying lm files to lm_outfolder:\n%s\n', cmd_mv6); fprintf('copying lm files to lm_outfolder:\n%s\n', cmd_mv6);
        %  system(cmd_mv6); pause(0.1);
        %  fprintf(fid_log, 'copying lm files to lm_outfolder:\n%s\n', cmd_mv7); fprintf('copying lm files to lm_outfolder:\n%s\n', cmd_mv7);
        %  system(cmd_mv7); pause(0.1); 

          % wipe the temporary dir on the server
          if wipe_directories
            cmd_clean = ['cd ', handles.server_recon_data_dir,' && rm -r ', handles.server_temp_dir{kk},'/ ',]; 
            system(cmd_clean); 
            pause(0.01);
          end
        end
        recon_frame_done(kk) = 1; 
        recon_frame_running(kk) = 0; 
      end % if ch_img
    end % end of if recon_frame_done(kk) < 0.5
  end % for loop over images


  if sum(recon_frame_done(:)) == length(handles.reconFrames)
    not_done = false; 
    fprintf(fid_log, 'recon frames done\n'); fprintf('recon frames done\n');
  end
    
 
  for N = 1:length(handles.reconFrames)
    num_recons = sum(recon_frame_running(:));
    N_lm = N; 
    if subsample_on
      N_lm = 1; 
    end
    if (recon_frame_running(N) == 0) && (lm_data_ready(N_lm) == 1) && (recon_frame_done(N) == 0) && (num_recons < handles.max_par_recon) && not_done
      % start recon frame xxxx
      %run_uex_recon_x(handles, handles.reconFrames(N));
      
      % run the recon script
      recon_frame_running(N) = 1; 
      
      if subsample_on
            if N==1
              cmd_runrecon = [handles.server_recon_data_dir,'/',handles.server_temp_dir{1},'/run_lmrecon_explorer_f',num2str(handles.reconFrames(N)),' & wait '];
              recon_frame_done(N) = 1;
              recon_frame_running(N) = 0;   
            else
             cmd_runrecon = [handles.server_recon_data_dir, '/',handles.server_temp_dir{1},'/run_lmrecon_explorer_f',num2str(handles.reconFrames(N)),' & ']; 
            end 
 
      else
            cmd_runrecon = [handles.server_recon_data_dir,'/',handles.server_temp_dir{N},'/run_lmrecon_explorer_f',num2str(handles.reconFrames(N)),' & ']; 
      end
	  %cmd_runrecon = ['ssh ', handles.user, '@', handles.server_name, ' "cd ', handles.recon_path_server,' ; ',cmd_combine,' ; cd ', handles.server_recon_data_dir,'/',outfolder_server_temp,' ; ./run_lmrecon_explorer_f',num2str(m),'"; ', cmd_imgreturn, ' & ']; 
      fprintf(fid_log, 'starting recon for frame %d\n', handles.reconFrames(N)); fprintf('starting recon for frame %d\n', handles.reconFrames(N));
      system(cmd_runrecon); 
      pause(0.1); 
     
      %server_temp_dir{N} = fout;
    end
  end
  

  % pause for 5 seconds, then look again.
  pause(5); 


end


fprintf(fid_log, 'done all recons of all frames and iterations.'); %toBeDeleted
disp('done all recons of all frames and iterations.'); %toBeDeleted

fclose(fid_log); %toBeDeleted
%delete(handles_name); 
pause(0.1); 
quit











function make_explorer_lmacc_config_server(handles, frame, server_dir)

use_scat_corr = handles.scatter_onoff.Value;

fname_cfg = [handles.lm_outfolder, 'lmacc_scanner_parameter_f', num2str(frame), '.cfg']; 

fid1 = fopen(fname_cfg,'w'); 

str = 'detector_ring_diameter = 786.0\n\n';
fprintf(fid1,str);
str = ''; 


str = 'crystal_size = 2.85, 2.85, 2.85, 18.1\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'crystal_gap_size = 0.0, 0.0, 0.0\n\n'; 
fprintf(fid1,str); 
str = ''; 

 
str = 'crystal_array_size = 35, 679\n\n'; 
fprintf(fid1,str); 
str = ''; 

 
str = 'number_of_detector_modules = 24, 1\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'TOF_information = 505, 39.0625\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'number_of_radial_bins = 549\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'image_size = 239, 239, 679\n\n';
fprintf(fid1,str); 
str = ''; 


str = 'voxel_size = 2.85, 2.85, 2.85\n\n'; 
fprintf(fid1,str); 
str = ''; 

%sensitivity_path_server_temp = [handles.server_recon_data_dir, '/',server_dir, handles.sensitivity_image_name]; 
sensitivity_path_server_temp = [handles.server_recon_data_dir, '/',server_dir, remove_space_name(handles.sensitivity_image_name)];
%sensitivity_path_server_temp = handles.sensitivity_path_server_backup; 
%str = ['sensitivity = ',handles.sensitivity_path_server, '\n\n'];
str = ['sensitivity = ',sensitivity_path_server_temp, '\n\n'];  
fprintf(fid1,str); 
str = ''; 


%str = 'iPSF_model = /run/media/meduser/data/lmrecon_explorer/aker_679x679_rd339.pmat\n\n'; 


str = ['iPSF_model = ',handles.tker_path_server, ', ', handles.aker_path_server, '\n\n']; %/run/media/meduser/data/lmrecon_explorer/tker_57121x57121.pmat, /run/media/meduser/data/lmrecon_explorer/aker_679x679_rd339.pmat\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'iterative_algorithm_type = 0\n\n';
fprintf(fid1,str); 
str = ''; 

str = 'warmup_setting = -1, 1\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_strength = 1e-3\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_model_type = 0\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_potential_function_type = 2   #0: Quadratic, 1: Hyperbola, 2: Fair, 3: Huber, 4: hw\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_neighborhood_properties = 1, 1    #size, isotropic  #0: 3x3x3, 1:isotropic  # 1st(0) or 2nd(1), aniso(0) or iso(1)\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_buildin_parameter_list = 1e-9\n\n'; 
fprintf(fid1,str); 
str = '';

if use_scat_corr > 0.5
  str = ['iteration_setting = 13, 1, 1    #  subsets, save step, iteration number\n\n']; 
else
  str = ['iteration_setting = 13, 1, ',num2str(handles.osem_iter),'    #  subsets, save step, iteration number\n\n']; 
end
fprintf(fid1,str); 
str = '';

if use_scat_corr > 0.5
  str = ['initial_guess = ', handles.server_recon_data_dir, '/', server_dir, '/img_guess_next_iter_f', num2str(frame), '\n\n'];
  fprintf(fid1,str); 
  str = '';
end

str = ['input_raw_data_file = ',handles.server_recon_data_dir,'/',server_dir,'/lm_reorder_f',num2str(frame),'_prompts\n\n']; 
fprintf(fid1,str); 
str = ''; 



str = 'input_raw_data_format_type = 0\n\n'; 
fprintf(fid1,str); 
str = ''; 

if use_scat_corr > 0.5
  str = ['reconstruction_output_setting = ',handles.server_recon_data_dir,'/',server_dir,', ./lmrecon_output_f',num2str(frame),'\n\n']; 
else
  str = ['reconstruction_output_setting = ',handles.server_recon_data_dir,'/',server_dir,', ./',handles.fname_recon,'_f',num2str(frame),'\n\n']; 
end
fprintf(fid1,str); 


fclose(fid1); 



function name_nospace = remove_space_name(name_in)

name_nospace = name_in; 
for k = 1:length(name_in)
  if strcmp(name_in(k), ' ')
    name_nospace(k) = '_'; 
  end
end


