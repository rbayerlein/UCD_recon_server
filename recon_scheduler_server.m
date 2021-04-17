function recon_scheduler_server

handles_name = '/home/rbayerlein/code/explorer-master/handles_scheduler.mat'; 

pause(0.1); 
h = load(handles_name); 
handles = h.handles
clear h; 
pause(0.1); 

delete(handles_name); 
pause(0.1); 

lm_counter = 1; 
img_counter = 1; 

lm_data_ready = zeros(length(handles.reconFrames)+1, 1); 

subsample_on = handles.subsample; 
subsample_done = false; 

if subsample_on
  reconFrames_lm = [0,1]; 
else
  reconFrames_lm = [handles.reconFrames, (handles.reconFrames(end)+1)];
end

lm_data_ready = zeros(length(reconFrames_lm), 1); 


recon_frame_done = zeros(length(handles.reconFrames), 1); 
recon_frame_running = zeros(length(handles.reconFrames), 1); 





not_done = true; 
recon_full = false; 
lm_done = false; 

while not_done
  
  % first check if the next in line lm files are all ready, by checking if all next ones are created. 
  if ~lm_done
    ch_lm = true; 
    for k = 1:8
      fname_lm_ch = [handles.lm_outfolder, 'lm_reorder_f', num2str(reconFrames_lm(lm_counter+1)), '_prompts.', num2str(k), '.lm']; 
      if ~exist(fname_lm_ch, 'file');
        ch_lm = false; 
      end
    end
    if ch_lm

      m = reconFrames_lm(lm_counter);
  	  %andles.server_temp_dir{lm_counter} = outfolder_server_temp; 
      outfolder_server_temp = handles.server_temp_dir{lm_counter}; 

  	  % delete data in the recon data drive on the server
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

      else
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
          fprintf(fid,str_delete_lm); 
          str_delete_mul_fac = ['rm ', handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.',num2str(lmk),'.mul_fac\n\n'];
          fprintf(fid,str_delete_mul_fac); 
          str_delete_add_fac = ['rm ', handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f',num2str(m),'_prompts.',num2str(lmk),'.add_fac\n\n'];
          fprintf(fid, str_delete_add_fac); 
        end

        str = ['export OMP_NUM_THREADS=',num2str(handles.num_threads_server),'\n\n'];
      
        fprintf(fid,str);
        str = '';
        str = [handles.recon_path_server, 'lmrecon_tof  ',handles.server_recon_data_dir,'/', outfolder_server_temp,'/lmacc_scanner_parameter_f',num2str(m),'.cfg\n\n'];
     
        fprintf(fid,str);
        str = ''; 

        fclose(fid);

        ssss = ['chmod +x ', toRun_sh]; 
        system(ssss); 
        pause(0.1);

      end  
        
        % move the data to recon server 
  	  cmd_mvdata = ['cp -r ', '"',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.1.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.1.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.1.add_fac" "',...
  		handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.2.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.2.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.2.add_fac" "',...
  		handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.3.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.3.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.3.add_fac" "',...
  		handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.4.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.4.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.4.add_fac" "',...
  		handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.5.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.5.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.5.add_fac" "',...
  		handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.6.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.6.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.6.add_fac" "',...
  		handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.7.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.7.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.7.add_fac" "',...
  		handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.8.lm" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.8.mul_fac" "',handles.lm_outfolder,'lm_reorder_f',num2str(m),'_prompts.8.add_fac" "',...
  		handles.sensitivity_path_server_backup,'" "', toRun_sh,'" "', handles.lm_outfolder,'lmacc_scanner_parameter_f',num2str(m),'.cfg" ', handles.server_recon_data_dir, '/',outfolder_server_temp,'/ '];  
  	  system(cmd_mvdata); 
  	  pause(0.1);

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
  	  system(cmd_deletedata); 
  	  pause(0.1); 

      lm_data_ready(lm_counter) = 1; 
      lm_counter = lm_counter + 1;
    end
    
  end

  if sum(lm_data_ready(:)) == (length(lm_data_ready)-1)
    lm_done = true; 
  end


  % check how many images are completed
  img_ind_start = find(recon_frame_done == 0, 1, 'first'); 
  img_ind_end = min([(img_ind_start+3*handles.max_par_recon), length(handles.reconFrames)]); 
   
  for kk = img_ind_start:img_ind_end
    if recon_frame_done(kk) < 0.5
      ch_img = true;
      for ii = handles.osem_iter:handles.osem_iter

        fname_img_ch_local = [handles.lm_outfolder, handles.fname_recon, '_f', num2str(handles.reconFrames(kk)), '.intermediate.',num2str(ii)]; 
        
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
      end
      
      if ch_img
        if ~subsample_on
          % wipe the temporary dir on the server
          cmd_clean = ['cd ', handles.server_recon_data_dir,' && rm -r ', handles.server_temp_dir{kk},'/ ',]; 
          system(cmd_clean); 
          pause(0.01);
        end
        recon_frame_done(kk) = 1; 
        recon_frame_running(kk) = 0; 
        img_counter = img_counter + 1; 
      end
    end
  end

  num_recons = sum(recon_frame_running(:));

  if sum(recon_frame_done(:)) == length(handles.reconFrames)
    not_done = false; 
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
              subsample_done = true;
              recon_frame_done(N) = 1;
              recon_frame_running(N) = 0;   
            else
             cmd_runrecon = [handles.server_recon_data_dir, '/',handles.server_temp_dir{1},'/run_lmrecon_explorer_f',num2str(handles.reconFrames(N)),' & ']; 
            end 
 
          else
            cmd_runrecon = [handles.server_recon_data_dir,'/',handles.server_temp_dir{N},'/run_lmrecon_explorer_f',num2str(handles.reconFrames(N)),' & ']; 
          end
	  %cmd_runrecon = ['ssh ', handles.user, '@', handles.server_name, ' "cd ', handles.recon_path_server,' ; ',cmd_combine,' ; cd ', handles.server_recon_data_dir,'/',outfolder_server_temp,' ; ./run_lmrecon_explorer_f',num2str(m),'"; ', cmd_imgreturn, ' & ']; 

	  system(cmd_runrecon); 
	  pause(0.1); 
     
      %server_temp_dir{N} = fout;
    end
  end
  

  % pause for 5 seconds, then look again.
  pause(5); 


end


quit











function make_explorer_lmacc_config_server(handles, frame, server_dir)

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


str = 'TOF_information = 460, 39.0625\n\n'; 
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


str = ['iteration_setting = 13, 1, ',num2str(handles.osem_iter),'    #  subsets, save step, iteration number\n\n']; 
fprintf(fid1,str); 
str = '';

%     #  KEM   2

%# initial_guess = ./lmrecon_tof_475x475x1355_200MBq_20min_PSF.intermediate.1


str = ['input_raw_data_file = ',handles.server_recon_data_dir,'/',server_dir,'/lm_reorder_f',num2str(frame),'_prompts\n\n']; 
fprintf(fid1,str); 
str = ''; 



str = 'input_raw_data_format_type = 0\n\n'; 
fprintf(fid1,str); 
str = ''; 

str = ['reconstruction_output_setting = ',handles.server_recon_data_dir,'/',server_dir,', ./',handles.fname_recon,'_f',num2str(frame),'\n\n']; 
fprintf(fid1,str); 

fclose(fid1); 



function name_nospace = remove_space_name(name_in)

name_nospace = name_in; 
for k = 1:length(name_in)
  if strcmp(name_in(k), ' ')
    name_nospace(k) = '_'; 
  end
end


