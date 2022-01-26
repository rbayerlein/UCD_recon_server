function attn_fp(frame_number_from_recon_scheduler_server)
% function to perform forward projection (fp) of attenuation (attn) 
% coefficients for correct normalization of randoms and scatter corrections
% rbayerlein, 10/29/2021

%% HARD CODED PARAMETERS 
% this is for parameters like frame number that might be relevant to be
% used later on. So the infrastructure is created here.
frame_num = frame_number_from_recon_scheduler_server;  % can get this from handler later

num_threads = 48; % for forward projector

mumap_fname = 'CTAC_201_mumap_kVp-140_size-239x239x679_vox-2.85x2.85x2.85.img';

lm_fp_exp = '/home/rbayerlein/code/explorer-master/read_lm/lmrecon_exploer/app/lm_fp_exp';

AddAttn2add_fac = '/home/rbayerlein/code/explorer-master/read_lm/AddAttn2add_fac/bin/AddAttn2add_fac';

%% read in handler and parameters
handles_name = '/home/rbayerlein/code/explorer-master/handles_scheduler.mat'; 
pause(0.1); 
h = load(handles_name); 
handles = h.handles;
clear h; 
pause(0.1); 

handles.server_temp_dir

outfolder_server_temp = handles.server_temp_dir{frame_num+1};

%% check if lm files exist already
lm_not_done = true;
timeStamps = ones(1,8);
timeStamps=timeStamps*(-1);
file_exists = zeros(1,8);
file_finished = zeros(1,8);

while lm_not_done
    for i=1:8
        file2check = [handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f', num2str(frame_num), '_prompts.', num2str(i), '.lm'];
        fprintf('now checking file %d\n\t%s\n', i, file2check);
        if file_exists(1,i) == 0 % file not registered as 'existing'
            if exist(file2check, 'file') % check if exists
                file_exists(1,i) = 1;   % register as existing
                FileInfo = dir(file2check); %get time stamp
                timeStamps(1,i) = FileInfo.datenum; % save time stamp
            end
        else    % if file has been registered as 'existing'
            % compare timeStamps
            FileInfo = dir(file2check); % get time stamp
            currentTimeStamp = FileInfo.datenum; 
            if currentTimeStamp == timeStamps(1,i) % compare time stamp
                file_finished(1,i) = 1;
                fprintf('lm file %d done.', i);
            else
                timeStamps(1,i) = currentTimeStamp;
            end
        end
    end
    if sum(file_finished) == 8
        disp('all lm files done');
        lm_not_done = false;
    else
        pause(5.0);
    end
end
        
%% create config files for attn fp
cfg_folder = [handles.server_recon_data_dir, '/', outfolder_server_temp, '/cfg_attn_fp_exp'];
if ~exist(cfg_folder, 'dir')
    mkdir(cfg_folder);
end
cfg_fname = [cfg_folder, '/copy_cfg_lmfp.sh'];

% create script that writes the config files
CreateConfigFiles(cfg_fname, handles, outfolder_server_temp, frame_num, mumap_fname);
pause(0.1);
% run that script
cmd = ['bash ', cfg_fname];
system(cmd);
pause(0.1);

%% perform attenuation forward projection

for i = 1:8
    file2BeCreated = [handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f', num2str(frame_num), '_prompts.', num2str(i), '.attn_fac'];
    cmd = [lm_fp_exp, ' ', cfg_folder, '/lmrecon_attn_fp_exp_f', num2str(i), '.cfg ', num2str(num_threads)];
    if exist(file2BeCreated, 'file')
        fprintf('file %d (%s) exists already. Omitting.\n', i, file2BeCreated);
    else
        fprintf('now running forward projection for file %d: \n%s\n', i, cmd);
        system(cmd);
    end
end

files_done = zeros(1,8);
not_done = true;
while not_done
    for i = 1:8
        file2check = [handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f', num2str(frame_num), '_prompts.', num2str(i), '.attn_fac'];
        fprintf('\nchecking file %d...', i);
        if exist(file2check, 'file')
            files_done(1,i) = 1;
            disp('exists');
        else
            cmd_tch = ['touch ', file2check];
            system(cmd_tch);
            disp('will be created');
        end
    end
    if sum(files_done) == 8
        not_done = false;
        disp('all attn files created');
    else
        fprintf('\nchecking again ...\n');
        pause(2.0);
    end
end


%% run AddAttn2add_fac in parallel

crys_eff = [handles.dcm_dir_init_ucd_server, '/crys_eff_679x840'];
plane_eff = [handles.dcm_dir_init_ucd_server, '/plane_eff_679x679'];

for i = 1:8
    cmd = [AddAttn2add_fac, ' ', handles.server_recon_data_dir,'/',outfolder_server_temp,'/lm_reorder_f', num2str(frame_num), '_prompts.', num2str(i), '.lm ',...
        handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f', num2str(frame_num), '_prompts.', num2str(i), '.attn_fac ', crys_eff, ' ', plane_eff, '&'];
    fprintf('now running command %s\n', cmd); 
    system(cmd);
    pause(0.2);
end

% check if add fac files have been updated: look for a check file created
% by AddAttn2add_fac
cf_exists = zeros(1,8);
%check if file exists (i.e. process has started)
not_done = true;
while not_done
    for i = 1:8
         file2check = [handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f', num2str(frame_num), '_prompts.', num2str(i), '.checkfile'];
         if exist(file2check, 'file')
             cf_exists(1,i) = 1;
         end
    end
    if sum(cf_exists) == 8
        not_done = false;
        fprintf('all checkfiles created');
    end
    pause(0.01);
end

%now check, if it has been deleted again (i.e. process is finished)
files_done = zeros(1,8);
not_done = true;
while not_done
     for i = 1:8
        % check if temp file has been deleted indicating the process is done
        file2check = [handles.server_recon_data_dir,'/',outfolder_server_temp, '/lm_reorder_f', num2str(frame_num), '_prompts.', num2str(i), '.checkfile'];
        fprintf('\nchecking file %d...', i);
        if ~exist(file2check, 'file')
            files_done(1,i) = 1;
            disp('checkfile deleted. done.');
        else
            disp('not done.');
        end
    end
    if sum(files_done) == 8
        not_done = false;
        disp('all attn files created');
    else
        pause(2.0);
        fprintf('\nchecking again ...\n');
    end   
end

fprintf('done with attenuation forward projection.\n');
pause(0.1);
end




%% function that creates the bash cript which in turn creates the 8 different config files for the attenuation forward projection
function CreateConfigFiles(cfg_fname, handles, outfolder_server_temp, frame_num, mumap_fname)

fid_cfg = fopen(cfg_fname, 'w');
fprintf(fid_cfg, '#!/bin/bash \n');
fprintf(fid_cfg, 'study_dir=%s\n', handles.study_dir_server);
fprintf(fid_cfg, 'server_temp_dir=%s/%s\n', handles.server_recon_data_dir,outfolder_server_temp);
fprintf(fid_cfg, 'cfg_folder=${server_temp_dir}/cfg_attn_fp_exp \n');
fprintf(fid_cfg, 'lm_folder=${server_temp_dir} \n');
fprintf(fid_cfg, 'if [[ ! -d $cfg_folder ]]; then \n	mkdir $cfg_folder\n	chmod -R 775 $cfg_folder\nfi\n');
fprintf(fid_cfg, 'if [[ ! -d $lm_folder ]]; then \n	mkdir $lm_folder\n	chmod -R 775 $lm_folder\nfi\n');

fprintf(fid_cfg, 'for (( m=1; m<=8; m++ )) do \n{\n	echo "');
fprintf(fid_cfg, 'detector_ring_diameter = 786.0 \n');
fprintf(fid_cfg, 'crystal_size = 2.85, 2.85, 2.85, 18.1\n');
fprintf(fid_cfg, 'crystal_gap_size = 0.0, 0.0, 0.0 \n');
fprintf(fid_cfg, 'crystal_array_size = 35, 679  \n');
fprintf(fid_cfg, 'number_of_detector_modules = 24, 1 \n');
fprintf(fid_cfg, 'number_of_radial_bins = 549 \n');
fprintf(fid_cfg, 'voxel_size = %0.2f, %0.2f, %0.2f\n', handles.pet_voxel_size);
fprintf(fid_cfg, 'image_size = %d, %d, %d \n', handles.pet_img_size);

fprintf(fid_cfg, 'iterative_algorithm_type = 0 \n');

fprintf(fid_cfg, 'initial_guess =  ${study_dir}/UCD/Image/%s\n', mumap_fname);
fprintf(fid_cfg, 'warmup_setting = 1, 1 \n');
fprintf(fid_cfg, 'iteration_setting = 1, 1, 1    #  1, 10, 500 \n');
fprintf(fid_cfg, 'regularizer_strength = 0 # mandatory line, otherwise core dumps \n');
fprintf(fid_cfg, 'regularizer_model_type = 0    #0: pairwise MRF, 1: patch \n');
fprintf(fid_cfg, 'regularizer_potential_function_type = 2   #0: Quadratic, 1: Hyperbola, 2: Fair, 3: Huber, 4: hw \n');
fprintf(fid_cfg, 'regularizer_neighborhood_properties = 0, 1    #size, isotropic  #0: 3x3x3, 1:isotropic  # 1st(0) or 2nd(1), aniso(0) or iso(1) \n');
fprintf(fid_cfg, 'regularizer_buildin_parameter_list = 1e-10 \n');

fprintf(fid_cfg, 'input_raw_data_file = ${server_temp_dir}/lm_reorder_f%d_prompts.$[$m].lm\n', frame_num);

fprintf(fid_cfg, 'input_raw_data_format_type = 0 \n');

fprintf(fid_cfg, 'reconstruction_output_setting = ${lm_folder}, ./lm_reorder_f%d_prompts.$[$m].attn_fac\n', frame_num);

fprintf(fid_cfg, '	 " > ${cfg_folder}/lmrecon_attn_fp_exp_f$m.cfg \n');

fprintf(fid_cfg, '}\ndone \n');

fclose(fid_cfg);
end

