function scat_corr_recon(frame_num)
%% get access to handles
disp('reading in parameters from handler');
handles_name = '/home/rbayerlein/code/explorer-master/handles_scheduler.mat'; 
%handles_name = '/home/rbayerlein/Code/Recon/server/explorer-master/handles_scheduler.mat'; 
pause(0.1); 
h = load(handles_name);
handles = h.handles;
clear h; 
pause(0.1); 

%% set up log file
fname_log = [handles.lm_outfolder, 'scatter_recon.log'];
if exist(fname_log, 'file')
    delete(fname_log)
end
fid_log = fopen(fname_log, 'w');
cmd_chmod = ['chmod 775 ', fname_log]; 
system(cmd_chmod); 
pause(0.1);
fprintf(fid_log, 'done reading in parameters from handler\n'); 
fprintf(fid_log, 'starting scatter corrected list-mode reconstruction\n'); disp('starting scatter corrected list-mode reconstruction');

%% input parameters
use_scat_corr=handles.scatter_onoff.Value;
fprintf(fid_log, 'correct for scatter: %d (should be 1, otherwise this script should not have been invoked!)\n', double(use_scat_corr));

user = 'rbayerlein';
server_name = 'exp-sim-001';

%simset install dir
install_dir = '/home/rbayerlein/code/Recon/scatter_correction/';    % indepentent of handler file
amap_processing_dir = [install_dir, 'amap_processing'],
umap_processing_dir = [install_dir, 'umap_processing'];
scripts_dir = [install_dir, 'scripts'];
data_processing_dir = [install_dir, 'data_processing'];

fprintf(fid_log, 'frame number: %d\n', frame_num);

% u-map
dir_ctac = handles.AC_path_server;
fprintf(fid_log, 'sens image folder path: %s\n', dir_ctac);
tx_fov = 500; % mm

% a-map
outfolder_server_temp = handles.server_temp_dir{frame_num+1};
dir_amaps_general = [handles.server_recon_data_dir, '/', outfolder_server_temp, '/', handles.fname_recon, '_f', num2str(frame_num), '.intermediate.']; 
fprintf(fid_log, 'activity image folder path: %sY (Y = iteration number)\n', dir_amaps_general);

% amap image dimensions (CT will be downscaled to match amap)
img_size = [239,239,679];
vox_size = [2.85,2.85,2.85]; % mm

bin_phg = [install_dir, 'simset/2.9.2/bin/phg'];
bin_hist2lm = [install_dir, 'simset/data_processing/hist2lm/bin/hist2lm'];
bin_lm2blocksino =[install_dir, 'simset/data_processing/lm2blocksino_old/bin/lm2blocksino'];
bin_scatter_add_fac =[install_dir, 'simset/data_processing/scatter_add_fac/bin/scatter_add_fac'];

% LUT
fname_lut = [bin_lm2blocksino(1:strfind(bin_lm2blocksino,'/bin')), 'include/index_blockpairs_transaxial_2x91x60_int16'];

% u-map
hu_limits = [-1000,2000]; % lower and upper HU limits
fname_water = 'u_water_511keV.dat'; % 511 keV water density LUT

% a-map
n_bins = 22100; % number of available bins - same as number of materials

num_threads_avail = 44;    % actually it's 48 but should leave a few empty to prevent lagging. 
num_threads_per_frame = num_threads_avail;



%% main program


for iter = 1 : handles.osem_iter
%% recon part
    
    %start recon
    cmd_recon = [handles.recon_path_server, 'lmrecon_tof  ',handles.server_recon_data_dir,'/', outfolder_server_temp,'/lmacc_scanner_parameter_f',num2str(frame_num),'.cfg &'];
    msg = ['starting recon for iteration ', num2str(iter)],
    fprintf(fid_log, '%s\n', msg); disp(msg);
    system(cmd_recon);
    pause(0.1);
    
    % check if done and continue
    recon_file_name = [handles.server_recon_data_dir,'/', outfolder_server_temp, '/lmrecon_output_f', num2str(frame_num), '.intermediate.1'];
    msg = ['checking if done: %s', recon_file_name];
    fprintf(fid_log, '%s\n', msg); disp(msg);
    
    not_done = true;
    waitcounter=0; % just for control
    waittime = 5.0;
    while not_done  
        msg = ['Estimated time elapsed since start: ', num2str(waitcounter*waittime)];
        fprintf(fid_log, '%s\n', msg); disp(msg);
        
        if exist(recon_file_name, 'file')
            not_done = false;
            fprintf(fid_log, 'recon file exists. continue with next step...\n');
            disp('recon file exists. continue with next step...');
        end
        waitcounter = waitcounter+1;
        pause(waittime);    % wait 5 seconds then look again
    end
        
    % copy recon output files 
    fprintf(fid_log, 'copying recon output files\n'); disp('copying recon output files');
    % 1)
    cmd_cp_1 = ['cp ', recon_file_name, ' ', handles.server_recon_data_dir,'/', outfolder_server_temp, '/', handles.fname_recon, '_f', num2str(frame_num), '.intermediate.', num2str(iter)];
    system(cmd_cp_1);
    pause(0.1);
    % 2)    
    cmd_cp_2 = ['cp ', recon_file_name, ' ', handles.server_recon_data_dir,'/', outfolder_server_temp, '/img_guess_next_iter_f', num2str(frame_num)];
    system(cmd_cp_2);
    pause(0.1);
    
    %remove recon output files before starting new recon
    fprintf(fid_log, 'removing recon output file\n'); disp('removing recon output file');
    % .intermediate file
    cmd_rm = ['rm ', recon_file_name];
    system(cmd_rm);
    pause(0.1);
    %.img file
    recon_file_name_img = [recon_file_name(1:strfind(recon_file_name, '.intermediate.')), 'img' ]
    cmd_rm_img = ['rm ', recon_file_name_img];
    system(cmd_rm_img);
    pause(0.1);   
    
%% simulation part
    
    % skip scatter correction part if it's the last iteration
    if iter == handles.osem_iter
        fprintf(fid_log, 'finished recon of last iteration. no scatter correction neccessary. deleting image guess.');
        cmd_rm = ['rm ', handles.server_recon_data_dir,'/', outfolder_server_temp, '/img_guess_next_iter_f', num2str(frame_num)];
        system(cmd_rm);
        pause(0.1);
        break; % break out of big for loop
    end

    
    %base_dir for simulation output
    sim_outfolder_temp = handles.lm_outfolder;
    slash = strfind(sim_outfolder_temp, '/');
    if slash(end) == length(sim_outfolder_temp)
        slash_pos=slash(end-1);
        sim_outfolder_temp = sim_outfolder_temp(1:slash(end)-1);
    else
        slash_pos = slash(end);
    end
    dir_base = ['/mnt/data/rbayerlein/simset', sim_outfolder_temp(slash_pos:length(sim_outfolder_temp)), '_it', num2str(iter)];
    fprintf(fid_log, 'dir_base: %s\n', dir_base); disp(dir_base);

    % create main folder structures
    if not(isfolder(dir_base))
        CreateFolder(dir_base);
    end
    % create on sim node
    cmd_0 = ['ssh ', user, '@', server_name, ' "mkdir ', dir_base,'"']; 
    system(cmd_0);
    pause(0.1);
    fprintf(fid_log, 'created folder on sim node: \n%s\n', cmd_0);
    
    dir_umap = strcat(dir_base,'/umap');
    dir_amaps_out = strcat(dir_base,'/amap');
    dir_phg_params = strcat(dir_base,'/phg_params');
    dir_hist = strcat(dir_base,'/history');
    dir_phg_log = strcat(dir_base,'/phg_log');

    CreateFolder(dir_umap);
    CreateFolder(dir_amaps_out);
    CreateFolder(dir_phg_params);
    CreateFolder(dir_hist);
    CreateFolder(dir_phg_log);

    fprintf(fid_log, 'Created simset output sub-folders\n');
    
    % calculate number of events to simulate
    lm_filename = [handles.server_recon_data_dir, '/', outfolder_server_temp, '/lm_reorder_f',num2str(frame_num),'_prompts.lm'];
    if ~exist(lm_filename, 'file')
        warn_msg=['WARNING! lm file %s does not exist. This script will not finish.', lm_filename];
        disp(msg);
        fprintf(fid_log, '%s\n', warn_msg);
    end
    s = dir(lm_filename)
    num_to_simulate = double(s.bytes);   % simulate 10 times as may as in lm file. an event has 10 byte. So file size equals num_to_simulate.
    if num_to_simulate < 1e8
        num_to_simulate = 1e8;
    elseif num_to_simulate > 5e10
        num_to_simulate = 5e10;
    end
    fprintf(fid_log, 'lm file size (%s): %1.1f bytes\n', lm_filename, s.bytes);
    fprintf(fid_log, 'number of events to simulate: %d\n', num_to_simulate);

    % a-map for this run
    dir_amaps_in = [dir_amaps_general, num2str(iter)];
    fprintf(fid_log, '%s\n', dir_amaps_in); disp(dir_amaps_in);
    
    % simualation start
    system('echo $(date): Simset Simulation Script...');
    msg = ['Simualation script started: ', string(datetime)];
    msg = [msg{1}, msg{2}];
    fprintf(fid_log,'%s\n', msg);
   
    % generate and save material index map

    % start parallel pool
    parpool(maxNumCompThreads-2); % leave 2 threads out to minimize SSH lag

    fprintf(fid_log,'Generating material map...\n'); disp('Generating material map...');

    cd(sprintf('%s', umap_processing_dir));
    [~,basename_out] = fileparts(dir_ctac)
    fname_matmap = strcat(dir_umap,'/',basename_out,'.matmap')
    data_matmap = Ct2Matmap(dir_ctac,hu_limits,fname_water,img_size,vox_size,tx_fov);
    fwrite(fopen(fname_matmap,'w'),data_matmap,'int32');
    cd(sprintf('%s', scripts_dir));

    % shut down parallel pool
    delete(gcp('nocreate'))    
    
    
    % generate and save each activity map (w/ corresponding activity table) and phg_params
    fprintf(fid_log, 'Generating activity map(s) and phg_params...\n'); disp('Generating activity map(s) and phg_params...');

    rng('shuffle'); % initialize RNG
    num_to_simulate_per_thread = num_to_simulate/num_threads_per_frame;


    mkdir(sprintf('%s/f%d',dir_amaps_out,iter)); % equivalent to mkdir -p in BASH
    mkdir(sprintf('%s/f%d',dir_phg_params,iter)); % equivalent to mkdir -p in BASH
    mkdir(sprintf('%s/f%d',dir_hist,iter)); % equivalent to mkdir -p in BASH
    mkdir(sprintf('%s/f%d',dir_phg_log,iter)); % equivalent to mkdir -p in BASH

    cd(sprintf('%s', amap_processing_dir));
    fname_amaps_in = dir_amaps_in; 
    fprintf(fid_log, 'fname_amaps_in: %s\n', fname_amaps_in);
    
    fname_amaps_out = sprintf('%s/f%d/amap.%d',dir_amaps_out,iter, iter); 
    fprintf(fid_log, 'fname_amaps_out: %s\n', fname_amaps_out);
    
    fname_phg_act_tables = sprintf('%s/f%d/phg_act_table.%d',dir_amaps_out,iter, iter); 
    fprintf(fid_log, 'fname_phg_act_tables: %s\n', fname_phg_act_tables);
    
    SaveImg2Amap(fname_amaps_in,fname_amaps_out,fname_phg_act_tables,n_bins,img_size,vox_size);
    
    disp(sprintf(fid_log, 'amap_in: %s\n', fname_amaps_in));
    disp(sprintf(fid_log, 'amap_out: %s\n', fname_amaps_out));
    disp(sprintf(fid_log, 'amap_in: %s\n', fname_amaps_in));
       
    cd(sprintf('%s', scripts_dir));
    
    fprintf(fid_log, 'creating phg_params file for each chopped frame...\n'); disp('creating phg_params file for each chopped frame.');
    % create a phg_params file for each chopped frame
    for k = 1:num_threads_per_frame

        % generate and save each *.detparms

        fname_detector_params = sprintf('%s/f%d/f%d_part%02d.detparms.%d',dir_phg_params,iter,iter,k-1,iter); 
        fid_detector_params = fopen(fname_detector_params,'w');
        fprintf(fid_log, 'fname_detector_params: %s\n', fname_detector_params);

        fname_hist = sprintf('%s/f%d/f%d_part%02d.hist.%d',dir_hist,iter,iter,k-1,iter);
        fprintf(fid_log, 'fname_hist: %s\n', fname_hist);
        
        WriteUxDetParms(fid_detector_params,fname_hist);
        fclose(fid_detector_params);

        % generate and save each *.phg_params

        fname_phg_params = sprintf('%s/f%d/f%d_part%02d.phg_params.%d',dir_phg_params,iter,iter,k-1,iter);
        fprintf(fid_log, 'fname_phg_params: %s\n', fname_phg_params);
        
        fid_phg_params = fopen(fname_phg_params,'w');
        x_rng = randi(1e6); % random seed
        WritePhgOptions(fid_phg_params,num_to_simulate_per_thread,x_rng);
        WriteObjectGeometryValues(fid_phg_params,img_size,vox_size);
        WriteUxTargetCylinderInformation(fid_phg_params);
        WriteFilePaths(fid_phg_params,fname_matmap,fname_amaps_out,fname_phg_act_tables,fname_detector_params);
        fclose(fid_phg_params);

    end

    % setting up simulation scripts

    % create Bash scripts for each frame
    fprintf(fid_log, 'Setting up simulation scripts... \n'); disp('Setting up simulation scripts...');

    fname_phg_sh = sprintf('%s/f%d/run_f%d.sh',dir_phg_params,iter,iter);
    fid_phg_sh = fopen(fname_phg_sh,'w');
    fprintf(fid_phg_sh,'%s\n','#!/bin/bash');
    fprintf(fid_phg_sh,'\n');
    fprintf(fid_phg_sh,'%s%d\n','iter_num=',iter);
    fprintf(fid_phg_sh,'\n');
    fprintf(fid_phg_sh,'%s%d\n','num_threads=',num_threads_per_frame);
    fprintf(fid_phg_sh,'\n');
    fprintf(fid_phg_sh,'%s%s\n','fname_simset=',bin_phg);
    fprintf(fid_phg_sh,'\n');
    fprintf(fid_phg_sh,'%s%s\n','dir_phg_params=',dir_phg_params);
    fprintf(fid_phg_sh,'%s%s\n','dir_simset_log=',dir_phg_log);
    fprintf(fid_phg_sh,'\n');
    fprintf(fid_phg_sh,'%s\n','echo "$(date): Starting simulations for iteration #${iter_num}..."');
    fprintf(fid_phg_sh,'\n');
    fprintf(fid_phg_sh,'%s\n','parallel --will-cite -j ${num_threads} "${fname_simset} ${dir_phg_params}/f$(printf %d ${iter_num})/{} > ${dir_simset_log}/f$(printf %d ${iter_num})/{}.log" ::: *.phg_params.$iter_num');
    fprintf(fid_phg_sh,'\n');
    fprintf(fid_phg_sh,'%s\n','echo "$(date): Finished simulations for iteration #${iter_num}."');
    fclose(fid_phg_sh);

    % double-check in log file
    fprintf(fid_log,'%s\n','#!/bin/bash');
    fprintf(fid_log,'\n');
    fprintf(fid_log,'%s%d\n','iter_num=',iter);
    fprintf(fid_log,'\n');
    fprintf(fid_log,'%s%d\n','num_threads=',num_threads_per_frame);
    fprintf(fid_log,'\n');
    fprintf(fid_log,'%s%s\n','fname_simset=',bin_phg);
    fprintf(fid_log,'\n');
    fprintf(fid_log,'%s%s\n','dir_phg_params=',dir_phg_params);
    fprintf(fid_log,'%s%s\n','dir_simset_log=',dir_phg_log);
    fprintf(fid_log,'\n');
    fprintf(fid_log,'%s\n','echo "$(date): Starting simulations for iteration #${iter_num}..."');
    fprintf(fid_log,'\n');
    fprintf(fid_log,'%s\n','parallel --will-cite -j ${num_threads} "${fname_simset} ${dir_phg_params}/f$(printf %d ${iter_num})/{} > ${dir_simset_log}/f$(printf %d ${iter_num})/{}.log" ::: *.phg_params.$iter_num');
    fprintf(fid_log,'\n');
    fprintf(fid_log,'%s\n','echo "$(date): Finished simulations for iteration #${iter_num}."');

    % copy everything to sim node 001

    cmd_dir_umap = ['scp -r ', dir_umap, ' ', user, '@', server_name, ':', dir_umap]; 
    system(cmd_dir_umap);
    pause(0.1);
    cmd_dir_amaps_out = ['scp -r ', dir_amaps_out, ' ', user, '@', server_name, ':', dir_amaps_out]; 
    system(cmd_dir_amaps_out);
    pause(0.1);
    cmd_dir_phg_params = ['scp -r ', dir_phg_params, ' ', user, '@', server_name, ':', dir_phg_params]; 
    system(cmd_dir_phg_params);
    pause(0.1);
    cmd_dir_hist = ['scp -r ', dir_hist, ' ', user, '@', server_name, ':', dir_hist]; 
    system(cmd_dir_hist);
    pause(0.1);
    cmd_dir_phg_log = ['scp -r ', dir_phg_log, ' ', user, '@', server_name, ':', dir_phg_log]; 
    system(cmd_dir_phg_log);
    pause(0.1);

    
 fprintf(fid_log, 'Starting simulations on server node sim-001. This may take a while...\n'); disp('Starting simulations on server node sim-001. This may take a while...');
   
    
    
    
    
    
    
    
end % for loop over iterations

        
    
%% sketch for main program
% v for loop
% v run recon by invoking skript run_lmrecon in the outfolder_server_temp
% v check regularly if this is finished with while loop and check if output file exists and do while not exists. wait 5sec then look again.
% v COPY NOT RENAME output file from recon: 1) img_guess_next_iter_fX and 2) handles.fname_recon, '_f', num2str(handles.reconFrames(kk)), '.intermediate.ITERATION_NUM
% v calculate number to simulate from lm file and take 10 times as many, but at least 1bn (1e9)
% v only run simulation if current iteration number is smaller than handles.osem_iter
% - invoke simulation and conversion and read lm and everything
% - clean up simualtion directory on sim node at the end of each iteration!




%% clean up
%delete(handles_name); 
pause(0.1); 
fprintf(fid_log, 'done.');


fclose(fid_log);
quit;
end % function 

function name_nospace = remove_space_name(name_in)
    name_nospace = name_in; 
    for k = 1:length(name_in)
        if strcmp(name_in(k), ' ')
            name_nospace(k) = '_'; 
        end
    end
end


