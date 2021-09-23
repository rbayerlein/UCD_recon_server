function scat_corr_recon(frame_num)

%% set up log file
fname_log = 'scatter_recon.log';
if exist(fname_log, 'file')
    delete(fname_log)
end
fid_log = fopen(fname_log, 'w');
cmd_chmod = ['chmod 775 ', fname_log]; 
system(cmd_chmod); 
pause(0.1);
fprintf(fid_log, 'starting scatter corrected list-mode reconstruction\n'); disp('starting scatter corrected list-mode reconstruction');


%% input parameters
fprintf(fid_log, 'reading in parameters from handler\n'); disp('reading in parameters from handler');

handles_name = '/home/rbayerlein/code/explorer-master/handles_scheduler.mat'; 
%handles_name = '/home/rbayerlein/Code/Recon/server/explorer-master/handles_scheduler.mat'; 
pause(0.1); 
h = load(handles_name);
handles = h.handles;
clear h; 
pause(0.1); 
use_scat_corr=handles.scatter_onoff.Value;
fprintf(fid_log, 'correct scatter: %d (should be 1, otherwise this script should not have been invoked!)\n', double(use_scat_corr));

user = 'rbayerlein';
server_name = 'exp-sim-001';

%simset install dir
install_dir = '/home/rbayerlein/code/Recon/scatter_correction/';    % indepentent of handler file

%base_dir for simulation output
slash = strfind(handles.lm_outfolder, '/');
if slash(end) == length(handles.lm_outfolder)
    slash_pos=slash(end-1);
else
    slash_pos = slash(end);
end
dir_base = ['/mnt/data/rbayerlein/simset', handles.lm_outfolder(slash_pos:length(handles.lm_outfolder))];
fprintf(fid_log, 'dir_base: %s\n', dir_base); disp(dir_base);

fprintf(fid_log, 'frame number: %d\n', frame_num);

% u-map
dir_ctac = handles.AC_path_server;
fprintf(fid_log, 'sens image folder path: %s\n', dir_ctac);
tx_fov = 500; % mm

% a-map
outfolder_server_temp = handles.server_temp_dir{frame_num+1};
dir_amaps_in = [handles.server_recon_data_dir, '/', outfolder_server_temp, '/lmrecon_output_f', num2str(frame_num), '.intermediate.1']; 
fprintf(fid_log, 'activity image folder path: %s\n', dir_amaps_in);

% amap image dimensions (CT will be downscaled to match amap)
img_size = [239,239,679];
vox_size = [2.85,2.85,2.85]; % mm

bin_phg = [install_dir, 'simset/2.9.2/bin/phg'];
bin_hist2lm = [install_dir, 'simset/data_processing/hist2lm/bin/hist2lm'];
bin_lm2blocksino =[install_dir, 'simset/data_processing/lm2blocksino_old/bin/lm2blocksino'];
bin_scatter_add_fac =[install_dir, 'simset/data_processing/scatter_add_fac/bin/scatter_add_fac'];

% LUT
fname_lut = [bin_lm2blocksino(1:strfind(bin_lm2blocksino,'/bin')), 'include/index_blockpairs_transaxial_2x91x60_int16'];

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
    
    % caclulate number of events to simulate
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
    
    
    
    
end % for loop over iterations

        
    
%% sketch for main program
% v for loop
% v run recon by invoking skript run_lmrecon in the outfolder_server_temp
% v check regularly if this is finished with while loop and check if output file exists and do while not exists. wait 5sec then look again.
% v COPY NOT RENAME output file from recon: 1) img_guess_next_iter_fX and 2) handles.fname_recon, '_f', num2str(handles.reconFrames(kk)), '.intermediate.ITERATION_NUM
% v calculate number to simulate from lm file and take 10 times as many, but at least 1bn (1e9)
% - invoke simulation and conversion and read lm and everything
% - only run simulation if current iteration number is smaller than handles.osem_iter




%% clean up
%delete(handles_name); 
pause(0.1); 
fprintf(fid_log, 'done.');


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


