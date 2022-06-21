% function to invoke the 5d sinogram creation process and add them together
function getSino5dFromLmFiles(lm_folder_in, crys_eff_file, frame_num, outfolder)

%% parameters and input values
lm2blocksino5d_P_and_D='/home/rbayerlein/code/Recon/scatter_correction/simset/data_processing/lm2blocksino_P_and_D/bin/lm2blocksino_P_and_D';
lm_file=[lm_folder_in, '/lm_reorder_f0_prompts.'];
dim='5D';
d_type='uint16'; % data type of the sinogram entries
img_size_5d = [91,60,112,112,27];
sino5d_r = zeros(img_size_5d);
sino5d_p = zeros(img_size_5d);

%output names
fname_out=[outfolder, '/block_sino_f', num2str(frame_num),'_'];
f_out_p=[fname_out,'prompts.sino5d'];
f_out_r=[fname_out, 'randoms.sino5d'];

%% create sinos 5d
not_done = true;
files_exist = zeros(8,1);

for i = 1 : 8
    cmd_run = [lm2blocksino5d_P_and_D, ' ', lm_file, num2str(i), '.lm ', crys_eff_file, ' ', dim, ' ', outfolder, ' &']; % run in bkgr   
    file2becreated_r = [fname_out, 'randoms.', num2str(i), '.sino5d'];
    file2becreated_p = [fname_out, 'prompts.', num2str(i), '.sino5d'];
    if exist(file2becreated_r, 'file') && exist(file2becreated_p, 'file')
        fprintf('sino exists for file %d. -> omitting.\n', i);
        files_exist(i)=1;
        % check if checkfile exists and create one if not
        temp = [lm_file, num2str(i), '.checkfile'];
        if ~exist(temp, 'file')
            fprintf('creating checkfile for existing data\n');
            cmd_tch = ['touch ', temp]
            system(cmd_tch);
            pause(0.5);
        end
    else
        temp = [lm_file, num2str(i), '.checkfile'];
        if exist(temp, 'file') % delete checkfile if exists
            fprintf('removing checkfile for non-existing data\n');
            cmd_del=['rm ', temp]
            system(cmd_del);
            pause(0.5);
        end
        fprintf('execute sino5d creation using command:\n%s\n', cmd_run);
        system(cmd_run);
    end
    pause(0.5);
end
if sum(files_exist)==8
    not_done = false;
end

%% check if done:
files_done = zeros(8,1);
pause_time = 5.0;
ct = 0;

while not_done
    for i = 1 : 8
        temp = [lm_file, num2str(i), '.checkfile'];
        if exist(temp, 'file')
            files_done(i)=1;
        end
    end
    if sum(files_done) == 8
        not_done = false;
        fprintf('all sinograms created');
    else
        if rem(ct*pause_time,30) == 0 % print out number of files done every 30sec
            fprintf('files done: %d (time elapsed: %d s)\n', sum(files_done), ct*pause_time);
        end
        ct = ct + 1;
        pause(pause_time);
    end
end
fprintf('all files created\n');

%% add them together
%randoms
fprintf('adding randoms...');
if exist(f_out_r, 'file')
    fprintf('file exists:\n%s -> omitting\n', f_out_r);
else
    for i = 1 : 8
        tic
        fprintf('now adding %d\n', i);
        ftemp = [fname_out, 'randoms.', num2str(i), '.sino5d'];
        sino_temp = fread(fopen(ftemp, 'rb'), inf, d_type);
        sino_temp = reshape(sino_temp, img_size_5d);
        sino5d_r = sino5d_r + sino_temp;

        p_end=toc;
        fprintf('\t(processing time for last file number: %d sec.)\n', p_end);
    end
    % save output
    fprintf('saving image\n%s\n', f_out_r);
    fid_out=fopen(f_out_r, 'wb');
    fwrite(fid_out, sino5d_r, d_type);
    clear fid_out;
end

%prompts
fprintf('adding prompts...');
if exist(f_out_p, 'file')
    fprintf('file exists:\n%s -> omitting\n', f_out_p);
else
    for i = 1 : 8
        tic
        fprintf('now adding prompts %d\n', i);
        ftemp = [fname_out, 'prompts.', num2str(i), '.sino5d'];
    %     ftemp = [fname, num2str(i), '.sino5d'];
        sino_temp = fread(fopen(ftemp, 'rb'), inf, d_type);
        sino_temp = reshape(sino_temp, img_size_5d);
        sino5d_p = sino5d_p + sino_temp;

        p_end=toc;
        fprintf('\t(processing time for last file number: %d sec.)\n', p_end);
    end
    % save output
    fprintf('saving image\n%s\n', f_out_p);
    fid_out=fopen(f_out_p, 'wb');
    fwrite(fid_out, sino5d_p, d_type);
    clear fid_out;
end
fprintf('sino5d creation from list-mode files is finished.\n');
end