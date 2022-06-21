% macro script that does the same as the function
% global_sinogram_scalind_5d.m, but it allows to run individual sections
% after the first execution. Parameters and output values are saved in the
% workspace

%% user-configurable paramters
iteration_num=3;
fname_pd = '/home/rbayerlein/ssd/HumSubj_NA_SC_360s_2it_global_scaling/block_sino_f0_pd.sino5d';
simulation_basename = '/home/rbayerlein/data/explorer/20210714/NERVO_ATILLIO_7696939_144717_Raw/UCD/Scatter_Data_bigger_data_set_20220404/f3.3';

TOFbyTOF_scaling=false; 
% determines whether to use global scaling; true: no global scaling but tof bin by tof bin. 
% false: global scaling, i.e. use all sinograms and all tof bins combined and get ONE scaling factor

data_type_pd = 'uint16';
data_type_t = 'uint16';
data_type_s = 'uint16';

basename_ts = simulation_basename;
basename_scaled = ['f', num2str(iteration_num), '_scatters'];

fname_t = [basename_ts, '_trues.sino5d'];
fname_s = [basename_ts, '_scatters.sino5d'];

fprintf('simulation_basename: %s \nbasename_scaled: %s \nfname_t: %s \nfname_s: %s\n', simulation_basename, basename_scaled, fname_t, fname_s);

num_tof_bins = 27;

img_size_sino_5d = [91,60,112,112,num_tof_bins];

% for plotting
axAi = 50; %70; %30; % 39
axBi = 50; %70; %30; % 39
tof = 14;

%% main code

% reading data
fprintf('reading data from sinograms (takes 60-80s)\n');
fprintf('p-d...\n');
fid_pd = fopen(fname_pd, 'rb');
data_pd = fread(fid_pd,data_type_pd);
data_pd = reshape(data_pd, img_size_sino_5d);

fprintf('trues...\n');
fid_t = fopen(fname_t, 'rb');
data_t_5d = fread(fid_t,data_type_t);
data_t_5d = reshape(data_t_5d, img_size_sino_5d);

fprintf('scatters...\n');
fid_s = fopen(fname_s, 'rb');
data_s_5d = fread(fid_s ,data_type_s);
data_s_5d = reshape(data_s_5d, img_size_sino_5d);
fprintf('all sinograms read in.\n');

disp('setting those T+S entries to zero where P-D is also zero');
data_s_5d(data_pd==0)=0;
data_t_5d(data_pd==0)=0;

%% finding scaling factor(s)
disp('Finding global scaling factor (takes up to 3 minutes)...');
tic
x0 = 1; % scaling - initial guess

% parpool(maxNumCompThreads-2); % leave 2 threads out to minimize SSH lag

if TOFbyTOF_scaling
    x_optimal = zeros(num_tof_bins, 1);%#ok<UNRCH>
    for tofbin = 1 : num_tof_bins
        fprintf('tof bin %d\n', tofbin);
        data_pd_slice = data_pd(:,:,:,:,tofbin);
        data_s_t_slice = data_s_5d(:,:,:,:,tofbin) + data_t_5d(:,:,:,:,tofbin);
%         func = @(x) GetImsse(data_pd_slice, data_s_t_slice, x);
%         x_val = fminsearch(func,x0);
        x_val = sum(sum(sum(sum(data_pd_slice)))) / sum(sum(sum(sum(data_s_t_slice))));
        if x_val<0    % apply non-zero negativity constraint to prevent errors
            x_val=0;
        end
        x_optimal(tofbin,1) = x_val;
    end
else
%     func = @(x) GetImsse(data_pd, data_s_5d + data_t_5d, x); 
%     x_optimal = fminsearch(func, x0);
    x_optimal=sum(sum(sum(sum(sum(data_pd))))) / sum(sum(sum(sum(sum(data_s_5d+data_t_5d)))));
    if x_optimal <0     % apply non-zero negativity constraint to prevent errors
        x_optimal = 0;
    end
end

toc

% shut down parallel pool
% delete(gcp('nocreate'))


%% scaling scatter sinogram
sino_s_scaled_5d = zeros(img_size_sino_5d);
if TOFbyTOF_scaling
    for tofbin = 1 : num_tof_bins
        sino_s_scaled_5d(:,:,:,:,tofbin) = data_s_5d(:,:,:,:,tofbin) * x_optimal(tofbin,1);
    end
else
    sino_s_scaled_5d = data_s_5d * x_optimal; 
    fprintf('global scaling factor: %d\n', x_optimal);
end
delete(gcp('nocreate'))

%% preparation of plotting: Summing up tof bins
axAi = 17; %70; %30; % 39
axBi = 17;
angle = 12;
graph_pd = zeros( img_size_sino_5d(1),1);
graph_t = zeros( img_size_sino_5d(1),1);
graph_s = zeros( img_size_sino_5d(1),1);

for tt =12:12
    for z=0:0
        graph_pd = graph_pd+data_pd(:,angle,axAi+z,axBi+z, tt);
        if TOFbyTOF_scaling
            graph_t = graph_t+data_t_5d(:,angle,axAi+z,axBi+z, tt) *x_optimal(tt,1);
        else
            graph_t = graph_t+data_t_5d(:,angle,axAi+z,axBi+z, tt) *x_optimal;
        end
    %     graph_s = graph_s+data_s_5d(:,10,axAi,axBi, tt) *x_optimal(tt,1);
        graph_s = graph_s+sino_s_scaled_5d(:,angle,axAi+z,axBi+z, tt);

    end
end
%graph_s(graph_s>1.1*graph_pd) = graph_pd(graph_s>1.1*graph_pd);


%% plotting
figure;
plot(graph_pd, '-b'); hold on;
plot(graph_t+graph_s, '--b');
plot(graph_pd-graph_s, '-r');
% plot(graph_t, '--r');
plot(graph_s, '-k');
ylabel('Counts'); xlabel('Bin number');
% legend('Experimental (P-D)','Simulation (T+S)','P-D-S','Simulation (T)','Simulation (S)');
% legend('uncorrected', 'corrected', 'scatters');
legend('uncorrected (P-D)', 'scaled T+S', 'corrected (P-D-S)', 'scatter');
set(gca,'FontSize',11,'FontWeight','bold');

figure;
fprintf('plotting results...');
plot(data_pd(:,1,axAi,axBi, tof),'-b'); hold on;
if TOFbyTOF_scaling
    plot((data_t_5d(:,1,axAi,axBi,tof)+data_s_5d(:,1,axAi,axBi, tof))*x_optimal(tof,1),'--b');
else
    plot((data_t_5d(:,1,axAi,axBi,tof)+data_s_5d(:,1,axAi,axBi, tof))*x_optimal,'--b');
end
plot(data_pd(:,1,axAi,axBi, tof)-sino_s_scaled_5d(:,1,axAi,axBi, tof),'-r');
if TOFbyTOF_scaling
    plot(data_t_5d(:,1,axAi,axBi, tof)*x_optimal(tof,1),'--r');
else
    plot(data_t_5d(:,1,axAi,axBi, tof)*x_optimal,'--r');
end

plot(sino_s_scaled_5d(:,1,axAi,axBi, tof),'-k');
ylabel('Counts'); xlabel('Bin number');
legend('Experimental (P-D)','Simulation (T+S)','P-D-S','Simulation (T)','Simulation (S)');
title('Direct plane');
set(gca,'FontSize',11,'FontWeight','bold');

%% write out

[FILEPATH,NAME,EXT] = fileparts(fname_s);
fname_s_scaled = strcat(FILEPATH,'/',basename_scaled,'_scaled',EXT);
fname_scale_factor = strcat(FILEPATH,'/',NAME,'.scale_fac');
fprintf('Writing\t%s...\nand \t\t%s\n',fname_s_scaled, fname_scale_factor);
% fwrite(fopen(fname_s_scaled,'w'),sino_s_scaled_5d, 'single');   % use single instead of double to limit file size.
% fwrite(fopen(fname_scale_factor,'w'),x_optimal,'double');


%% cleaning up
fclose(fid_pd);
fclose(fid_t);
fclose(fid_s);

