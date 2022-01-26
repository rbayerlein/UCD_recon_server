function [fname_s_scaled] = sinogram_scaling_5d(iteration_num, fname_pd, simulation_basename)
%% user-configurable paramters

data_type_pd = 'double';
data_type_t = 'uint16';
data_type_s = 'uint16';

% fname_pd = '/mnt/ssd/rbayerlein/explorer/recon_data_temp/phantom_6beds_78rings_from188_34rings_overlap_220_4it_scatter_corr_20210916/block_sino_f0_pd.raw';
% basename_ts = '/mnt/data/rbayerlein/simset/phantom_20210916_1e10_after_crash_fix_3rd_to_4th_iteration/history/f00000'
basename_ts = simulation_basename;
basename_scaled = ['f', num2str(iteration_num), '_scatters'];

fname_t = [basename_ts, '_trues.sino5d'];
fname_s = [basename_ts, '_scatters.sino5d'];

fprintf('simulation_basename: %s \nbasename_scaled: %s \nfname_t: %s \nfname_s: %s', simulation_basename, basename_scaled, fname_t, fname_s);

% use one of the other setup for mask

fname_ct_fproj = ''; % mask

img_size = [91,60,112,112];
img_size_sino_5d = [91,60,112,112,27];

num_tof_bins = 27;

MUD = 4;
n_ax_blk_per_unit = 14;

% for plotting
ssrb_i = 81; % 112;
axAi = 35; %70; %30; % 39
axBi = 35; %70; %30; % 39
axAj = 59; %84;
axBj = 7; %29;
tof = 14;

% adaptive ssrb
cv_threshold = 0.1;

%% main code

if img_size(3) ~= img_size(4)
    error('img_size(3) is not equal to img_size(4)!');
end

% reading data

data_pd = fread(fopen(fname_pd),data_type_pd);
data_pd = reshape(data_pd,img_size);
data_t_5d = fread(fopen(fname_t),data_type_t);
data_t_5d = reshape(data_t_5d,img_size_sino_5d);
data_s_5d = fread(fopen(fname_s),data_type_s);
data_s_5d = reshape(data_s_5d,img_size_sino_5d);

data_s = zeros(img_size);
data_t = zeros(img_size);

if isempty(fname_ct_fproj)
    disp('fname_ct_fproj is empty. Using default mask...');
    data_mask = ones(img_size);
else
    data_mask = fread(fopen(fname_ct_fproj),data_type_mask);
    data_mask = reshape(data_mask,img_size);
    % udpated for AF (air = 1)
    data_mask(data_mask > threshold) = 0;
    data_mask(data_mask <= threshold) = 1;
end

% initialize intermediate sinograms (mich, ssrb)
sino_pd = zeros(img_size);
sino_t = zeros(img_size);
sino_s = zeros(img_size);
sino_pd_ssrb = zeros(img_size(1),img_size(2),2*img_size(3)-1);
sino_t_ssrb = zeros(img_size(1),img_size(2),2*img_size(3)-1);
sino_s_ssrb = zeros(img_size(1),img_size(2),2*img_size(3)-1);

% get masked sinograms (mich, ssrb)
disp('getting masked sinograms (takes up to 30 sec)');
for axB = 1:img_size(4)
    for axA = 1:img_size(3)
        
        % apply MUD
        if abs(idivide(int32(axA-1),14) - idivide(int32(axB-1),14)) > 4
            continue;
        end
        
        % add tof bins together to get 4d sinogram for scaling
        for tb = 1 : num_tof_bins
            data_t(:,:, axA, axB) = data_t(:,:, axA, axB) + data_t_5d(:,:,axA,axB,tb);
            data_s(:,:, axA, axB) = data_s(:,:, axA, axB) + data_s_5d(:,:,axA,axB,tb);
        end
        
        % apply binary mask
        masked_pd = data_pd(:,:,axA,axB) .* data_mask(:,:,axA,axB);
        masked_t = data_t(:,:,axA,axB) .* data_mask(:,:,axA,axB);
        masked_s = data_s(:,:,axA,axB) .* data_mask(:,:,axA,axB);
        
        % add to mich
        sino_pd(:,:,axA,axB) = sino_pd(:,:,axA,axB) + masked_pd;
        sino_t(:,:,axA,axB) = sino_t(:,:,axA,axB) + masked_t;
        sino_s(:,:,axA,axB) = sino_s(:,:,axA,axB) + masked_s;
        
        % add to ssrb
        sino_pd_ssrb(:,:,axA+axB-1) = sino_pd_ssrb(:,:,axA+axB-1) + masked_pd;
        sino_t_ssrb(:,:,axA+axB-1) =  sino_t_ssrb(:,:,axA+axB-1) + masked_t;
        sino_s_ssrb(:,:,axA+axB-1) = sino_s_ssrb(:,:,axA+axB-1) + masked_s;    

    end
end

disp('setting those T+S entries to zero where P-D is also zero');
% set those T+S entries to zero where P-D is also zero. 
sino_t(sino_pd==0) = 0;
sino_s(sino_pd==0) = 0;

% set those T+S entries to zero where P-D is also zero. 
disp('setting those T+S entries in ssrb sinos to zero where P-D is also zero');
sino_t_ssrb(sino_pd_ssrb==0)=0;
sino_s_ssrb(sino_pd_ssrb==0)=0;        
        
% find scaling factor (mich)
disp('Finding scaling factor (mich)...');
x_optimal = zeros(img_size(3),img_size(4));
for axB = 1:img_size(4)
    for axA = 1:img_size(3)
        % see https://www.mathworks.com/matlabcentral/answers/430640-scale-one-scattered-dataset-to-fit-another-scattered-dataset for more info
        x0 = 1; % scaling - initial guess
        fun = @(x) GetImsse(sino_pd(:,:,axA,axB),sino_t(:,:,axA,axB)+sino_s(:,:,axA,axB),x); % anonymous function
        x_optimal(axA,axB) = fminsearch(fun,x0); % find optimal scaling
        if x_optimal(axA,axB) < 0
            x_optimal(axA,axB) = 0; % apply non-zero negativity constraint to prevent errors
        end
    end
end

disp('Applying MUD (x_optimal)...');
for axB = 1:img_size(4)
    for axA = 1:img_size(3)
        UiA0 = idivide(axA-1,int32(n_ax_blk_per_unit));
        UiB0 = idivide(axB-1,int32(n_ax_blk_per_unit));
        if abs(UiA0-UiB0) > MUD
            x_optimal(axA,axB) = 0;
        end
    end
end

% find scaling factor (ssrb)
disp('Finding scaling factor (ssrb)...');
vec_optimal_ssrb = zeros(2*img_size(3)-1,1);
for k = 1:numel(vec_optimal_ssrb)
    % see https://www.mathworks.com/matlabcentral/answers/430640-scale-one-scattered-dataset-to-fit-another-scattered-dataset for more info
    x0 = 1; % scaling - initial guess
    fun = @(x) GetImsse(sino_pd_ssrb(:,:,k),sino_t_ssrb(:,:,k)+sino_s_ssrb(:,:,k),x); % anonymous function
    vec_optimal_ssrb(k) = fminsearch(fun,x0); % find optimal scaling
    if vec_optimal_ssrb(k) < 0
        vec_optimal_ssrb(k) = 0; % apply non-zero negativity constraint to prevent errors
    end
end

% convert scaling factor from vector to map (ssrb)
x_optimal_ssrb = zeros(img_size(3),img_size(4));


disp('Re-mapping SSRB scaling factors (x_optimal_ssrb)...');
for axB = 1:img_size(4)
    for axA = 1:img_size(3)
        UiA0 = idivide(axA-1,int32(n_ax_blk_per_unit));
        UiB0 = idivide(axB-1,int32(n_ax_blk_per_unit));
        if abs(UiA0-UiB0) > MUD
            x_optimal_ssrb(axA,axB) = 0; % apply MUD to data
        else
            x_optimal_ssrb(axA,axB) = vec_optimal_ssrb(axA+axB-1);
        end
    end
end

% debug - compare diagonal to ssrb scaling line profile
% figure();
% subplot(1,2,1); imshow(x_optimal,[],'InitialMagnification','fit');
% subplot(1,2,2); plot([1:112]/112,diag(x_optimal),'b',[1:223]/223,vec_optimal_ssrb,'r');
% figure();
% subplot(1,2,1); imshow(x_optimal,[],'InitialMagnification','fit');
% subplot(1,2,2); imshow(x_optimal_ssrb,[],'InitialMagnification','fit');

% debug to show sino (mich)
% figure();
% for axB = 1:img_size(4)
%     for axA = 1:img_size(3)
%         if axA == axB
%             subplot(1,2,1); imshow(sino_pd(:,:,axA,axB),[0,max(sino_pd(:,:,axA,axB),[],'all')],'InitialMagnification','fit'); title("P-D axA = " + axA + ", axB = " + axB); colorbar;
%             subplot(1,2,2); imshow(x_optimal(axA,axB)*(sino_t(:,:,axA,axB)+sino_s(:,:,axA,axB)),[],'InitialMagnification','fit'); title("T+S axA = " + axA + ", axB = " + axB); colorbar;
%             pause;
%         end
%     end
% end

% get scaled sinograms (mich)
sino_s_scaled = zeros(img_size);
sino_s_scaled_5d = zeros(img_size_sino_5d);

for axB = 1:img_size(4)
    for axA = 1:img_size(3)
        sino_s_scaled(:,:,axA,axB) = x_optimal(axA,axB) * sino_s(:,:,axA,axB);
        sino_s_scaled_5d(:,:,axA,axB,:) = x_optimal(axA,axB) * single(data_s_5d(:,:,axA,axB,:));
    end
end
clear data_s_5d;

% get scaled sinograms (ssrb)
sino_s_scaled_ssrb = zeros(img_size);
for axB = 1:img_size(4)
    for axA = 1:img_size(3)
        sino_s_scaled_ssrb(:,:,axA,axB) = x_optimal_ssrb(axA,axB) * sino_s(:,:,axA,axB);
    end
end

%% plot

% figure(); subplot(1,2,1); imshow(x_optimal,[0,max(x_optimal,[],'all')],'InitialMagnification','fit'); colorbar; ylabel('AxBlkA'); xlabel('AxBlkB'); title('4-D block sinogram scaling factors'); set(gca,'FontSize',11,'FontWeight','bold');

% figure(); imshow(x_optimal_ssrb,[0,max(x_optimal,[],'all')],'InitialMagnification','fit'); colorbar; ylabel('AxBlkA'); xlabel('AxBlkB'); title('SSRB block sinogram scaling factors'); set(gca,'FontSize',11,'FontWeight','bold');
% 
% figure();
% plot(sino_pd_ssrb(:,1,ssrb_i),'-b'); hold on;
% plot(sino_pd_ssrb(:,1,ssrb_i)-sino_s_ssrb(:,1,ssrb_i)*vec_optimal_ssrb(ssrb_i),'-r');
% plot(sino_s_ssrb(:,1,ssrb_i)*vec_optimal_ssrb(ssrb_i),'-k');
% ylabel('Counts'); xlabel('Bin number');
% legend('Experimental (P-D)','P-D-S','Simulation (S)');
% title('SSRB');
% set(gca,'FontSize',11,'FontWeight','bold');
% 
% figure();
% plot(sino_pd(:,1,axAi,axBi),'-b'); hold on;
% plot(sino_pd(:,1,axAi,axBi)-sino_s_scaled(:,1,axAi,axBi),'-r');
% plot(sino_s_scaled(:,1,axAi,axBi),'-k');
% ylabel('Counts'); xlabel('Bin number');
% legend('Experimental (P-D)','P-D-S','Simulation (S)');
% title('Direct plane');
% set(gca,'FontSize',11,'FontWeight','bold');
% 
% figure();
% plot(sino_pd(:,1,axAj,axBj),'-b'); hold on;
% plot(sino_pd(:,1,axAj,axBj)-sino_s_scaled(:,1,axAj,axBj),'-r');
% plot(sino_s_scaled(:,1,axAj,axBj),'-k');
% ylabel('Counts'); xlabel('Bin number');
% legend('Experimental (P-D)','P-D-S','Simulation (S)');
% title('Oblique plane');
% set(gca,'FontSize',11,'FontWeight','bold');

subplot(1,2,2);
plot(sino_pd_ssrb(:,1,ssrb_i),'-b'); hold on;
plot((sino_t_ssrb(:,1,ssrb_i)+sino_s_ssrb(:,1,ssrb_i))*vec_optimal_ssrb(ssrb_i),'--b');
plot(sino_pd_ssrb(:,1,ssrb_i)-sino_s_ssrb(:,1,ssrb_i)*vec_optimal_ssrb(ssrb_i),'-r');
plot(sino_t_ssrb(:,1,ssrb_i)*vec_optimal_ssrb(ssrb_i),'--r');
plot(sino_s_ssrb(:,1,ssrb_i)*vec_optimal_ssrb(ssrb_i),'-k');
ylabel('Counts'); xlabel('Bin number');
legend('Experimental (P-D)','Simulation (T+S)','P-D-S','Simulation (T)','Simulation (S)');
title('SSRB');
set(gca,'FontSize',11,'FontWeight','bold');


subplot(1,2,1);
plot(sino_pd(:,1,axAi,axBi),'-b'); hold on;
plot((sino_t(:,1,axAi,axBi)+sino_s(:,1,axAi,axBi))*x_optimal(axAi,axBi),'--b');
plot(sino_pd(:,1,axAi,axBi)-sino_s_scaled(:,1,axAi,axBi),'-r');
plot(sino_t(:,1,axAi,axBi)*x_optimal(axAi,axBi),'--r');
plot(sino_s_scaled(:,1,axAi,axBi),'-k');
ylabel('Counts'); xlabel('Bin number');
legend('Experimental (P-D)','Simulation (T+S)','P-D-S','Simulation (T)','Simulation (S)');
title('Direct plane');
set(gca,'FontSize',11,'FontWeight','bold');

% figure();
% plot(sino_pd(:,1,axAj,axBj),'-b'); hold on;
% plot((sino_t(:,1,axAj,axBj)+sino_s(:,1,axAj,axBj))*x_optimal(axAj,axBj),'--b');
% plot(sino_pd(:,1,axAj,axBj)-sino_s_scaled(:,1,axAj,axBj),'-r');
% plot(sino_t(:,1,axAj,axBj)*x_optimal(axAj,axBj),'--r');
% plot(sino_s_scaled(:,1,axAj,axBj),'-k');
% ylabel('Counts'); xlabel('Bin number');
% legend('Experimental (P-D)','Simulation (T+S)','P-D-S','Simulation (T)','Simulation (S)');
% title('Oblique plane');
% set(gca,'FontSize',11,'FontWeight','bold');

%% adaptive SSRB

k_max = 41; % k-th anti diagonal that is affected by MUD
% offset = zeros(k0_max+1,1); % initialization - for removing 0 values due to MUD - use if no LUT - old code
offset = [29:-1:16,15*ones(1,15),14:-1:2]'; % LUT
cv_mud = zeros(numel(vec_optimal_ssrb),1); % coefficient of variation of the anti-diagonal
for k = size(x_optimal,1)-1:-1:-size(x_optimal,1)+1 % k-th anti-diagonal
    anti_diag = diag(fliplr(x_optimal),k);
    if numel(anti_diag) == 0
        error('No elements detected for anti-diagonal!');
    end
    % offset(i0+1) = find(anti_diag,1); % use if no LUT - old code
    if abs(k) <= k_max
        anti_diag = anti_diag(offset(abs(k)+1):end-offset(abs(k)+1)+1);
    end
    cv_mud(k+size(x_optimal_ssrb,1)) = std(anti_diag)/mean(anti_diag);
end

x_optimal_adaptive = x_optimal; % copy
for axB = 1:img_size(4)
    for axA = 1:img_size(3)
        UiA0 = idivide(axA-1,int32(n_ax_blk_per_unit));
        UiB0 = idivide(axB-1,int32(n_ax_blk_per_unit));
        if abs(UiA0-UiB0) > MUD
            continue;
        end
        if cv_mud(axA+axB-1) > cv_threshold
            x_optimal_adaptive(axA,axB) = vec_optimal_ssrb(axA+axB-1);
        end
    end
end

% figure(); imshow(x_optimal_adaptive,[0,max(x_optimal,[],'all')],'InitialMagnification','fit'); colorbar; ylabel('AxBlkA'); xlabel('AxBlkB'); title('Adaptive block sinogram scaling factors'); set(gca,'FontSize',11,'FontWeight','bold');

%% write out

[FILEPATH,NAME,EXT] = fileparts(fname_s);
fname_s_scaled = strcat(FILEPATH,'/',basename_scaled,'_scaled',EXT);
% fname_s_scaled_ssrb = strcat(FILEPATH,'/',NAME,'_scaled_ssrb',EXT);
% fname_s_scaled_adaptive = strcat(FILEPATH,'/',NAME,'_scaled_adaptive',EXT);
fname_scale_factor = strcat(FILEPATH,'/',NAME,'.scale_fac');
% fname_scale_factor_ssrb = strcat(FILEPATH,'/',NAME,'_ssrb.scale_fac');
% fname_scale_factor_adaptive = strcat(FILEPATH,'/',NAME,'_adaptive.scale_fac');
fprintf('Writing %s...\n',fname_s_scaled);
fwrite(fopen(fname_s_scaled,'w'),sino_s_scaled_5d, 'single');   % use single instead of double to limit file size.
% fwrite(fopen(fname_s_scaled_ssrb,'w'),sino_s_scaled_ssrb,'double');
% fwrite(fopen(fname_s_scaled_adaptive,'w'),sino_s_scaled_adaptive,'double');
fwrite(fopen(fname_scale_factor,'w'),x_optimal,'double');
% fwrite(fopen(fname_scale_factor_ssrb,'w'),x_optimal_ssrb,'double');
% fwrite(fopen(fname_scale_factor_adaptive,'w'),x_optimal_adaptive,'double');
% fclose('all');
