function tof_cal_linesource



% first, get unit vector of line source position from reconstructed images

% read dicom image of one acquisition
udir = uigetdir(); 

udir_lmdata = uigetdir(); 

toff_sw_fname = 'toff.sw'; 
fid = fopen(toff_sw_fname,'r'); 
t_off = fread(fid, inf, 'double');
t_off = t_off';  

thist = linspace(-500, 500, 300); 
toffhist = histc(t_off, thist); 
figure
hold on
plot(thist, toffhist); 
ylim([0 2000]); 
hold off

t_off = zeros(size(t_off)); 
t_off = 39.0625 .* (rand(1, length(t_off)) - 0.5); 
t_off_new = t_off; 


fclose(fid); 


acq_start = 1; 
acq_end = 1; 

max_frames = 100; 

max_z_diff = 2000;

unit_select1 = 3; 
unit_select2 = 3; 

module_select1 = 0; 
module_select2 = 12; 

trans_block_select1 = 2; 
trans_block_select2 = 62; 

ax_block_select1 = 49; 
ax_block_select2 = 49; 

trans_crys_select1 = 17; %14 to 21
trans_crys_select2 =  437; 

ax_crys_select1 = 0; 
ax_crys_select2 = 0; 

choose_zdiff = false; 
choose_unit = true; 
choose_unit2 = true; 
choose_module = false; 
choose_module2 = false; 
choose_block = false; 
choose_block2 = false;
choose_crystal = false;
choose_crystal2 = false;  


lst = dir(udir)

img_temp = ''; 
i = 1;

while ~contains(img_temp, '.dcm') && i < length(lst)
  img_temp = lst(i).name;
  i = i + 1; 
end 

img_temp = fullfile(udir, img_temp)
dcm_info = dicominfo(img_temp)

 
ring_d = 786;
ax_fov = 1944; 
trans_fov = double(dcm_info.PixelSpacing(1)) * double(dcm_info.Rows)
num_crys = 840; % explorer
num_crys_rings_wgap = 679; 
num_crys_rings_nogap = 672; 
pix_size = [double(dcm_info.PixelSpacing(1)), double(dcm_info.PixelSpacing(1)), double(dcm_info.SliceThickness)]


img_temp = ''; 
i = 1;

img_store = zeros(dcm_info.Rows, dcm_info.Rows, dcm_info.NumberOfSlices); 
zpos_store = zeros(1, dcm_info.NumberOfSlices); 
counter = 1; 

for i = 1:length(lst) 
  img_temp = lst(i).name;
  i = i + 1;
  if contains(img_temp, '.dcm')
    %img_temp_f = fullfile(udir, img_temp); 
    %dcm_info_temp = dicominfo(img_temp_f);
    %img_temp = dicomread(dcm_info_temp);
    %img_temp = zeros(600, 600); 
    %img_temp = double(img_temp) .* double(dcm_info_temp.RescaleSlope); 
    %img_store(:,:,counter) =  img_temp; 
    %zpos_store(1, counter) = double(dcm_info_temp.SliceLocation); 
    %counter = counter + 1; 
  end 
end 

%[q, iSort] = sort(zpos_store); 

%zpos_store = zpos_store(iSort); 

%img_store = img_store(:,:,iSort); 

%zpos_store = zpos_store - min(zpos_store); % range 0 to ~1940 mm
%min(zpos_store)
%max(zpos_store)

%ind_end1 = find(zpos_store > ((max(zpos_store)/2) - 350/2), 1, 'first');
%ind_end2 = find(zpos_store < ((max(zpos_store)/2) + 350/2), 1, 'last'); 

%img_end1 = img_store(:,:,ind_end1);
%img_end2 = img_store(:,:,ind_end2); 

%[y1, yi1] = max(img_end1(:));
%[y2, yi2] = max(img_end2(:));

%[r1, c1] = ind2sub(size(img_end1), yi1); 
%[r2, c2] = ind2sub(size(img_end2), yi2);

%p1 = [c1, r1]; 
%p1 = p1 - (double(dcm_info.Rows) / 2); 
%p1 = p1 .* double(dcm_info.PixelSpacing(1));  
%p1(2) = -1 .* p1(2); 
%p2 = [c2, r2]; 
%p2 = p2 - (double(dcm_info.Rows) / 2); 
%p2 = p2 .* double(dcm_info.PixelSpacing(1)); 
%p2(2) = -1 .* p2(2); 

%p1 = [p1, zpos_store(ind_end1)]
%p2 = [p2, zpos_store(ind_end2)]


% oct 2019
p1 = [1, -48, 647]; 
p2 = [3, -47, 1296]; 

% ge68 line nov 2019
p1 = [4.69, -42.2, 796.96]; 
p2 = [4.7, -46.9, 1146.2]; 



% get line source vector
p = p2 - p1; 
lp = sqrt((p2(1) - p1(1))^2 + (p2(2) - p1(2))^2 + (p2(3) - p1(3))^2); 
p_vec = p ./ lp; 
p_vec = p_vec'; 


%figure
%imagesc(img_store(:,:,ind_end1))

%figure
%imagesc(img_store(:,:,ind_end2))

pause(0.1)





% get positions of each crystal in x,y,z space
theta = (0:(num_crys-1)) ./ num_crys; 
theta = theta .* (2*pi); 
theta = theta - theta(18); 
ry = ring_d .* cos(theta); 
rx = -1.* ring_d .* sin(theta); 
zz = linspace(0, 1944, num_crys_rings_wgap); 


% histogram parameters
tx = 39.0625; 
thist_bins_edges = (-25.25*tx):(tx/2):(25.26*tx); 
thist_bins_centers = (-25*tx):(tx/2):(25.01*tx); 

xhist_bins_edges = -20.25:0.5:20.26;
xhist_bins_centers = -20:0.5:20.1; 

tx_histo = zeros(length(xhist_bins_centers), length(thist_bins_centers));





tol = 0.000; 
min_counts = 1e6;
max_counts = 25e7;  
counts = 0; 

tres_block_module = zeros(1, 14); 

%for NN = 42:55

for NN = 1:5

NN

toff_mean = zeros(size(t_off)); 
toff_mean_counts = zeros(size(t_off)); 



tres_store = zeros(1, (acq_end-acq_start+1)); 
frame_counter = 0; 
acq_counter = 1; 

tx_histo = tx_histo .* 0;
counts = 0;

for n = acq_start:2:acq_end
    n
    
    frame_counter = 0;
    
    while (frame_counter < max_frames) && (counts < max_counts)
        frame_counter
        %[3, 4, 7, 8]
        for nNode = 1:8
            %fname_lm_p = [udir_lmdata, '/lm_',num2str(nNode),'_', num2str(n),'/lm_reorder_f',num2str(frame_counter),'_prompts.',num2str(nNode),'.lm'];
            fname_lm_p = [udir_lmdata, '/lm_reorder_f',num2str(frame_counter),'_prompts.',num2str(nNode),'.lm'];
            %fname_lm_p = [udir_lmdata, '/lm_reorder_f',num2str(frame_counter),'_prompts.',num2str(nNode),'.lm'];
            
            fid_p = fopen(fname_lm_p, 'r');
            
            if fid_p > 0
                
                % load processed listmode data
                d_p = fread(fid_p, inf, 'int16');
                d_p = reshape(d_p, 5, length(d_p)/5);
                fclose(fid_p);
                
                % restrict to central 65 cm as specified in NEMA NU 2 2018
                slice_p = (zz(d_p(4,:) + 1) + zz(d_p(2,:)+1)) / 2;
                val_p = ones(size(slice_p));
                val_p(slice_p > ((1944/2) + 550/2)) = 0;
                val_p(slice_p < ((1944/2) - 550/2)) = 0;
                d_p =  d_p(:, val_p >  0.5);
                
                % max ring difference
                if choose_zdiff
                    z_diff = abs(d_p(4,:) - d_p(2,:));
                    val_p = ones(size(z_diff));
                    val_p(z_diff > max_z_diff) = 0;
                    d_p =  d_p(:, val_p  > 0.5);
                end
                
                % unit select
                if choose_unit
                    unit1 = floor(d_p(2,:) ./ 85);
                    val_p = ones(size(unit1));
                    val_p(unit1 ~= unit_select1) = 0;
                    if choose_unit2
                        unit2 = floor(d_p(4,:) ./ 85);
                        val_p(unit2 ~= unit_select2) = 0;
                    end
                    d_p = d_p(:, val_p >  0.5);
                end
                
                % module select
                if choose_module
                    mod1 = floor(d_p(1,:) ./ 35);
                    val_p = ones(size(mod1));
                    val_p(mod1 ~= module_select1) = 0;
                    if choose_module2
                        mod2 = floor(d_p(3,:) ./ 35);
                        val_p(mod2 ~= module_select2) = 0;
                    end
                    d_p = d_p(:, val_p > 0.5);
                end
                
                % block select
                if choose_block
                    block_trans1 = floor(d_p(1,:) ./ 7);
                    ax_nogap_p1 = d_p(2,:) - floor(d_p(2,:) ./ 85);
                    block_ax1 = floor(ax_nogap_p1 ./ 6);
                    val_p = ones(size(block_trans1));
                    val_p(block_trans1 ~= trans_block_select1) = 0;
                    val_p(block_ax1 ~= ax_block_select1) = 0;
                    if choose_block2
                        block_trans2 = floor(d_p(3,:) ./ 7);
                        ax_nogap_p2 = d_p(4,:) - floor(d_p(4,:) ./ 85);
                        block_ax2 = floor(ax_nogap_p2 ./ 6);
                        val_p(block_trans2 ~= trans_block_select2) = 0;
                        val_p(block_ax2 ~= ax_block_select2) = 0;
                    end
                    d_p = d_p(:,val_p > 0.5);
                end
                
                %if xxx
                %d_p_temp = d_p(:,block_trans1 == block_select1(1));
                %block_trans2_temp = floor(d_p_temp(3,:) ./ 7);
                %block_pair_counts = zeros(1, 120);
                %for pp = 1:length(block_trans2_temp)
                %	block_pair_counts(block_trans2_temp(pp)+1) = block_pair_counts(block_trans2_temp(pp)+1) + 1;
                %end
                %unit_disp = [unit1(1:15); unit2(1:15)]
                %figure
                %plot(block_pair_counts)
                %pause
                %[mm, mmind] = max(block_pair_counts(:));
                %block_select2(1) = mmind - 1
                %choose_block = true;
                %end
                
                
                x_p1 = rx(d_p(1,:)+1);
                x_p2 = rx(d_p(3,:)+1);
                y_p1 = ry(d_p(1,:)+1);
                y_p2 = ry(d_p(3,:)+1);
                z_p1 = zz(d_p(2,:)+1);
                z_p2 = zz(d_p(4,:)+1);
                
                ax_nogap_p1 = d_p(2,:) - floor(d_p(2,:) ./ 85);
                ax_nogap_p2 = d_p(4,:) - floor(d_p(4,:) ./ 85);
                crys_abs_p1 = d_p(1,:) + 840.*ax_nogap_p1;
                crys_abs_p2 = d_p(3,:) + 840.*ax_nogap_p2;
                
                v_p = [];
                v_p(1,:) = x_p2 - x_p1;
                v_p(2,:) = y_p2 - y_p1;
                v_p(3,:) = z_p2 - z_p1;
                
                l_v_p = sqrt((v_p(1,:).^2) + (v_p(2,:).^2) + (v_p(3,:).^2));
                
                v_p(1,:) = v_p(1,:) ./ l_v_p;
                v_p(2,:) = v_p(2,:) ./ l_v_p;
                v_p(3,:) = v_p(3,:) ./ l_v_p;
                
                pcrossv = [];
                pcrossv(1,:) = v_p(2,:).*p_vec(3,1) - v_p(3,:).*p_vec(2,1);
                pcrossv(2,:) = v_p(3,:).*p_vec(1,1) - v_p(1,:).*p_vec(3,1);
                pcrossv(3,:) = v_p(1,:).*p_vec(2,1) - v_p(2,:).*p_vec(1,1);
                
                
                l_pcrossv = sqrt(pcrossv(1,:).^2 + pcrossv(2,:).^2 + pcrossv(3,:).^2);
                
                
                pcrossv(1,:) = pcrossv(1,:) ./ l_pcrossv;
                pcrossv(2,:) = pcrossv(2,:) ./ l_pcrossv;
                pcrossv(3,:) = pcrossv(3,:) ./ l_pcrossv;
                
                
                pdotv = (v_p(1,:) .* p_vec(1,1)) + (v_p(2,:) .* p_vec(2,1)) + (v_p(3,:) .* p_vec(3,1));
                pdotv_squared = pdotv .^ 2;
                f1 = [];
                f1(1,:) = v_p(1,:) - (pdotv*p_vec(1));
                f1(2,:) = v_p(2,:) - (pdotv*p_vec(2));
                f1(3,:) = v_p(3,:) - (pdotv*p_vec(3));
                
                f2 = ((x_p1(1,:) - p1(1)) .* f1(1,:)) + ((y_p1(1,:) - p1(2)) .* f1(2,:)) + ((z_p1(1,:) - p1(3)) .* f1(3,:));
                f2 = f2 ./ (pdotv_squared - 1);
                pca = [];
                pca(1,:) = x_p1(1,:) + (f2.*v_p(1,:));
                pca(2,:) = y_p1(1,:) + (f2.*v_p(2,:));
                pca(3,:) = z_p1(1,:) + (f2.*v_p(3,:));
                
                
                
                t_cor_p = sqrt( ((x_p1(1,:) - pca(1,:)).^2) + ((y_p1(1,:) - pca(2,:)).^2) + ((z_p1(1,:) - pca(3,:)).^2) ) - sqrt( ((x_p2(1,:) - pca(1,:)).^2) + ((y_p2(1,:) - pca(2,:)).^2) + ((z_p2(1,:) - pca(3,:)).^2) );
                t_cor_p = t_cor_p ./ 0.3; %0.3 mm / ps = c
                
                t_p = d_p(5,:).* 39.0625;
                t_p = t_p - t_off(crys_abs_p1+1) + t_off(crys_abs_p2+1) - t_cor_p;
                %t_p(1:100)
                %pause
                
                r_p = ((x_p1(1,:) - p1(1)).*pcrossv(1,:)) + ((y_p1(1,:) - p1(2)).*pcrossv(2,:)) + ((z_p1(1,:) - p1(3)).*pcrossv(3,:));
                
                %r_p(1:100)
                %pause
                
                val_p = ones(size(r_p));
                val_p(r_p < xhist_bins_edges(1)) = 0;
                val_p(r_p > xhist_bins_edges(end)) = 0;
                val_p(t_p < thist_bins_edges(1)) = 0;
                val_p(t_p > thist_bins_edges(end)) = 0;
                
                %d_p = d_p(:,val_p > 0.5);
                r_p = r_p(1, val_p > 0.5);
                t_p = t_p(1, val_p > 0.5);
                
                crys_abs_p1 = crys_abs_p1(1, val_p > 0.5);
                crys_abs_p2 = crys_abs_p2(1, val_p > 0.5);
                
                
                counts = counts + length(t_p);
                
                
                for k = 1:length(t_p)
                    
                    
                    toff_mean(crys_abs_p1(k)+1) = toff_mean(crys_abs_p1(k)+1) + t_p(k);
                    toff_mean(crys_abs_p2(k)+1) = toff_mean(crys_abs_p2(k)+1) - t_p(k);
                    
                    toff_mean_counts(crys_abs_p1(k)+1) = toff_mean_counts(crys_abs_p1(k)+1) + 1; 
                    toff_mean_counts(crys_abs_p2(k)+1) = toff_mean_counts(crys_abs_p2(k)+1) + 1;
                    
                    bin_t = find(t_p(k) > thist_bins_edges, 1, 'last');
                    bin_x = find(r_p(k) > xhist_bins_edges, 1, 'last');
                    tx_histo(bin_x, bin_t) = tx_histo(bin_x, bin_t) + 1;
                    
                end
            end
            
        end
        
        % calc fwhm
        if counts > 4000 && mod(frame_counter, 10) == 0
            tx_histo_sum = zeros(1, size(tx_histo,2));
            for ct = 1:size(tx_histo,2)
                cleft = tx_histo(1,ct);
                cright = tx_histo(end,ct);
                for cx = 1:size(tx_histo, 1)
                    rpos = xhist_bins_centers(cx);
                    tx_histo_sum(ct) = tx_histo_sum(ct) + ( tx_histo(cx, ct) - ((20 - rpos) / 40)*cleft + ((20 + rpos) / 40)*cright );
                end
            end
            
            tx_histo_sum = tx_histo_sum - mean([tx_histo_sum(1), tx_histo_sum(end)]);
            
            
            [mm, mind] = max(tx_histo_sum);
            pk_pts_timing = tx_histo_sum((mind-1):(mind+1));
            p_timing = polyfit((mind-1):(mind+1), pk_pts_timing, 2);
            xfine_t = linspace((mind-1), (mind+1), 1000);
            fit_t = p_timing(1).*(xfine_t.^2) + p_timing(2).*xfine_t + p_timing(3);
            max_fit_timing = max(fit_t(:));
            xleft = find(tx_histo_sum > (max_fit_timing / 2), 1, 'first') - 1;
            xright = find(tx_histo_sum > (max_fit_timing / 2), 1, 'last');
            
            
            fwhm = interp1([tx_histo_sum(xright),tx_histo_sum(xright+1)] ,[xright, xright+1],(max_fit_timing / 2)) - interp1([tx_histo_sum(xleft),tx_histo_sum(xleft+1)] ,[xleft, xleft+1],(max_fit_timing / 2));
            
            fwhm = fwhm * (thist_bins_centers(2) - thist_bins_centers(1))
            if frame_counter < 1
                tres_store(1, acq_counter) = fwhm;
            else
                
                diff_fwhm = abs(fwhm - tres_store(1, acq_counter)) / tres_store(1, acq_counter);
                tres_store(1, acq_counter) = fwhm;
            end
            counts
        end
        
        
        frame_counter = frame_counter + 1;
        
    end
    acq_counter = acq_counter + 1;
    
end



max(toff_mean_counts)
t_off_new = toff_mean ./ toff_mean_counts;

t_off = t_off + t_off_new; 







figure
plot(thist_bins_centers, tx_histo_sum); 
title(num2str(NN))
pause(0.1); 


tres_store = tres_store(1:acq_counter-1); 

fwhm




end




t_off_map = reshape(t_off, 840, 84*8); 

figure
imagesc(t_off_map); 


fb_out = 't_off_new.sw'; 
fid_out = fopen(fb_out, 'w'); 
fwrite(fid_out, t_off, 'double'); 
fclose(fid_out); 

 

         


	

 






