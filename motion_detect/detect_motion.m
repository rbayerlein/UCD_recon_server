function motion_frames = detect_motion(cod_fname, thr, t_sampling, roi_num, num_roi)


%fname = 'C:/Documents/EXPLORER/explorer_motion/cod_xyz.cg'; 

%cod_fname = '/mnt/data/rbayerlein/code/listmode_xyz_cod/cod_xyz.cg'; 
fid = fopen(cod_fname, 'r'); 
dd = fread(fid, inf, 'double'); 

dd = reshape(dd, 4, length(dd)/4); 
dd(1,:) = dd(1,:) ./ dd(4,:);
dd(2,:) = dd(2,:) ./ dd(4,:);
dd(3,:) = dd(3,:) ./ dd(4,:); 


XY = dd(1:2,roi_num:num_roi:end); 

fclose(fid); 

thr = 1.0; 
max_km = 6; 
kmm_maxiter = 10; 
min_frame_length = 10; 

num_frame_pts = round(min_frame_length/t_sampling);

xtemp = XY(1,:); 
ytemp = XY(2,:); 

if size(XY, 1) > 2.5
	ztemp = XY(:,3); 
end

mean_baseline_x = mean(xtemp(1:20)); 
mean_baseline_y = mean(ytemp(1:20)); 

std_baseline_x = std(xtemp(1:20)); 
std_baseline_y = std(ytemp(1:20)); 

thr = 2 * max([std_baseline_x, std_baseline_y])

xtemp = xtemp - mean_baseline_x;
ytemp = ytemp - mean_baseline_y; 

if size(xtemp, 2) > size(xtemp, 1)
	xtemp = xtemp';  
	ytemp = ytemp'; 
	%ztemp = ztemp'; 
end

figure
plot(xtemp)

figure
plot(ytemp); 


figure
plot(xtemp, ytemp, '.b'); 

pause(0.5)


% Method 1: iteratively find number of regions
% smooth XY traces
xtemp1 = smooth(xtemp, 3); 
ytemp1 = smooth(ytemp, 3);
%ztemp = smooth(ztemp, 3); 

max_km = 10; 
V = [xtemp1', ytemp1']; 
for kkmm = 1:max_km
	%idx = kmeans(V, kkmm, 'EmptyAction', 'error', 'MaxIter', kmm_maxiter);

	% calculate quality metric
	%q(kkmm) = xxxc; 
end

% choose best kkmm based on quality metric



% Method 2: Find number of motion frames by looking at position jumps.

% savitzy golay  smooth 

order = 5; 
framel = 13; 



xtemp2 = sgolayfilt(xtemp, order, framel);  
ytemp2 = sgolayfilt(ytemp, order, framel); 

% iterative smoothing
%xtemp2 = smooth(xtemp, 3); 
%ytemp2 = smooth(ytemp, 3); 
%xtemp2 = medfilt1(xtemp2, 3); 
%ytemp2 = medfilt1(ytemp2, 3); 
%xtemp2 = smooth(xtemp2, 3); 
%ytemp2 = smooth(ytemp2, 3); 
%xtemp2 = medfilt1(xtemp2, 3); 
%ytemp2 = medfilt1(ytemp2, 3); 


figure
plot(xtemp2)

figure
plot(ytemp2); 



inds_x = abs(xtemp2((num_frame_pts+2):end) - xtemp2(1:(end-num_frame_pts-1))) > thr;
inds_y = abs(ytemp2((num_frame_pts+2):end) - ytemp2(1:(end-num_frame_pts-1))) > thr;


inds_all = inds_x + inds_y;
inds_all(inds_all > 1.5) = 1;



%inds_all = [0; inds_all];

% find groups of 3-5 nearly consecutive numbers,  indicates consistent jump
inds_new = [];
for kn = 1:(length(inds_all)-num_frame_pts)
	temp_sum = sum(inds_all(kn:(kn+num_frame_pts-1)));
	if temp_sum > (num_frame_pts-1.5)
		inds_new = [inds_new, kn];
	end
end

inds_new = inds_new + num_frame_pts;
inds_new = [1, inds_new]


inds_new2 = 1;

for kk = 2:(length(inds_new))

	if inds_new(kk-1) > (inds_new(kk) - 2)
		%%%%
	else
		inds_new2 = [inds_new2, inds_new(kk)]
	end
end
inds_new2 =  [inds_new2,  length(xtemp2)];  



inds_all = inds_new2

frame_length = t_sampling .* (inds_all(2:end) - inds_all(1:(end-1))); 





num_motion_frames = length(inds_all) + 1

%inds_all = inds_all + 2; 
motion_frames_init = inds_all .* t_sampling; 


%km = num_motion_frames; 
%km = 3; % hardcode for now. 

%whos xtemp


figure
hold on
plot(xtemp, 'k'); 
plot(inds_all(1:end), xtemp2(inds_all(1:end)), '.r', 'markersize', 20); 
hold off


return

pause




V = [xtemp', ytemp']; 

%idx = kmeans(V, km, 'EmptyAction', 'error', 'MaxIter', 10); 



idx

Cidx = unique(idx); 

whos Cidx

num_motion_frames = length(Cidx);




% make frames continuous 

% do median filtering with span of xxx
med_span = 5; 
idx_med_smooth = medfilt1(idx, med_span); 

whos idx_med_smooth

%idx_new = idx; 
%idx_new(ceil(med_span/2):(end - ceil(med_span/2))) = idx_med_smooth; 
%idx = idx_new; 
idx = idx_med_smooth; 

for nn = 2:length(idx)-1
	if idx(nn-1) == idx(nn+1)
		idx(nn) = idx(nn-1);
	end
end






cmap = {'r', 'g', 'b', 'c', 'm', 'k','r', 'g', 'b', 'c', 'm', 'k','r', 'g', 'b', 'c', 'm', 'k'}; 
figure
hold on
plot(1:length(xtemp), xtemp, 'k'); 
for i = 1:length(xtemp)
   cc = cmap{idx(i)}; 
   pt = ['.',cc]; 
   plot(i, xtemp(i), pt); 
end
hold off



frame_start_end(1,1) = 0; 
cur_frame = idx(1); 
prev_frame = idx(1);
frame_order = 1; 

for nn = 1:length(idx)
	cur_frame = idx(nn); 
	cur_time = (nn-1) .* t_sampling;
	if cur_frame ~= prev_frame
		frame_start_end(frame_order, 2) = cur_time; 
		frame_start_end(frame_order+1, 1) = cur_time; 
		frame_order = frame_order + 1; 
	end
	prev_frame = cur_frame; 

end

frame_start_end(end, 2) = cur_time; 

motion_frames = frame_start_end








%inds_x = find( (abs(xtemp2(5:end) - xtemp2(1:(end-4))) > thr) ); 
%inds_y = find( abs(ytemp2(5:end) - ytemp2(1:(end-4))) > thr ); 


%whos inds_x
%whos inds_y
 

% what is min separation of x,y inds --> 3?
%for xk = 1:length(inds_x)

%	inds_check = inds_y - inds_x(xk); 
%	i_over = find(inds_check < 3); 
%
%	if ~isempty(i_over)
%		inds_y = [inds_y(1:(i_over-1)); inds_y((i_over+1):end)]; 
%	end
%
%
%end


%inds_all = [inds_x; inds_y]; 
%inds_all = sort(inds_all, 'ascend'); 
%inds_all = [0; inds_all]



%inds_all_new = []; 
%for fi = 1:length(frame_length)
%	if frame_length(fi) > min_frame_length
%		inds_all_new = [inds_all_new; inds_all(fi + 1)]; 
%	end
%end

%inds_all = inds_all_new
