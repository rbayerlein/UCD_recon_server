function sino_trace

fdir = uigetdir(); 

sino_size = [91 60 223]; 

frame_select = 1:5; 

slice_select = [100:110]; 

numFrames = (frame_select(end) - frame_select(1) + 1); 

sino_store = zeros([sino_size, numFrames]); 
slice_sum_store = zeros(sino_size(3), numFrames); 

counter = 1; 
for fn = frame_select

	fname_in_p = [fdir, '/sinogram_block_f', num2str(fn), '_prompts.raw'];
	fname_in_r = [fdir, '/sinogram_block_f', num2str(fn), '_randoms.raw']; 

	fid_in_p = fopen(fname_in_p, 'r'); 
	fid_in_r = fopen(fname_in_r, 'r'); 

	sino_p = fread(fid_in_p, inf, 'int32');
	sino_r = fread(fid_in_r, inf, 'int32'); 


	slice_sum_store(:,counter) = squeeze(sum(sum(sino_p, 1),2)) - squeeze(sum(sum(sino_r, 1),2)); 

	sino_store(:,:,:,counter) = sino_p - sino_r;

	counter = counter + 1;  

end

figure
hold on
for sn = slice_select

	plot(1:numFrames, slice_sum_store(slice_select == sn, :)); 
end
hold off

