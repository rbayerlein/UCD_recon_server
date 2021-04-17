function plot_singles_sinos(fdir)


num_blocks = 120*8*14; 
num_bins_module = 21*12*8*8; 


inds_half_sino = []; 




fname_s1 = [fdir, 'singles.1.raw'];
fname_s2 = [fdir, 'singles.2.raw'];
fname_s3 = [fdir, 'singles.3.raw'];
fname_s4 = [fdir, 'singles.4.raw'];
fname_s5 = [fdir, 'singles.5.raw'];
fname_s6 = [fdir, 'singles.6.raw'];
fname_s7 = [fdir, 'singles.7.raw'];
fname_s8 = [fdir, 'singles.8.raw']; 

fname_ps1 = [fdir, 'prompts_sino.1.raw'];
fname_ps2 = [fdir, 'prompts_sino.2.raw'];
fname_ps3 = [fdir, 'prompts_sino.3.raw'];
fname_ps4 = [fdir, 'prompts_sino.4.raw'];
fname_ps5 = [fdir, 'prompts_sino.5.raw'];
fname_ps6 = [fdir, 'prompts_sino.6.raw'];
fname_ps7 = [fdir, 'prompts_sino.7.raw'];
fname_ps8 = [fdir, 'prompts_sino.8.raw']; 


fname_rs1 = [fdir, 'randoms_sino.1.raw'];
fname_rs2 = [fdir, 'randoms_sino.2.raw'];
fname_rs3 = [fdir, 'randoms_sino.3.raw'];
fname_rs4 = [fdir, 'randoms_sino.4.raw'];
fname_rs5 = [fdir, 'randoms_sino.5.raw'];
fname_rs6 = [fdir, 'randoms_sino.6.raw'];
fname_rs7 = [fdir, 'randoms_sino.7.raw'];
fname_rs8 = [fdir, 'randoms_sino.8.raw']; 



fid_s1 = fopen(fname_s1, 'r'); 
fid_s2 = fopen(fname_s2, 'r'); 
fid_s3 = fopen(fname_s3, 'r'); 
fid_s4 = fopen(fname_s4, 'r'); 
fid_s5 = fopen(fname_s5, 'r'); 
fid_s6 = fopen(fname_s6, 'r'); 
fid_s7 = fopen(fname_s7, 'r'); 
fid_s8 = fopen(fname_s8, 'r'); 

fid_ps1 = fopen(fname_ps1, 'r'); 
fid_ps2 = fopen(fname_ps2, 'r'); 
fid_ps3 = fopen(fname_ps3, 'r'); 
fid_ps4 = fopen(fname_ps4, 'r'); 
fid_ps5 = fopen(fname_ps5, 'r'); 
fid_ps6 = fopen(fname_ps6, 'r'); 
fid_ps7 = fopen(fname_ps7, 'r'); 
fid_ps8 = fopen(fname_ps8, 'r'); 


fid_rs1 = fopen(fname_rs1, 'r'); 
fid_rs2 = fopen(fname_rs2, 'r'); 
fid_rs3 = fopen(fname_rs3, 'r'); 
fid_rs4 = fopen(fname_rs4, 'r'); 
fid_rs5 = fopen(fname_rs5, 'r'); 
fid_rs6 = fopen(fname_rs6, 'r'); 
fid_rs7 = fopen(fname_rs7, 'r'); 
fid_rs8 = fopen(fname_rs8, 'r'); 



s1 = fread(fid_s1, inf, 'double');
s2 = fread(fid_s2, inf, 'double');
s3 = fread(fid_s3, inf, 'double');
s4 = fread(fid_s4, inf, 'double');
s5 = fread(fid_s5, inf, 'double');
s6 = fread(fid_s6, inf, 'double');
s7 = fread(fid_s7, inf, 'double');
s8 = fread(fid_s8, inf, 'double'); 



psino1 = fread(fid_ps1, inf, 'double');
psino2 = fread(fid_ps2, inf, 'double');
psino3 = fread(fid_ps3, inf, 'double');
psino4 = fread(fid_ps4, inf, 'double');
psino5 = fread(fid_ps5, inf, 'double');
psino6 = fread(fid_ps6, inf, 'double');
psino7 = fread(fid_ps7, inf, 'double');
psino8 = fread(fid_ps8, inf, 'double'); 


rsino1 = fread(fid_rs1, inf, 'double');
rsino2 = fread(fid_rs2, inf, 'double');
rsino3 = fread(fid_rs3, inf, 'double');
rsino4 = fread(fid_rs4, inf, 'double');
rsino5 = fread(fid_rs5, inf, 'double');
rsino6 = fread(fid_rs6, inf, 'double');
rsino7 = fread(fid_rs7, inf, 'double');
rsino8 = fread(fid_rs8, inf, 'double'); 



inds_read = ones(21*12*8*8, 1); 
counter = 1; 
for k1 = 0:7
	for k2 = 0:7
		for k3 = 1:(21*12)
			if k2 >= k1
				inds_read(counter) = 1;
			else
				inds_read(counter) = 0; 
			end
			counter = counter + 1; 
		end
	end
end
counter = counter - 1; 





s1 = reshape(s1, (num_blocks+1), length(s1)/(num_blocks+1));
s2 = reshape(s2, (num_blocks+1), length(s2)/(num_blocks+1));
s3 = reshape(s3, (num_blocks+1), length(s3)/(num_blocks+1));
s4 = reshape(s4, (num_blocks+1), length(s4)/(num_blocks+1));
s5 = reshape(s5, (num_blocks+1), length(s5)/(num_blocks+1));
s6 = reshape(s6, (num_blocks+1), length(s6)/(num_blocks+1));
s7 = reshape(s7, (num_blocks+1), length(s7)/(num_blocks+1));
s8 = reshape(s8, (num_blocks+1), length(s8)/(num_blocks+1)); 



psino1 = reshape(psino1, (num_bins_module+1), length(psino1)/(num_bins_module+1));
psino2 = reshape(psino2, (num_bins_module+1), length(psino2)/(num_bins_module+1));
psino3 = reshape(psino3, (num_bins_module+1), length(psino3)/(num_bins_module+1));
psino4 = reshape(psino4, (num_bins_module+1), length(psino4)/(num_bins_module+1));
psino5 = reshape(psino5, (num_bins_module+1), length(psino5)/(num_bins_module+1));
psino6 = reshape(psino6, (num_bins_module+1), length(psino6)/(num_bins_module+1));
psino7 = reshape(psino7, (num_bins_module+1), length(psino7)/(num_bins_module+1));
psino8 = reshape(psino8, (num_bins_module+1), length(psino8)/(num_bins_module+1));   



rsino1 = reshape(rsino1, (num_bins_module+1), length(rsino1)/(num_bins_module+1));
rsino2 = reshape(rsino2, (num_bins_module+1), length(rsino2)/(num_bins_module+1));
rsino3 = reshape(rsino3, (num_bins_module+1), length(rsino3)/(num_bins_module+1));
rsino4 = reshape(rsino4, (num_bins_module+1), length(rsino4)/(num_bins_module+1));
rsino5 = reshape(rsino5, (num_bins_module+1), length(rsino5)/(num_bins_module+1));
rsino6 = reshape(rsino6, (num_bins_module+1), length(rsino6)/(num_bins_module+1));
rsino7 = reshape(rsino7, (num_bins_module+1), length(rsino7)/(num_bins_module+1));
rsino8 = reshape(rsino8, (num_bins_module+1), length(rsino8)/(num_bins_module+1)); 



for k = 1:size(psino1,2)
	psino1(2:end,k) = psino1(2:end,k) .* inds_read;
	psino2(2:end,k) = psino2(2:end,k) .* inds_read;
	psino3(2:end,k) = psino3(2:end,k) .* inds_read;
	psino4(2:end,k) = psino4(2:end,k) .* inds_read;
	psino5(2:end,k) = psino5(2:end,k) .* inds_read;
	psino6(2:end,k) = psino6(2:end,k) .* inds_read;
	psino7(2:end,k) = psino7(2:end,k) .* inds_read;
	psino8(2:end,k) = psino8(2:end,k) .* inds_read;
	
	rsino1(2:end,k) = rsino1(2:end,k) .* inds_read;
	rsino2(2:end,k) = rsino2(2:end,k) .* inds_read;
	rsino3(2:end,k) = rsino3(2:end,k) .* inds_read;
	rsino4(2:end,k) = rsino4(2:end,k) .* inds_read;
	rsino5(2:end,k) = rsino5(2:end,k) .* inds_read;
	rsino6(2:end,k) = rsino6(2:end,k) .* inds_read;
	rsino7(2:end,k) = rsino7(2:end,k) .* inds_read;
	rsino8(2:end,k) = rsino8(2:end,k) .* inds_read;
	
end
	
	 



t_s1_plot = s1(1,:) ./ 1000;
t_s2_plot = s2(1,:) ./ 1000;
t_s3_plot = s3(1,:) ./ 1000;
t_s4_plot = s4(1,:) ./ 1000;
t_s5_plot = s5(1,:) ./ 1000;
t_s6_plot = s6(1,:) ./ 1000;
t_s7_plot = s7(1,:) ./ 1000;
t_s8_plot = s8(1,:) ./ 1000;

 
t_psino1_plot = psino1(1,:) ./ 1000;
t_psino2_plot = psino2(1,:) ./ 1000;
t_psino3_plot = psino3(1,:) ./ 1000;
t_psino4_plot = psino4(1,:) ./ 1000;
t_psino5_plot = psino5(1,:) ./ 1000;
t_psino6_plot = psino6(1,:) ./ 1000;
t_psino7_plot = psino7(1,:) ./ 1000;
t_psino8_plot = psino8(1,:) ./ 1000; 

t_rsino1_plot = rsino1(1,:) ./ 1000; 
t_rsino2_plot = rsino2(1,:) ./ 1000; 
t_rsino3_plot = rsino3(1,:) ./ 1000; 
t_rsino4_plot = rsino4(1,:) ./ 1000; 
t_rsino5_plot = rsino5(1,:) ./ 1000; 
t_rsino6_plot = rsino6(1,:) ./ 1000; 
t_rsino7_plot = rsino7(1,:) ./ 1000; 
t_rsino8_plot = rsino8(1,:) ./ 1000; 

s1_plot = sum(s1(2:end,:), 1);
s2_plot = sum(s2(2:end,:), 1);
s3_plot = sum(s3(2:end,:), 1);
s4_plot = sum(s4(2:end,:), 1);
s5_plot = sum(s5(2:end,:), 1);
s6_plot = sum(s6(2:end,:), 1);
s7_plot = sum(s7(2:end,:), 1);
s8_plot = sum(s8(2:end,:), 1); 


psino1_plot = sum(psino1(2:end, :), 1);
psino2_plot = sum(psino2(2:end, :), 1);
psino3_plot = sum(psino3(2:end, :), 1);
psino4_plot = sum(psino4(2:end, :), 1);
psino5_plot = sum(psino5(2:end, :), 1);
psino6_plot = sum(psino6(2:end, :), 1);
psino7_plot = sum(psino7(2:end, :), 1);
psino8_plot = sum(psino8(2:end, :), 1); 

rsino1_plot = sum(rsino1(2:end, :), 1);
rsino2_plot = sum(rsino2(2:end, :), 1);
rsino3_plot = sum(rsino3(2:end, :), 1);
rsino4_plot = sum(rsino4(2:end, :), 1);
rsino5_plot = sum(rsino5(2:end, :), 1);
rsino6_plot = sum(rsino6(2:end, :), 1);
rsino7_plot = sum(rsino7(2:end, :), 1);
rsino8_plot = sum(rsino8(2:end, :), 1); 


tsino1_plot = psino1_plot - rsino1_plot; 
tsino2_plot = psino2_plot - rsino2_plot; 
tsino3_plot = psino3_plot - rsino3_plot; 
tsino4_plot = psino4_plot - rsino4_plot; 
tsino5_plot = psino5_plot - rsino5_plot; 
tsino6_plot = psino6_plot - rsino6_plot; 
tsino7_plot = psino7_plot - rsino7_plot; 
tsino8_plot = psino8_plot - rsino8_plot; 


smean_plot = (s1_plot + s2_plot + s3_plot + s4_plot + s5_plot + s6_plot + s7_plot + s8_plot) ./ 8; 
sdiff1 = s1_plot - smean_plot;
sdiff2 = s2_plot - smean_plot;
sdiff3 = s3_plot - smean_plot;
sdiff4 = s4_plot - smean_plot;
sdiff5 = s5_plot - smean_plot;
sdiff6 = s6_plot - smean_plot;
sdiff7 = s7_plot - smean_plot;
sdiff8 = s8_plot - smean_plot; 

figure
hold on
subplot(2,4,1);
plot(t_s1_plot, s1_plot, 'k');
subplot(2,4,2);
plot(t_s2_plot, s2_plot, 'g');
subplot(2,4,3);
plot(t_s3_plot, s3_plot, 'b');
subplot(2,4,4);
plot(t_s4_plot, s4_plot, 'c');
subplot(2,4,5);
plot(t_s5_plot, s5_plot, 'm');
subplot(2,4,6);
plot(t_s6_plot, s6_plot, 'k');
subplot(2,4,7);
plot(t_s7_plot, s7_plot, 'y');
subplot(2,4,8);
plot(t_s8_plot, s8_plot, 'r'); 
hold off




figure
hold on
subplot(2,4,1);
plot(t_psino1_plot, tsino1_plot, 'k');
subplot(2,4,2);
plot(t_psino2_plot, tsino2_plot, 'g');
subplot(2,4,3);
plot(t_psino3_plot, tsino3_plot, 'b');
subplot(2,4,4);
plot(t_psino4_plot, tsino4_plot, 'c');
subplot(2,4,5);
plot(t_psino5_plot, tsino5_plot, 'm');
subplot(2,4,6);
plot(t_psino6_plot, tsino6_plot, 'k');
subplot(2,4,7);
plot(t_psino7_plot, tsino7_plot, 'y');
subplot(2,4,8);
plot(t_psino8_plot, tsino8_plot, 'r'); 
hold off



figure
hold on
subplot(2,4,1);
plot(t_s1_plot, sdiff1, 'k');
subplot(2,4,2);
plot(t_s2_plot, sdiff2, 'g');
subplot(2,4,3);
plot(t_s3_plot, sdiff3, 'b');
subplot(2,4,4);
plot(t_s4_plot, sdiff4, 'c');
subplot(2,4,5);
plot(t_s5_plot, sdiff5, 'm');
subplot(2,4,6);
plot(t_s6_plot, sdiff6, 'k');
subplot(2,4,7);
plot(t_s7_plot, sdiff7, 'y');
subplot(2,4,8);
plot(t_s8_plot, sdiff8, 'r'); 
hold off


















