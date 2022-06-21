function plot_singles_sinos(fdir)


num_blocks = 120*8*14; 
num_bins_module = 21*12*8*8; 



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



t_s_plot = s1(1,:) ./ 1000; 
s1_plot = sum(s1(2:end,:), 1);
s2_plot = sum(s2(2:end,:), 1);
s3_plot = sum(s3(2:end,:), 1);
s4_plot = sum(s4(2:end,:), 1);
s5_plot = sum(s5(2:end,:), 1);
s6_plot = sum(s6(2:end,:), 1);
s7_plot = sum(s7(2:end,:), 1);
s8_plot = sum(s8(2:end,:), 1); 

figure
hold on
plot(t_s_plot, s1_plot, 'r');
plot(t_s_plot, s1_plot, 'g');
plot(t_s_plot, s1_plot, 'b');
plot(t_s_plot, s1_plot, 'c');
plot(t_s_plot, s1_plot, 'm');
plot(t_s_plot, s1_plot, 'k');
plot(t_s_plot, s1_plot, 'y');
plot(t_s_plot, s1_plot, 'r'); 
hold off



















