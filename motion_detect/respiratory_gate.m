function respiratory_gate

close all


outfolder = '/mnt/data/rbayerlein/data/motion/20191223/GREEN CRAIG_8239801_131514/PET/RawData/1.2.156.112605.18587648329783.191223211515.9.6660.122659/t2/'

fname_cod = [outfolder, 'cod_xyz.cg']; 

fid = fopen(fname_cod, 'r'); 

dcod = fread(fid, inf, 'double'); 


% ts = str2double(get(handles.t_sampling, 'String')); 
ts = 0.5; 
roi_select = 1;
flim = [0.05 0.5]; % 5 - 16 breaths / minute 

cod_xyz = reshape(dcod, 4, length(dcod)/4); 

cod_xyz(1,:) = cod_xyz(1,:) ./ cod_xyz(4,:); 
cod_xyz(2,:) = cod_xyz(2,:) ./ cod_xyz(4,:); 
cod_xyz(3,:) = cod_xyz(3,:) ./ cod_xyz(4,:); 

cod_xyz = cod_xyz(1:3,:); 

cod_xyz = cod_xyz(:, roi_select:2:end); 

cod_xyz(1,:) = cod_xyz(1,:) - mean(cod_xyz(1,1:50));
cod_xyz(2,:) = cod_xyz(2,:) - mean(cod_xyz(2,1:50));
cod_xyz(3,:) = cod_xyz(3,:) - mean(cod_xyz(3,1:50)); 


tplot = 0:(size(cod_xyz, 2)-1); 
tplot = tplot .* ts; 


figure
hold on
plot(tplot, cod_xyz(1,:), 'r'); 
plot(tplot, cod_xyz(2,:), 'g');
plot(tplot, cod_xyz(3,:), 'b');
hold off


X = cod_xyz(1,:)'; 
Y = cod_xyz(2,:)'; 
Z = cod_xyz(3,:)'; 



YZ = Y + Z; 

figure
hold on
plot(tplot, Y, 'g'); 
plot(tplot, Z, 'b'); 
plot(tplot, YZ, 'k'); 
hold off

pause



Fs = 1/ts;                                           
Fn = Fs/2;                                    
L = length(X)  ;                         
t = (1:L) .* ts; 
                                              
FT_X = fft(X)/L;
FT_Y = fft(Y)/L;
FT_Z = fft(Z)/L;     
                                         
Fv = linspace(0, 1, fix(L/2)+1)*Fn;                  
Iv = 1:length(Fv);                                              
                                                

FLen = 48;                                                          
b_filt = fir1(FLen, flim, 'bandpass');                               

X_filt = fftfilt(b_filt, X);
Y_filt = fftfilt(b_filt, Y); 
Z_filt = fftfilt(b_filt, Z); 



pts_cut = FLen / 2; 

X_filt(1:(end-pts_cut)) = X_filt((pts_cut+1):end); 
Y_filt(1:(end-pts_cut)) = Y_filt((pts_cut+1):end); 
Z_filt(1:(end-pts_cut)) = Z_filt((pts_cut+1):end); 

figure
subplot(2,1,1)
plot(t, X)
title('Original Data')
grid
subplot(2,1,2)
plot(t, X_filt)
title('Filtered Data')
grid



figure
subplot(2,1,1)
plot(t, Y)
title('Original Data')
grid
subplot(2,1,2)
plot(t, Y_filt)
title('Filtered Data')
grid



figure
subplot(2,1,1)
plot(t, Z)
title('Original Data')
grid
subplot(2,1,2)
plot(t, Z_filt)
title('Filtered Data')
grid



%figure
%hold on
%plot(t, X_filt, 'r'); 
%plot(t, Y_filt, 'g'); 
%plot(t, Z_filt, 'b'); 
%hold off




min_mt_free = 3; % 3 seconds
max_diff = 0.3; 

patch_size = round(min_mt_free / ts); 


num_patch = length(X) - patch_size + 1; 

tkeep = zeros(size(t)); 

for k = 1:(length(X)-patch_size+1)

  % Xtemp = X(k:(k+patch_size-1));
   Ytemp = Y(k:(k+patch_size-1));
  % Ztemp = X(k:(k+patch_size-1)); 


  if abs(max(Ytemp) - min(Ytemp)) < max_diff
     tkeep(k:(k+patch_size-1)) = 1; 
  end
end


t_break = [0]; 
mt_state = 0; 
for n = 1:length(tkeep)
   if tkeep(n) ~= mt_state
      t_break = [t_break, t(n)]; 
      mt_state = abs(mt_state - 1); 
      if mt_state == 0
         t_break(end) = t_break(end) - ts; 
      end
   end
end


t_break
t_break = [t_break, t(end)]; 


% save t_break to a text file
tb_text = [];
st = 0; 
time_on = 0; 
for nn = 1:(length(t_break) - 1)
  frame_len = t_break(nn+1) - t_break(nn); 
  txt_temp = ['1,', num2str(frame_len), ',']; 
  tb_text = [tb_text, txt_temp]; 
  if st > 0.5
    time_on = time_on + frame_len; 
  end

  st = abs(st - 1); 
end
time_on / t(end)
tb_text = tb_text(1:(end-1)); 

tb_text



Ykeep = Y(tkeep > .5); 
tkeep_p = t(tkeep > .5); 


figure
subplot(2,1,1)
plot(t, Y)
title('Original Data')
grid
subplot(2,1,2)
hold on
plot(t, Y)
plot(tkeep_p, Ykeep, '.r', 'markersize', 6); 
title('Filtered Data')
grid
hold off





























































 
