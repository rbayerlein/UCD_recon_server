function varargout = respiratory_gui(varargin)
% RESPIRATORY_GUI MATLAB code for respiratory_gui.fig
%      RESPIRATORY_GUI, by itself, creates a new RESPIRATORY_GUI or raises the existing
%      singleton*.
%
%      H = RESPIRATORY_GUI returns the handle to a new RESPIRATORY_GUI or the handle to
%      the existing singleton*.
%
%      RESPIRATORY_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in RESPIRATORY_GUI.M with the given input arguments.
%
%      RESPIRATORY_GUI('Property','Value',...) creates a new RESPIRATORY_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before respiratory_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to respiratory_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help respiratory_gui

% Last Modified by GUIDE v2.5 20-Mar-2020 18:05:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @respiratory_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @respiratory_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before respiratory_gui is made visible.
function respiratory_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to respiratory_gui (see VARARGIN)

% Choose default command line output for respiratory_gui
handles.output = hObject;



% Update handles structure
guidata(hObject, handles);

% UIWAIT makes respiratory_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = respiratory_gui_OutputFcn(hObject, eventdata, handles) 


% Get default command line output from handles structure
varargout{1} = handles.output;


set(handles.outfolder_name, 'String', ''); 
set(handles.dynamic_frames, 'String', ''); 
set(handles.AC_onoff, 'Value', 0);
set(handles.AC_onoff, 'Enable', 'off');
set(handles.norm_onoff, 'Value', 0);
set(handles.norm_onoff, 'Enable', 'off');
set(handles.randoms_onoff, 'Value', 0);
set(handles.randoms_onoff, 'Enable', 'off');
set(handles.scatter_onoff, 'Value', 0);
set(handles.scatter_onoff, 'Enable', 'off');

guidata(hObject, handles); 


% --- Executes on button press in choose_listmode.
function choose_listmode_Callback(hObject, eventdata, handles)

% get listmode file
[fname_choose, pathname_choose] = uigetfile('.raw');

% extract base lm name
ind_ext = strfind(fname_choose, '.raw') - 2; 
fname_choose_base = fname_choose(1:ind_ext);

% get os type
if strcmp(pathname_choose(1), '/') > 0.5
  handles.os_type = 'unix'; 
  disp('unix')
else
  handles.os_type = 'windows'; 
  disp('windows')
end

handles.fname_choose_base = fname_choose_base; 
handles.path_choose = pathname_choose; 
handles.install_dir = cd;
inst_str = 'explorer-master'; 
ind_inst = strfind(handles.install_dir, inst_str);
handles.install_dir = handles.install_dir(1:(ind_inst+length(inst_str)-1)); 
handles.install_dir


% define variables
 
% run time variables
handles.run_recon = '0'; 
handles.write_lm = '0'; 
handles.write_sino = '0';
handles.write_sino_block = '0';  
handles.write_histo_img = '0';

% recon variables 
handles.num_threads = 30;
handles.tker_path = [handles.install_dir, '/reconstruction/tker_57121x57121.pmat']; 
handles.aker_path = [handles.install_dir,'/reconstruction/aker_679x679_rd339.pmat']; 
handles.recon_path = [handles.install_dir, '/reconstruction/lmrecon_explorer/app/']; 
handles.sensitivity_path = ''; 

% image variables
handles.mu_map = []; 
handles.attn_voxel_size = []; 
handles.attn_img_size = []; 
handles.CT_kVp = 0; 
handles.pet_img_size = [239 239 679]; 
handles.pet_voxel_size = [2.85 2.85 2.85]; 





% calc scan length
fname_dicom_header = [handles.fname_choose_base, '1.dcm']; 
fpath_dicom_header = fullfile(handles.path_choose, fname_dicom_header); 
handles.scan_length = calc_scan_length(fpath_dicom_header); 
dyn_frames_init = ['1,',num2str(handles.scan_length)]; 
set(handles.dynamic_frames, 'String', dyn_frames_init);

% get scan date
% TBD

% get norm
% get norm dates
norm_path_base = [handles.install_dir, '/reconstruction/normalization/'];
% find closest but not past
norm_date = '20200101'; 
handles.norm_path_str = [norm_path_base, norm_date]; 
handles.sensitivity_path = [handles.norm_path_str, '/sensimage_239s679-voxelsize2_85mm-sum-brd0_4_v1.smap']; 
if exist(handles.sensitivity_path, 'file')
  s_out = ['Using norm from ', norm_date]; 
  disp(s_out); 
  set(handles.norm_onoff, 'Enable', 'on'); 
  set(handles.norm_onoff, 'Value', 1); 
  set(handles.norm_path, 'String', norm_date); 
end
%guidata(hObject, handles); 

% get AC image
ind_slash = strfind(handles.path_choose, '/'); 
if isempty(ind_slash)
  ind_slash = strfind(handles.path_choose, '\'); 
end

handles.dcm_dir_init = handles.path_choose(1:ind_slash(end-3)); 
handles.dcm_dir_init = [handles.dcm_dir_init, '/Image']; 
handles.AC_path = uigetdir(handles.dcm_dir_init, 'Choose AC image folder'); 
if (handles.AC_path == 0)
  set(handles.AC_onoff, 'Value', 0); 
else
  %pac = genpath('../reconstruction/lmrecon_senimg'); 
  %addpath(pac); 
  %rmpath(pac); 
  %rmpath(pim);
  set(handles.make_sens_img, 'Enable', 'on'); 
  set(handles.ct_path, 'String', handles.AC_path((length(handles.dcm_dir_init)+1):end)); 
end

% check for sensitivity image
fpath_sens = [handles.dcm_dir_init, '/', handles.AC_path(length(handles.dcm_dir_init):end), '.sen_img'] 
if (exist(fpath_sens, 'file')) 
  disp('Found sensitivity image for selected AC image'); 
  handles.sensitivity_path = fpath_sens;
  set(handles.AC_onoff, 'Enable', 'on'); 
  set(handles.AC_onoff, 'Value', 1); 
  set(handles.norm_onoff, 'Value', 1);
  set(handles.norm_path, 'String', handles.AC_path((length(handles.dcm_dir_init)+1):end)); 
else
  % if no sensitivity image found in /image directory, use norm-only sensitivity image as default
  disp('No sensitivity image found'); 
end

guidata(hObject, handles);




function [val, frame_str_nospace] = check_dyn_frames(frame_str)

str_val = ['0123456789,. ']; 

val = 1; 
frame_str_nospace = ''; 
for str_n = 1:length(frame_str)
    str_i = frame_str(str_n); 
    if contains(str_val, str_i) < 0.5
        val = 0; 
    end
    if ~strcmp(str_i, ' ')
        frame_str_nospace = [frame_str_nospace, str_i]; 
    end
end


function nframes = getNumFrames(frame_str)

str_b = 1; 
nframes = 0; 

count_on = 1; 

for str_n = 1:length(frame_str)
    str_i = frame_str(str_n); 
    
    if strcmp(str_i, ',')
        str_e = str_n - 1; 
        numstr = frame_str(str_b:str_e); 
        fnum = str2double(numstr); 
        if count_on > 0
            nframes = nframes + fnum;
        end
        count_on = count_on * -1; 
        str_b = str_n + 1; 
        str_e = str_b; 
    else
        % 
    end
            
end



function scan_dur = calc_scan_length(fname_dcm)

try 
   dcm_info = dicominfo(fname_dcm); 
   scan_dur = double(dcm_info.ActualFrameDuration) / 1000; 
   scan_dur = floor(scan_dur);
catch ME
   scan_dur = 1; 
end





% --- Executes on button press in make_sens_img.
function make_sens_img_Callback(hObject, eventdata, handles)

use_CT = true; 
if (~exist(handles.AC_path, 'dir'))
  ans = questdlg('No CT image chosen, are you sure you want to continue?', '', 'Yes', 'No', 'No'); 
  if ~strcmp(ans, 'Yes')
    return; 
  else
  	use_CT = false; 
  	handles.mu_map = zeros(handles.pet_img_size); 
  end
  
end

ans = questdlg('Generating sensitivity image will take ~3 hours, continue?', '', 'Yes', 'No', 'Yes'); 
if ~strcmp(ans, 'Yes')
  return; 
end

disp('Ok... running sensitivity image backprojection')

% disable all buttons while running
set(handles.choose_listmode, 'Enable', 'off'); 
set(handles.start_pushbutton, 'Enable', 'off'); 
set(handles.make_sens_img, 'Enable', 'off'); 
set(handles.choose_norm, 'Enable', 'off'); 
set(handles.choose_CT, 'Enable', 'off'); 
set(handles.norm_onoff, 'Enable', 'off'); 
set(handles.AC_onoff, 'Enable', 'off'); 
set(handles.randoms_onoff, 'Enable', 'off'); 
set(handles.scatter_onoff, 'Enable', 'off'); 
guidata(hObject, handles); 



% read normalization file
nc_path = [handles.path_choose, handles.fname_choose_base,num2str(1), '.nc']

% check if .nc file exists
if ~exist(nc_path, 'file')
  errordlg('Invalid .nc file path, check raw data'); 
  return;
end

% check nc file against previous sensitivity images
% to do

% extract crystal efficiencies and plane efficiency michelogram
p_sens = genpath('../reconstruction/lmrecon_senimg'); 
addpath(p_sens); 

pim = genpath('../image'); 
addpath(pim); 

crys_eff_840x672 = get_cryseff(nc_path); 
plane_eff_672x672 = get_plaeff(nc_path); 

% convert to axial crystal gap format
[crys_eff_840x679, plane_eff_679x679] = add_axial_gap(crys_eff_840x672, plane_eff_672x672); 

% save efficiencies in the outfolder
handles.fsave_crys_eff = [handles.dcm_dir_init, '/crys_eff_679x840']; 
fid_out1 = fopen(handles.fsave_crys_eff, 'wb'); 
fwrite(fid_out1, crys_eff_840x679, 'float'); 
fclose(fid_out1); 

handles.fsave_plane_eff = [handles.dcm_dir_init, '/plane_eff_679x679']; 
fid_out2 = fopen(handles.fsave_plane_eff, 'wb'); 
fwrite(fid_out2, plane_eff_679x679, 'float'); 
fclose(fid_out2); 


% make mu-map
[handles.mu_map, pars] = make_mumap(handles.AC_path); 
handles.attn_img_size = pars.img_size; 
handles.attn_voxel_size = pars.vox_size; 
handles.CT_kVp = pars.CT_kVp;

% save mu-map
handles.mumap_path = [handles.dcm_dir_init, '/', handles.AC_path(length(handles.dcm_dir_init):end), '_mumap_kVp-',num2str(handles.CT_kVp), '_size-', num2str(handles.attn_img_size(1)),'x',num2str(handles.attn_img_size(2)), 'x', num2str(handles.attn_img_size(3)), '_vox-',num2str(handles.attn_voxel_size(1)), 'x',num2str(handles.attn_voxel_size(2)), 'x', num2str(handles.attn_voxel_size(3)), '.img']; 
fid_musave = fopen(handles.mumap_path, 'wb'); 
fwrite(fid_musave, handles.mu_map, 'float'); 
fclose(fid_musave); 

rmpath(pim); 
rmpath(p_sens); 

% make sensitivity image
handles.sensitivity_path = [handles.dcm_dir_init, handles.AC_path(length(handles.dcm_dir_init):end), '.sen_img']

handles.fname_runsens_sh = [handles.install_dir, '/reconstruction/lmrecon_senimg/run_sens.sh']; 

guidata(hObject, handles); 

disp('Starting 5 parallel runs...'); 

succ = make_runsens_sh(handles); 

if succ < 0
	disp('Error! quit sensitivity image generation'); 
	return;
end

pause(0.1);

system(handles.fname_runsens_sh);

% sum all 5 sensitivity images
sens_img = zeros(handles.pet_img_size); 
sens_img = sens_img(:); 
for k = 0:4
	fsens_temp = [handles.dcm_dir_init, handles.AC_path(length(handles.dcm_dir_init):end),'.',num2str(k), '.sen_img'];
	fid_senstemp = fopen(fsens_temp, 'rb'); 
	sens_temp = fread(fid_senstemp, inf, 'float'); 
	fclose(fid_senstemp); 
	sens_img = sens_img + sens_temp; 
	delete(fsens_temp); 
end

fid_sensout = fopen(handles.sensitivity_path, 'w'); 
fwrite(fid_sensout, sens_img, 'float'); 
fclose(fid_sensout); 


% disable all buttons while running
set(handles.choose_listmode, 'Enable', 'on'); 
set(handles.start_pushbutton, 'Enable', 'on'); 
set(handles.make_sens_img, 'Enable', 'on'); 
set(handles.choose_norm, 'Enable', 'on'); 
set(handles.choose_CT, 'Enable', 'on'); 
set(handles.norm_onoff, 'Enable', 'on'); 
set(handles.AC_onoff, 'Enable', 'on'); 
set(handles.randoms_onoff, 'Enable', 'on'); 
set(handles.scatter_onoff, 'Enable', 'on'); 
set(handles.AC_onoff, 'Value', 1); 


guidata(hObject, handles); 
 

disp('Done sensitivity map!')
 


% --- Executes on button press in start_pushbutton.
function start_pushbutton_Callback(hObject, eventdata, handles)

guidata(hObject, handles);

dyn_frames = get(handles.dynamic_frames, 'String'); 
[val, dyn_frames_nospace] = check_dyn_frames(dyn_frames); 

if val < 0.5
    errordlg('Invalid characters in dynamic frames input, please check'); 
    return
end 

num_frames = getNumFrames(dyn_frames_nospace);
handles.reconFrames = 0:(num_frames-1); 
guidata(hObject, handles); 

if get(handles.AC_onoff, 'Value') > 0.5
  
end
 

%handles.ACimg_raw = ACimg_raw; 

% make mips
%mip_sag = squeeze(mean(ACimg_raw, 1)); 
%mip_cor = squeeze(mean(ACimg_raw, 2));

%figure
%imagesc(mip_sag); 
%axis image; 
%colormap('gray'); 

%figure
%imagesc(mip_cor); 
%axis image; 
%colormap('gray'); 


if strcmp(handles.os_type, 'windows') > 0.5

  fname_bat = 'listmode_process_resp.bat';
  fid_bat = fopen(fname_bat, 'w'); 


  bat_str = '@echo off'; 
  fprintf(fid_bat, '%s\n', bat_str); 
  bat_str = ''; 


  outfolder = [handles.path_choose, get(handles.outfolder_name, 'String'), '\']


  [succ, x, y] = mkdir(outfolder);
  if succ < 0.5
    errordlg('Could not create output folder'); 
    return
  end

  for fnum = 1:8
   
    bat_str = ''; 
   
    infile_fname = [handles.fname_choose_base, num2str(fnum), '.raw']; 
    infile_fullpath = fullfile(handles.path_choose, infile_fname); 
    %infile_fullpath = ['"',infile_fullpath,'"']; 
   
    fname_config = ['..\read_lm\Reconstruction_Parameters_',num2str(fnum),'.txt'];
    fid_config = fopen(fname_config, 'w'); 
   
    fprintf(fid_config, '%s\n%s\n%s\n%s\n%s\n%s\n%s', outfolder, infile_fullpath, handles.run_recon, dyn_frames_nospace, handles.write_lm, handles.write_sino, handles.write_sino_block); 
   
    fclose(fid_config); 
   
    bat_str = ['start ..\read_lm\process_lm_sino_explorer_singlefile.exe ', fname_config]; 
   
   
    fprintf(fid_bat, '%s\n', bat_str);      
    
  end

  fclose(fid_bat); 

  cmd = fname_bat; 
  system(fname_bat); 

end




if strcmp(handles.os_type, 'unix') > 0.5

  outfolder = [handles.path_choose, get(handles.outfolder_name, 'String'), '/'];
  [succ, x, y] = mkdir(outfolder);
  if succ < 0.5
    errordlg('Could not create output folder'); 
    return
  end 
  handles.outfolder = outfolder; 
  guidata(hObject, handles); 

  sh_str = []; 
  for fnum = 1:8
   
    sh_str_temp = ''; 
   
    infile_fname = [handles.fname_choose_base, num2str(fnum), '.raw']; 
    infile_fullpath = fullfile(handles.path_choose, infile_fname); 
    %infile_fullpath = ['"',infile_fullpath,'"']; 
   
    
    fname_config = ['../read_lm/Reconstruction_Parameters_',num2str(fnum)];
    fid_config = fopen(fname_config, 'w'); 
   
    fprintf(fid_config, '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', outfolder, infile_fullpath, handles.run_recon, dyn_frames_nospace, handles.write_lm, handles.write_sino, handles.write_sino_block, handles.write_histo_img); 
   
    fclose(fid_config); 
   
    if fnum == -1 || fnum == -1
      sh_str_temp = ['../read_lm/process_lm_sino_explorer_singlefile_release2  ',fname_config, ' ; ']; 
    else
      sh_str_temp = ['../read_lm/process_lm_sino_explorer_singlefile_release2  ',fname_config, ' & '];
    end
    
    sh_str = [sh_str, sh_str_temp]; 
 
    
  end
  
  sh_str = [sh_str, ' wait; ']; 

  system(sh_str); 
  sh_str = []; 



  if 1
	run_uex_recon(handles); 
  end
	

end








function run_uex_recon(handles)

toRun_sh = ''; 

for m = handles.reconFrames

	combine_str = [handles.recon_path, 'combine_listmode ',handles.outfolder, 'lm_reorder_f',num2str(m),'_prompts.lm ',...
		handles.outfolder, 'lm_reorder_f',num2str(m),'_prompts.1.lm ',...
		handles.outfolder, 'lm_reorder_f',num2str(m),'_prompts.2.lm ',...
		handles.outfolder, 'lm_reorder_f',num2str(m),'_prompts.3.lm ',...
		handles.outfolder, 'lm_reorder_f',num2str(m),'_prompts.4.lm ',...
		handles.outfolder, 'lm_reorder_f',num2str(m),'_prompts.5.lm ',...
		handles.outfolder, 'lm_reorder_f',num2str(m),'_prompts.6.lm ',...
		handles.outfolder, 'lm_reorder_f',num2str(m),'_prompts.7.lm ',...
		handles.outfolder, 'lm_reorder_f',num2str(m),'_prompts.8.lm ']; 

	system(combine_str); 


	make_explorer_lmacc_config(handles, m)
	
	toRun_sh = [handles.outfolder,'run_lmrecon_explorer_f',num2str(m),'.sh']; 
	fid = fopen(toRun_sh,'w');
	str = ['export OMP_NUM_THREADS=',num2str(handles.num_threads),'\n\n'];
	fprintf(fid,str);
	str = '';
	str = [handles.recon_path, 'lmrecon_tof  ',handles.outfolder,'lmacc_scanner_parameter_f',num2str(m),'.cfg'];
	fprintf(fid,str);
	str = '';
	fclose(fid);

	start_cmd = ['sh ', toRun_sh]; 

	system(start_cmd); 

end





function make_explorer_lmacc_config(handles, frame)

fname_cfg = [handles.outfolder, 'lmacc_scanner_parameter_f', num2str(frame), '.cfg']; 

fid1 = fopen(fname_cfg,'w'); 

str = 'detector_ring_diameter = 786.0\n\n';
fprintf(fid1,str);
str = ''; 


str = 'crystal_size = 2.85, 2.85, 2.85, 18.1\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'crystal_gap_size = 0.0, 0.0, 0.0\n\n'; 
fprintf(fid1,str); 
str = ''; 

 
str = 'crystal_array_size = 35, 679\n\n'; 
fprintf(fid1,str); 
str = ''; 

 
str = 'number_of_detector_modules = 24, 1\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'TOF_information = 460, 39.0625\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'number_of_radial_bins = 549\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'image_size = 239, 239, 679\n\n';
fprintf(fid1,str); 
str = ''; 


str = 'voxel_size = 2.85, 2.85, 2.85\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = ['sensitivity = ',handles.sensitivity_path, '\n\n']; %/run/media/meduser/data/lmrecon_explorer/uih3_sensimage_239s679-voxelsize2_85mm-sum-brd0_4-single.smap\n\n'; 
fprintf(fid1,str); 
str = ''; 


%str = 'iPSF_model = /run/media/meduser/data/lmrecon_explorer/aker_679x679_rd339.pmat\n\n'; 


str = ['iPSF_model = ',handles.tker_path, ', ', handles.aker_path, '\n\n']; %/run/media/meduser/data/lmrecon_explorer/tker_57121x57121.pmat, /run/media/meduser/data/lmrecon_explorer/aker_679x679_rd339.pmat\n\n'; 
fprintf(fid1,str); 
str = ''; 


str = 'iterative_algorithm_type = 0\n\n';
fprintf(fid1,str); 
str = ''; 

str = 'warmup_setting = -1, 1\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_strength = 1e-3\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_model_type = 0\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_potential_function_type = 2   #0: Quadratic, 1: Hyperbola, 2: Fair, 3: Huber, 4: hw\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_neighborhood_properties = 1, 1    #size, isotropic  #0: 3x3x3, 1:isotropic  # 1st(0) or 2nd(1), aniso(0) or iso(1)\n\n'; 
fprintf(fid1,str); 
str = '';

str = 'regularizer_buildin_parameter_list = 1e-9\n\n'; 
fprintf(fid1,str); 
str = '';


str = 'iteration_setting = 13, 1, 3    #  subsets, save step, iteration number\n\n'; 
fprintf(fid1,str); 
str = '';

%     #  KEM   2

%# initial_guess = ./lmrecon_tof_475x475x1355_200MBq_20min_PSF.intermediate.1


str = ['input_raw_data_file = ',handles.outfolder,'lm_reorder_f',num2str(frame),'_prompts\n\n']; 
fprintf(fid1,str); 
str = ''; 



str = 'input_raw_data_format_type = 0\n\n'; 
fprintf(fid1,str); 
str = ''; 

str = ['reconstruction_output_setting = ',handles.outfolder,', ./lmrecon_tof_OSEM_f',num2str(frame),'\n\n']; 
fprintf(fid1,str); 

fclose(fid1); 




function combine_lm(handles, m)

fname1 = [handles.outfolder, 'lm_reorder_f', num2str(m), '_prompts.1.lm'];
fname2 = [handles.outfolder, 'lm_reorder_f', num2str(m), '_prompts.2.lm'];
fname3 = [handles.outfolder, 'lm_reorder_f', num2str(m), '_prompts.3.lm'];
fname4 = [handles.outfolder, 'lm_reorder_f', num2str(m), '_prompts.4.lm'];
fname5 = [handles.outfolder, 'lm_reorder_f', num2str(m), '_prompts.5.lm'];
fname6 = [handles.outfolder, 'lm_reorder_f', num2str(m), '_prompts.6.lm'];
fname7 = [handles.outfolder, 'lm_reorder_f', num2str(m), '_prompts.7.lm'];
fname8 = [handles.outfolder, 'lm_reorder_f', num2str(m), '_prompts.8.lm']; 

fid1 = fopen(fname1, 'r');
fid2 = fopen(fname2, 'r');
fid3 = fopen(fname3, 'r');
fid4 = fopen(fname4, 'r');
fid5 = fopen(fname5, 'r');
fid6 = fopen(fname6, 'r');
fid7 = fopen(fname7, 'r');
fid8 = fopen(fname8, 'r');


d1 = fread(fid1, inf, 'short'); 
d2 = fread(fid2, inf, 'short');
d3 = fread(fid3, inf, 'short');
d4 = fread(fid4, inf, 'short');
d5 = fread(fid5, inf, 'short');
d6 = fread(fid6, inf, 'short');
d7 = fread(fid7, inf, 'short');
d8 = fread(fid8, inf, 'short');


fclose(fid1); fclose(fid2); fclose(fid3); fclose(fid4); fclose(fid5); fclose(fid6); fclose(fid7); fclose(fid8);  


d1 = reshape(d1, 5, length(d1)/5);
d2 = reshape(d2, 5, length(d2)/5);
d3 = reshape(d3, 5, length(d3)/5);
d4 = reshape(d4, 5, length(d4)/5);
d5 = reshape(d5, 5, length(d5)/5);
d6 = reshape(d6, 5, length(d6)/5);
d7 = reshape(d7, 5, length(d7)/5);
d8 = reshape(d8, 5, length(d8)/5); 


num_events = size(d1,2) + size(d2,2) + size(d3,2) + size(d4,2) + size(d5,2) + size(d6,2) + size(d7,2) + size(d8,2); 

 

d_out = cat(2, d1, d2, d3, d4, d5, d6, d7, d8); 


r_ind = randperm(num_events); 


d_out = d_out(:, r_ind); 

fname_out = [handles.outfolder, 'lm_reorder_f', num2str(m), '_prompts']; 
fid_out = fopen(fname_out, 'w'); 
fwrite(fid_out, d_out, 'short'); 

fclose(fid_out); 
clear d_out; 
clear d1;
clear d2;
clear d3;
clear d4;
clear d5;
clear d6;
clear d7;
clear d8; 



function succ = make_runsens_sh(handles)

succ = 1; 

str1 = ['${',handles.sensitivity_path,'}'];
str2 = ['${',handles.mumap_path,'}']; 
str3 = ['${',handles.fsave_crys_eff,'}']; 
str4 = ['${',handles.fsave_plane_eff,'}']; 

str_allin = [str1,',',str2,',',str3,',',str4]; 

run_sh_fname = handles.fname_runsens_sh
run1 = [run_sh_fname(1:(end-3)), '.1.sh']
run2 = [run_sh_fname(1:(end-3)), '.2.sh']
run3 = [run_sh_fname(1:(end-3)), '.3.sh']
run4 = [run_sh_fname(1:(end-3)), '.4.sh']
run5 = [run_sh_fname(1:(end-3)), '.5.sh']


%%%% 1
fid_runsh = fopen(run1,'w'); 
if fid_runsh < 0
	succ = -1; 
	disp('Could not open file 1'); 
	return; 
end
fprintf(fid_runsh, '%s\n\n','#!/bin/bash'); 

fprintf(fid_runsh, '%s','P1=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.sensitivity_path); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P2=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.mumap_path); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P3=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.fsave_crys_eff); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P4=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.fsave_plane_eff); 
fprintf(fid_runsh, '%s\n\n', '\''');

ss = ['matlab -nodesktop -nojvm -r "cd ',handles.install_dir, '/reconstruction/lmrecon_senimg/; make_senimg0(${P1},${P2},${P3},${P4}); quit"'];

fprintf(fid_runsh, '%s\n', ss); 

fclose(fid_runsh); 

pause(0.1); 

xx = ['chmod +x ',run1]; 
system(xx);

pause(0.1);


%%%% 2
fid_runsh = fopen(run2,'w'); 
if fid_runsh < 0
	succ = -1; 
	disp('Could not open file 1'); 
	return; 
end
fprintf(fid_runsh, '%s\n\n','#!/bin/bash'); 

fprintf(fid_runsh, '%s','P1=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.sensitivity_path); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P2=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.mumap_path); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P3=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.fsave_crys_eff); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P4=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.fsave_plane_eff); 
fprintf(fid_runsh, '%s\n\n', '\''');

ss = ['matlab -nodesktop -nojvm -r "cd ',handles.install_dir, '/reconstruction/lmrecon_senimg/; make_senimg1(${P1},${P2},${P3},${P4}); quit"'];

fprintf(fid_runsh, '%s\n', ss); 

%fprintf(fid_runsh, '%s\n', 'matlab -nodesktop -nojvm -r "cd /run/media/meduser/data/explorer-master/reconstruction/lmrecon_senimg/; make_senimg1(${P1},${P2},${P3},${P4}); quit"'); 

fclose(fid_runsh); 
pause(0.1);
xx = ['chmod +x ',run2]; 
system(xx);

pause(0.1);



%%%% 2
fid_runsh = fopen(run3,'w'); 
if fid_runsh < 0
	succ = -1; 
	disp('Could not open file 1'); 
	return; 
end
fprintf(fid_runsh, '%s\n\n','#!/bin/bash'); 

fprintf(fid_runsh, '%s','P1=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.sensitivity_path); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P2=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.mumap_path); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P3=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.fsave_crys_eff); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P4=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.fsave_plane_eff); 
fprintf(fid_runsh, '%s\n\n', '\''');

ss = ['matlab -nodesktop -nojvm -r "cd ',handles.install_dir, '/reconstruction/lmrecon_senimg/; make_senimg2(${P1},${P2},${P3},${P4}); quit"'];

fprintf(fid_runsh, '%s\n', ss); 


%fprintf(fid_runsh, '%s\n', 'matlab -nodesktop -nojvm -r "cd /run/media/meduser/data/explorer-master/reconstruction/lmrecon_senimg/; make_senimg2(${P1},${P2},${P3},${P4}); quit"'); 

fclose(fid_runsh); 
pause(0.1);
xx = ['chmod +x ',run3]; 
system(xx);

pause(0.1);



%%%% 2
fid_runsh = fopen(run4,'w'); 
if fid_runsh < 0
	succ = -1; 
	disp('Could not open file 1'); 
	return; 
end
fprintf(fid_runsh, '%s\n\n','#!/bin/bash'); 

fprintf(fid_runsh, '%s','P1=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.sensitivity_path); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P2=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.mumap_path); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P3=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.fsave_crys_eff); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P4=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.fsave_plane_eff); 
fprintf(fid_runsh, '%s\n\n', '\''');

ss = ['matlab -nodesktop -nojvm -r "cd ',handles.install_dir, '/reconstruction/lmrecon_senimg/; make_senimg3(${P1},${P2},${P3},${P4}); quit"'];

fprintf(fid_runsh, '%s\n', ss); 


%fprintf(fid_runsh, '%s\n', 'matlab -nodesktop -nojvm -r "cd /run/media/meduser/data/explorer-master/reconstruction/lmrecon_senimg/; make_senimg3(${P1},${P2},${P3},${P4}); quit"'); 

fclose(fid_runsh); 
pause(0.1);
xx = ['chmod +x ',run4]; 
system(xx);

pause(0.1);



%%%% 5
fid_runsh = fopen(run5,'w'); 
if fid_runsh < 0
	succ = -1; 
	disp('Could not open file 1'); 
	return; 
end
fprintf(fid_runsh, '%s\n\n','#!/bin/bash'); 

fprintf(fid_runsh, '%s','P1=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.sensitivity_path); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P2=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.mumap_path); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P3=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.fsave_crys_eff); 
fprintf(fid_runsh, '%s\n', '\''');

fprintf(fid_runsh, '%s','P4=\'); 
fprintf(fid_runsh, '%s', ''''); 
fprintf(fid_runsh, '%s', handles.fsave_plane_eff); 
fprintf(fid_runsh, '%s\n\n', '\''');

ss = ['matlab -nodesktop -nojvm -r "cd ',handles.install_dir, '/reconstruction/lmrecon_senimg/; make_senimg4(${P1},${P2},${P3},${P4}); quit"'];

fprintf(fid_runsh, '%s\n', ss); 


%fprintf(fid_runsh, '%s\n', 'matlab -nodesktop -nojvm -r "cd /run/media/meduser/data/explorer-master/reconstruction/lmrecon_senimg/; make_senimg4(${P1},${P2},${P3},${P4}); quit"'); 

fclose(fid_runsh); 
pause(0.1);
xx = ['chmod +x ',run5]; 
system(xx);
pause(0.1);



fid_runsh = fopen(run_sh_fname, 'w'); 
if fid_runsh < 0
	succ = -1; 
	disp('Could not open file 1'); 
	return; 
end
fprintf(fid_runsh, '%s\n\n','#!/bin/bash'); 

str = [run1, ' & ', run2, ' & ', run3, ' & ', run4, ' & ', run5, '& wait;']; 

fprintf(fid_runsh, '%s\n', str); 
fclose(fid_runsh); 
pause(0.1);
xx = ['chmod +x ',run_sh_fname]; 
system(xx); 

pause(0.1);






function dynamic_frames_Callback(hObject, eventdata, handles)
% hObject    handle to dynamic_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of dynamic_frames as text
%        str2double(get(hObject,'String')) returns contents of dynamic_frames as a double


% --- Executes during object creation, after setting all properties.
function dynamic_frames_CreateFcn(hObject, eventdata, handles)
% hObject    handle to dynamic_frames (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end






% --- Executes on button press in AC_onoff.
function AC_onoff_Callback(hObject, eventdata, handles)

if (get(handles.AC_onoff, 'Value') > 0.5)
  set(handles.choose_CT, 'Enable', 'on'); 
else
  set(handles.choose_CT, 'Enable', 'off');
  handles.AC_path = 0; 
end
guidata(hObject, handles); 



% --- Executes on button press in choose_CT.
function choose_CT_Callback(hObject, eventdata, handles)




function ct_path_Callback(hObject, eventdata, handles)




% --- Executes during object creation, after setting all properties.
function ct_path_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in norm_onoff.
function norm_onoff_Callback(hObject, eventdata, handles)



% --- Executes on button press in choose_norm.
function choose_norm_Callback(hObject, eventdata, handles)




function norm_path_Callback(hObject, eventdata, handles)



% --- Executes during object creation, after setting all properties.
function norm_path_CreateFcn(hObject, eventdata, handles)


if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in randoms_onoff.
function randoms_onoff_Callback(hObject, eventdata, handles)



% --- Executes on button press in scatter_onoff.
function scatter_onoff_Callback(hObject, eventdata, handles)



function outfolder_name_Callback(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function outfolder_name_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end









 

