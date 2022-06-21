function varargout = respiratory_gui(varargin)
% RESPIRATORY_GUI MATLAB code for detecting motion in patient PET scans
%      RESPIRATORY_GUI, by itself, creates a new RESPIRATORY_GUI or raises the existing
%      singleton*.
%

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

% Last Modified by GUIDE v2.5 05-Feb-2020 15:40:02

% Created by Eric Berg, UC Davis (2020)

% Default variables defined in the reset_Callback function (bottom of script)

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

% reset 
handles = reset_Callback(handles.reset, eventdata, handles);

guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = respiratory_gui_OutputFcn(hObject, eventdata, handles) 

  % Get default command line output from handles structure
  varargout{1} = handles.output;

  % ??
  set(handles.outfolder_name, 'String', ''); 
  set(handles.dynamic_frames, 'String', ''); 
  set(handles.add_roi, 'Enable', 'off'); 
  guidata(hObject, handles); 


% --- Executes on button press in choose_listmode.
function choose_listmode_Callback(hObject, eventdata, handles)

  % reset 
  handles = reset_Callback(handles.reset, eventdata, handles);
  guidata(hObject, handles); 

  % get file name
  [fname_choose, pathname_choose] = uigetfile('.raw', 'Choose one listmode files', '/mnt/data/rbayerlein/data');
  ind_ext = strfind(fname_choose, '.raw') - 2; 
  fname_choose_base = fname_choose(1:ind_ext);
  handles.fname_choose_base = fname_choose_base; 
  handles.path_choose = pathname_choose; 

  % calc scan length
  fname_dicom_header = [handles.fname_choose_base, '1.dcm']; 
  fpath_dicom_header = fullfile(handles.path_choose, fname_dicom_header); 
  handles.scan_length = calc_scan_length(fpath_dicom_header); 
  num_frames = floor(handles.scan_length / str2double(get(handles.t_sampling, 'String'))); 

  dyn_frames_init = [num2str(num_frames), ',', get(handles.t_sampling, 'String')]; 
  set(handles.dynamic_frames, 'String', dyn_frames_init);

  guidata(hObject, handles);

  % get dicom image 
  ppath = addpath('../image'); 
  ind_slash = strfind(pathname_choose, '/'); 
  if isempty(ind_slash)
    ind_slash = strfind(pathname_choose, '\'); 
  end

  handles.dcm_dir_init = pathname_choose(1:ind_slash(end-3)); 
  handles.dcm_dir_init = [handles.dcm_dir_init, '/Image']; 
  path_choose_dcm = uigetdir(handles.dcm_dir_init, 'Choose dicom image folder for motion detection'); 


  rot_ang = 3; 

  [img_raw, dcmimg_info] = dcm2raw_fast(path_choose_dcm, rot_ang); 
  handles.dcm_info = dcmimg_info; 
  handles.img_raw = img_raw; 

  % make mips
  mip_cor = squeeze(max(img_raw, [], 1)); 
  mip_sag = squeeze(max(img_raw, [], 2)); 

  % correct MIPs to scanner global coordinates (patient position), verify left/right is correct
  % TO DO

  % interpolate to 256 x 830 (hard-coded in listmode processing for now)
  if size(mip_sag, 1) ~= 256 || size(mip_sag, 2) ~= 830
    mip_sag = interp_mip(mip_sag); 
    mip_cor = interp_mip(mip_cor); 
  end

  %  make one image with both mips
  mip_both = cat(1, mip_cor, mip_sag); 
  mip_both(mip_both <= 0) = 0.1; 
  axes(handles.mip_image); 
  imagesc(mip_both, [0 0.2*max(mip_sag(:))]); 
  %imagesc(mip_both, [0 median(mip_both(mip_both > 1))]); 
  %imagesc(log(mip_both)); 
  colormap(flipud(gray(1024)));
  set(gca, 'xticklabel', {''}, 'yticklabel', {''}, 'ticklength', [0 0]); 

  %set(handles.add_roi, 'enable', 'on');

  % get os type
  if strcmp(pathname_choose(1), '/') > 0.5
    handles.os_type = 'unix'; 
    disp('unix')
  else
    handles.os_type = 'windows'; 
    disp('windows')
  end

  guidata(hObject, handles);


% --- Executes on button press in add_roi.
function add_roi_Callback(hObject, eventdata, handles)

  axes(handles.mip_image); 
  guidata(hObject, handles); 
  roi_choose = get(handles.roi_type,'Value');

  if roi_choose == 1  % square
      
    if handles.num_roi < 0.5
      handles.r1_pos_all = [0 0 0 0];
      handles.r2_pos_all = [0 0 0 0];
    else
      handles.r1_pos_all(handles.num_roi,:) = getPosition(handles.r1);
      handles.r2_pos_all(handles.num_roi,:) = getPosition(handles.r2);

      setResizable(handles.r1, false);
      setResizable(handles.r2, false); 

      fcn_lock_r1 = makeConstrainToRectFcn('imrect',[handles.r1_pos_all(handles.num_roi,1) handles.r1_pos_all(handles.num_roi,1) + handles.r1_pos_all(handles.num_roi,3)],[handles.r1_pos_all(handles.num_roi,2) handles.r1_pos_all(handles.num_roi,2) + handles.r1_pos_all(handles.num_roi,4)]);
      setPositionConstraintFcn(handles.r1,fcn_lock_r1); 

      fcn_lock_r2 = makeConstrainToRectFcn('imrect',[handles.r2_pos_all(handles.num_roi,1) handles.r2_pos_all(handles.num_roi,1) + handles.r2_pos_all(handles.num_roi,3)],[handles.r2_pos_all(handles.num_roi,2) handles.r2_pos_all(handles.num_roi,2) + handles.r2_pos_all(handles.num_roi,4)]);
      setPositionConstraintFcn(handles.r2,fcn_lock_r2); 
      
      handles.r1_pos_all(handles.num_roi + 1,:) = [0 0 0 0];
      handles.r2_pos_all(handles.num_roi + 1,:) = [0 0 0 0];
    end
    
    guidata(hObject, handles);
    
    handles.r1_pos = [0 0 0 0]; 
    handles.r2_pos = [0 0 0 0]; 

    handles.r1 = imrect(gca, [20 20 400 200]);
    setColor(handles.r1, handles.cc{handles.num_roi + 1}); 
    fcn_r1 = makeConstrainToRectFcn('imrect',[0 830],[0 256]);
    setPositionConstraintFcn(handles.r1,fcn_r1); 

    handles.r2 = imrect(gca, [20 276 400 200]);
    setColor(handles.r2, handles.cc{handles.num_roi + 1});
    fcn_r2 = makeConstrainToRectFcn('imrect',[0 830],[256 512]);
    setPositionConstraintFcn(handles.r2,fcn_r2); 

    guidata(hObject, handles);

    handles.r1_pos = getPosition(handles.r1); 
    handles.r2_pos = getPosition(handles.r2); 
    setPosition(handles.r1, [handles.r2_pos(1) handles.r1_pos(2) handles.r2_pos(3) handles.r1_pos(4)]); 
    handles.r1_pos = getPosition(handles.r1);

    guidata(hObject, handles);

    addNewPositionCallback(handles.r1,@myfun1);
    addNewPositionCallback(handles.r2,@myfun2);

    guidata(hObject, handles);
      
      
  elseif roi_choose == 2  % sphere
    errordlg('Not ready yet'); 
    return; 
    
  elseif roi_choose == 3 %cylinder 
    errordlg('Not ready yet!'); 
    return; 
      
  else
    errordlg('Invalid ROI'); 
    return; 
  end

  handles.num_roi = handles.num_roi + 1; 

  set(handles.plot_cod, 'Enable', 'on'); 
  set(handles.detect_motion, 'Enable', 'on'); 

  if handles.num_roi == 1
      set(handles.roi_check_all, 'Enable', 'on'); 
      set(handles.roi_check_all, 'Value', 0); 
      set(handles.roi_check_1, 'Enable', 'on');
      set(handles.roi_check_1, 'Value', 1); 
  end

  if handles.num_roi == 2
      set(handles.roi_check_2, 'Enable', 'on');
      set(handles.roi_check_2, 'Value', 1); 
  end

  if handles.num_roi == 3
      set(handles.roi_check_3, 'Enable', 'on');
      set(handles.roi_check_3, 'Value', 1); 
  end

  if handles.num_roi == 4
      set(handles.roi_check_4, 'Enable', 'on');
      set(handles.roi_check_4, 'Value', 1); 
  end

  if handles.num_roi == 5
      set(handles.roi_check_5, 'Enable', 'on');
      set(handles.roi_check_5, 'Value', 1); 
  end

  if handles.num_roi == 6
      set(handles.roi_check_6, 'Enable', 'on');
      set(handles.roi_check_6, 'Value', 1); 
  end

  if handles.num_roi == 7
      set(handles.roi_check_7, 'Enable', 'on');
      set(handles.roi_check_7, 'Value', 1); 
  end

  if handles.num_roi == 8
      set(handles.roi_check_8, 'Enable', 'on');
      set(handles.roi_check_8, 'Value', 1); 
  end

  guidata(hObject, handles);



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

  %if handles.num_roi < 0.5
    %errordlg('Must draw at least one ROI'); 
    %return
  %end

  % make output folder
  outfolder = [handles.path_choose, get(handles.outfolder_name, 'String'), '/'];
  [succ, x, y] = mkdir(outfolder);
  if succ < 0.5
    errordlg('Could not create output folder'); 
    return
  end 
  handles.outfolder_path = outfolder;
  guidata(hObject, handles); 

  redo_lm = true; 

  fname_cod_xyz = [handles.outfolder_path, 'cod_xyz.cg']; 
  if exist(fname_cod_xyz, 'file'); 
    redo_lm_str = questdlg('Overwrite existing file?','COD_XYZ already exists','No'); 
    if strcmp(redo_lm_str, 'No')
      redo_lm = false; 
      %par_fname = [outfolder, 'Reconstruction_Parameters_1']; 
      %fid_par = fopen(par_fname, 'r'); 
      %il = 1;
      %tline = fgetl(fid_par); 
      %B{il} = tline; 
      %while ischar(tline)
      %  il = il+1; 
      %  tline = fgetl(fid_par); 
      %  B{il} = tline;
      %end
      %fclose(fid_par);
      %dyn_frame = B{4}; 

      %ic = strfind(dyn_frame, ','); 
      %ts_old = dyn_frame((ic+1):end)
      %set(handles.dynamic_frames, 'String', dyn_frame); 
      %et(handles.t_sampling, 'String', ts_old); 
    end
  end


  if redo_lm

    % make roi mask image using rois
    mask_img = zeros(256, 256, 830) - 1; 

    roi_size = 50; %100 mm ROI cubes
    roi_vox_size = [20 20 20];
    num_roi_t = floor(256/20); 
    num_roi_z = floor(830/20); 
    t_span = num_roi_t * 20; 
    z_span = num_roi_z * 20; 
     
    res_t = 256 - t_span; 
    res_z = 830 - z_span; 

    sp_t = floor(res_t / 2) + 1; 
    sp_z = floor(res_z / 2) + 1; 

    kr = 1; 

    for zr = 1:num_roi_z
       for yr = 1:num_roi_t
          for xr = 1:num_roi_t
             p1 = [(sp_t + 20*(yr-1)), (sp_t + 20*(xr-1)), (sp_z + 20*(zr-1))]; 
    	       mask_img(p1(1):(p1(1)+19), p1(2):(p1(2)+19), p1(3):(p1(3)+19)) = kr; 
    	       kr = kr + 1; 
          end
       end
    end 


    handles.num_patch = kr - 1; 
    handles.mask_img = mask_img; 



    % save mask image temporarily
    handles.fname_mask_img = fullfile(handles.dcm_dir_init, get(handles.outfolder_name, 'String'));
    handles.fname_mask_img = [handles.fname_mask_img, '_roi_mask.img'];
    handles.fname_mask_img_cpp = ['"', handles.fname_mask_img, '"']; 

    fid_mask = fopen(handles.fname_mask_img, 'w'); 
    fwrite(fid_mask, mask_img, 'int'); 
    fclose(fid_mask); 

    handles.fname_raw_img = fullfile(handles.dcm_dir_init, get(handles.outfolder_name, 'String')); 
    handles.fname_raw_img = [handles.fname_raw_img, '_raw_img.img']; 
    handles_fname_raw_img_cpp = ['"', handles.fname_raw_img, '"']; 


    fid_rawimg = fopen(handles.fname_raw_img, 'w'); 
    fwrite(fid_rawimg, handles.img_raw, 'float'); 
    fclose(fid_rawimg); 




    % look for sensitivity image
    found_sens = false; 
    fpath_sensimg = handles.dcm_dir_init; 
    lst_sens = dir(handles.dcm_dir_init); 
    for kk = 1:length(lst_sens)
      if contains(lst_sens(kk).name, '.sen_img'); 
        fpath_sensimg = [handles.dcm_dir_init, '/', lst_sens(kk).name]; 
        found_sens = true; 
        disp('found sensitivity image'); 
      end
    end


    fpath_sensimg



    % make 256x256x830 sens image
    sens_vox_size_orig = [2.85 2.85 2.85]; 
    sens_img_size_orig = [239 239 679]; 

    sens_vox_size_new = [2.344 2.344 2.344]; 
    sens_img_size_new = [256 256 830]; 

    fid_sens = fopen(fpath_sensimg, 'rb'); 
    sens_img = fread(fid_sens, inf, 'float'); 
    sens_img = sens_img ./ median(sens_img(:));
    sens_img = reshape(sens_img, sens_img_size_orig); 

    for zz = 1:size(sens_img,3)
      sens_img(:,:,zz) = rot90(sens_img(:,:,zz), 1); 
    end

    x_b1 = ([0:(sens_img_size_new(1)-1)] - ceil(sens_img_size_new(1)/2)) .* sens_vox_size_new(1);
    y_b1 = ([0:(sens_img_size_new(2)-1)] - ceil(sens_img_size_new(2)/2)) .* sens_vox_size_new(2);
    z_b1 = ([0:(sens_img_size_new(3)-1)] - ceil(sens_img_size_new(3)/2)) .* sens_vox_size_new(3); 



    x_b2 = ([0:(sens_img_size_orig(1)-1)] - ceil(sens_img_size_orig(1)/2)) .* sens_vox_size_orig(1);
    y_b2 = ([0:(sens_img_size_orig(2)-1)] - ceil(sens_img_size_orig(2)/2)) .* sens_vox_size_orig(2);
    z_b2 = ([0:(sens_img_size_orig(3)-1)] - ceil(sens_img_size_orig(3)/2)) .* sens_vox_size_orig(3);


    xlim_1 = find(x_b2 < x_b1(1), 1, 'last'); 
    xlim_2 = find(x_b2 > x_b1(end), 1, 'first'); 

    ylim_1 = find(y_b2 < y_b1(1), 1, 'last'); 
    ylim_2 = find(y_b2 > y_b1(end), 1, 'first');

    %zlim_1 = find(z_b2 < z_b1(1), 1, 'last'); 
    %zlim_2 = find(z_b2 > z_b1(end), 1, 'first');



    sens_img = sens_img(ylim_1:ylim_2, xlim_1:xlim_2, :); 
    x_b2 = x_b2(xlim_1:xlim_2);
    y_b2 = y_b2(ylim_1:ylim_2);
    %z_b2 = z_b2(zlim_1:zlim_2); 

    whos x_b2
    whos y_b2
    whos z_b2

    [X,Y,Z] = meshgrid(y_b2,x_b2,z_b2);
    [Xq,Yq,Zq] = meshgrid(y_b1,x_b1,z_b1);  



    sens_img_new = interp3(X,Y,Z,sens_img,Xq,Yq,Zq); 


    figure
    imagesc(sens_img_new(:,:,300)); 

     pause(0.5)



    whos sens_img_new

    fpath_cod_sensimg = [outfolder, 'cod_sensimg.img'];
    fid_sensout = fopen(fpath_cod_sensimg, 'w'); 
    fwrite(fid_sensout, sens_img_new, 'float'); 
    fclose(fid_sensout); 

    fpath_cod_sensimg_cpp = ['"', fpath_cod_sensimg, '"']; 


    
    guidata(hObject, handles); 



    if strcmp(handles.os_type, 'unix') > 0.5

      sh_str = []; 
      for fnum = 1:8
       
        sh_str_temp = ''; 
      
        infile_fname = [handles.fname_choose_base, num2str(fnum), '.raw']; 
        infile_fullpath = fullfile(handles.path_choose, infile_fname); 
        
        % fname_config = ['../read_lm/Reconstruction_Parameters_',num2str(fnum)];
        fname_config = [outfolder, 'Reconstruction_Parameters_',num2str(fnum)];
        fname_config_cpp = ['"', fname_config, '"']; 
        fid_config = fopen(fname_config, 'w'); 
        fprintf(fid_config, '%s\n%s\n%s\n%s\n%s\n%s\n%s\n%s', outfolder, infile_fullpath, handles.run_recon, dyn_frames_nospace, handles.write_lm, handles.write_sino, handles.write_sino_block, handles.write_histo_img); 
        fclose(fid_config); 

        sh_str_temp = ['../read_lm/process_lm_cod_xyz_explorer_singlefile_fast  ',fname_config_cpp, '  ', handles.fname_mask_img_cpp, '  ', fpath_cod_sensimg_cpp, ' & '];
        sh_str = [sh_str, sh_str_temp];  

      end
      
      sh_str = [sh_str, ' wait; ']; 

      system(sh_str); 
      sh_str = []; 

      
      update_cod_xyz(handles);
      %plot_cod_xyz_gui(handles); 
      %guidata(hObject, handles);

    end

  else
    handles.fname_mask_img = fullfile(handles.dcm_dir_init, get(handles.outfolder_name, 'String'));
    handles.fname_mask_img = [handles.fname_mask_img, '_roi_mask.img'];
    fid_mask = fopen(handles.fname_mask_img, 'r'); 
    handles.mask_img = fread(fid_mask, inf, 'int'); 
    handles.mask_img = reshape(handles.mask_img, 256, 256, 830); 
    handles.num_patch = max(handles.mask_img(:)); 
    fclose(fid_mask); 
  end

   
  handles.roi_state = get_roi_status(handles); 
  guidata(hObject, handles);


  fname_cod_xyz = [handles.outfolder_path, 'cod_xyz.cg']; 
  fid = fopen(fname_cod_xyz, 'r'); 
  cod_xyz = fread(fid, inf, 'double'); 
  fclose(fid); 

  if isempty(cod_xyz)|| mod(length(cod_xyz), 4) ~= 0
    errordlg('Invalid file size, must be multiple of 4'); 
    return
  end 

  cod_xyz = reshape(cod_xyz, 4, length(cod_xyz)/4); 

  cod_xyz(1,:) = cod_xyz(1,:) ./ cod_xyz(4,:); 
  cod_xyz(2,:) = cod_xyz(2,:) ./ cod_xyz(4,:); 
  cod_xyz(3,:) = cod_xyz(3,:) ./ cod_xyz(4,:); 

  handles.cod_xyz = cod_xyz(1:3,:); 
  handles.cod_xyz_counts = cod_xyz(4,:); 

  %set(handles.plot_cod, 'Enable', 'on'); 
  set(handles.add_roi, 'Enable', 'on'); 
  guidata(hObject, handles); 


function outfolder_name_Callback(hObject, eventdata, handles)


function outfolder_name_CreateFcn(hObject, eventdata, handles)

  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
      set(hObject,'BackgroundColor','white');
  end


function dynamic_frames_Callback(hObject, eventdata, handles)


function dynamic_frames_CreateFcn(hObject, eventdata, handles)

  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
      set(hObject,'BackgroundColor','white');
  end


function t_sampling_Callback(hObject, eventdata, handles)

  num_frames = floor(handles.scan_length / str2double(get(handles.t_sampling, 'String'))); 

  dyn_frames_init = [num2str(num_frames), ',', get(handles.t_sampling, 'String')]; 
  set(handles.dynamic_frames, 'String', dyn_frames_init);


% --- Executes during object creation, after setting all properties.
function t_sampling_CreateFcn(hObject, eventdata, handles)

  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
      set(hObject,'BackgroundColor','white');
  end


% --- Executes on mouse press over axes background.
function mip_image_ButtonDownFcn(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function mip_image_CreateFcn(hObject, eventdata, handles)


% --- Executes on selection change in roi_type.
function roi_type_Callback(hObject, eventdata, handles)

  guidata(hObject, handles); 


% --- Executes during object creation, after setting all properties.
function roi_type_CreateFcn(hObject, eventdata, handles)

  if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
      set(hObject,'BackgroundColor','white');
  end


% --- Executes on button press in plot_cod.
function plot_cod_Callback(hObject, eventdata, handles)

  % update last ROIs
  handles.r1_pos_all(handles.num_roi,:) = getPosition(handles.r1);
  handles.r2_pos_all(handles.num_roi,:) = getPosition(handles.r2);


  guidata(hObject, handles); 

  cod_temp = get_ROI_COD(handles);

  handles.roi_cod(handles.num_roi).x = cod_temp(1,:);
  handles.roi_cod(handles.num_roi).y = cod_temp(2,:);
  handles.roi_cod(handles.num_roi).z = cod_temp(3,:);

  guidata(hObject, handles); 

  handles.roi_state = get_roi_status(handles); 
  handles.motion_frames = []; 
  guidata(hObject, handles); 
  plot_cod_xyz_gui(handles);
  %set(handles.detect_motion, 'Enable', 'off'); 
  guidata(hObject, handles); 


% --- Executes on button press in detect_motion.
function detect_motion_Callback(hObject, eventdata, handles)

  handles.roi_state = get_roi_status(handles);
  guidata(hObject, handles); 
  handles.motion_frames = detect_motion(handles);
  guidata(hObject, handles); 

  plot_cod_xyz_gui(handles);
  guidata(hObject, handles);


% --- Executes on button press in roi_check_all.
function roi_check_all_Callback(hObject, eventdata, handles)

  all_val = get(handles.roi_check_all, 'Value');

  if all_val > 0.5
      if strcmp(handles.roi_check_1.Enable, 'on')
          set(handles.roi_check_1, 'Value', 1);
      end
      if strcmp(handles.roi_check_2.Enable, 'on')
          set(handles.roi_check_2, 'Value', 1);
      end
      if strcmp(handles.roi_check_3.Enable, 'on')
          set(handles.roi_check_3, 'Value', 1);
      end
      if strcmp(handles.roi_check_4.Enable, 'on')
          set(handles.roi_check_4, 'Value', 1);
      end
      if strcmp(handles.roi_check_5.Enable, 'on')
          set(handles.roi_check_5, 'Value', 1);
      end
      if strcmp(handles.roi_check_6.Enable, 'on')
          set(handles.roi_check_6, 'Value', 1);
      end
      if strcmp(handles.roi_check_7.Enable, 'on')
          set(handles.roi_check_7, 'Value', 1);
      end
      if strcmp(handles.roi_check_8.Enable, 'on')
          set(handles.roi_check_8, 'Value', 1);
      end
  end

  if all_val < 0.5
      if strcmp(handles.roi_check_1.Enable, 'on')
          set(handles.roi_check_1, 'Value', 0);
      end
      if strcmp(handles.roi_check_2.Enable, 'on')
          set(handles.roi_check_2, 'Value', 0);
      end
      if strcmp(handles.roi_check_3.Enable, 'on')
          set(handles.roi_check_3, 'Value', 0);
      end
      if strcmp(handles.roi_check_4.Enable, 'on')
          set(handles.roi_check_4, 'Value', 0);
      end
      if strcmp(handles.roi_check_5.Enable, 'on')
          set(handles.roi_check_5, 'Value', 0);
      end
      if strcmp(handles.roi_check_6.Enable, 'on')
          set(handles.roi_check_6, 'Value', 0);
      end
      if strcmp(handles.roi_check_7.Enable, 'on')
          set(handles.roi_check_7, 'Value', 0);
      end
      if strcmp(handles.roi_check_8.Enable, 'on')
          set(handles.roi_check_8, 'Value', 0);
      end
  end

  guidata(hObject, handles); 
    


% --- Executes on button press in roi_check_1.
function roi_check_1_Callback(hObject, eventdata, handles)


% --- Executes on button press in roi_check_2.
function roi_check_2_Callback(hObject, eventdata, handles)


% --- Executes on button press in roi_check_3.
function roi_check_3_Callback(hObject, eventdata, handles)


% --- Executes on button press in roi_check_4.
function roi_check_4_Callback(hObject, eventdata, handles)


% --- Executes on button press in roi_check_5.
function roi_check_5_Callback(hObject, eventdata, handles)


% --- Executes on button press in roi_check_6.
function roi_check_6_Callback(hObject, eventdata, handles)


% --- Executes on button press in roi_check_7.
function roi_check_7_Callback(hObject, eventdata, handles)


% --- Executes on button press in roi_check_8.
function roi_check_8_Callback(hObject, eventdata, handles)




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    **********     USER DEFINED FUNCTIONS     **************    %

function update_cod_xyz(handles)  

  if handles.num_roi < 0.5
    %errordlg('Must specify at least 1 ROI'); 
    %return
  end

  %cod_xyz_all = zeros(handles.num_roi, 1); 
  %cod_xyz_counts_all = zeros(handles.num_roi, 1);

  for n = 1:8
    fname_cod_xyz = [handles.outfolder_path, 'cod_xyz.',num2str(n),'.cg']; 
    fid1 = fopen(fname_cod_xyz, 'r'); 
    cod_xyz_temp = fread(fid1, inf, 'double'); 	
    fclose(fid1);  
    
    if n == 1
      cod_xyz_all = cod_xyz_temp; 
    else
      cod_xyz_all = cod_xyz_all + cod_xyz_temp; 
    end
    
    delete(fname_cod_xyz); 

  end

  fname_cod_xyz = [handles.outfolder_path, 'cod_xyz.cg'];
  fid_cod_xyz = fopen(fname_cod_xyz, 'w'); 
  fwrite(fid_cod_xyz,cod_xyz_all, 'double'); 
  fclose(fid_cod_xyz);



function cod_xyz = get_ROI_COD(handles)

  p1 = handles.r1_pos_all(handles.num_roi,:); 
  p2 = handles.r2_pos_all(handles.num_roi,:);
      
  x_span = [ceil(p1(2)), floor(p1(2)+p1(4))]
  y_span = [ceil(p2(2))-256, floor(p2(2)+p2(4)) - 256]
  z_span = [ceil(p1(1)), floor(p1(1)+p1(3))]

  x_span(x_span < 1) = 1; 
  y_span(y_span < 1) = 1;
  z_span(z_span < 1) = 1; 


  %roi_img = zeros(size(handles.mask_img)); 
  %roi_img(y_span(1):y_span(2), x_span(1):x_span(2), z_span(1):z_span(2)) = 1; 


  roi_img = handles.mask_img(y_span(1):y_span(2), x_span(1):x_span(2), z_span(1):z_span(2)); 

  %roi_img = roi_img .* handles.mask_img; 
  roi_img(roi_img < 0) = 0; 

  %roi_list = roi_img(roi_img > 0.5); 
  %roi_list = roi_list(:); 


  rois = sort(unique(roi_img(:)), 'ascend');

  if rois(1) == 0
    rois = rois(2:end); 
  end

  wts = ones(size(rois)); 
  for k = 1:length(rois)
    roi_counts_temp = handles.cod_xyz_counts(rois(k):handles.num_patch:end); 
    if k == 1
      roi_counts_all = roi_counts_temp; 
    else
      roi_counts_all = roi_counts_all + roi_counts_temp; 
    end
     
  end


  for k = 1:length(rois)
    cod_temp = handles.cod_xyz(:, rois(k):handles.num_patch:end); 
    counts_temp = handles.cod_xyz_counts(:, rois(k):handles.num_patch:end); 
    wt = counts_temp ./ roi_counts_all; 

    if k == 1
      cod_xyz = cod_temp .* wt; 
    else
      cod_xyz = cod_xyz + cod_temp .* wt; 
    end
  end


function plot_cod_xyz_gui(handles)

  ts = str2double(get(handles.t_sampling, 'String'));  

  % X
  maxm = -10000; 
  minm = 100000;
  axes(handles.x_cod);
  cla reset;
  hold on
  for k = 1:handles.num_roi

      cod_xyz_u0 = handles.roi_cod(k).x; 

      t_plot = 0:(size(cod_xyz_u0, 2));
      t_plot = ts .* t_plot(1:(end-1)); 
      
      if handles.roi_state(k) > 0.5
          plot(t_plot', cod_xyz_u0, handles.cc{k});
          m1 = max(cod_xyz_u0(:)); 
          m2 = min(cod_xyz_u0(:)); 
          if m1 > maxm
              maxm = m1; 
          end
          if m2 < minm
              minm = m2; 
          end
      end

  end
  yl = [minm-3 maxm+3];
  ylim(yl); 

  if ~isempty(handles.motion_frames)
      for j = 1:length(handles.motion_frames)
          plot([handles.motion_frames(j) handles.motion_frames(j)], [yl(1) yl(2)], 'k'); 
      end
  end  
      
  xlabel('sec'); 
  ylabel('mm');
  hold off


  % Y
  maxm = -10000; 
  minm = 100000;
  axes(handles.y_cod);
  cla reset;
  hold on
  for k = 1:handles.num_roi

      cod_xyz_u0 = handles.roi_cod(k).y; 

      t_plot = 0:(size(cod_xyz_u0, 2));
      t_plot = ts .* t_plot(1:(end-1)); 
      
      if handles.roi_state(k) > 0.5
          plot(t_plot', cod_xyz_u0, handles.cc{k});
          m1 = max(cod_xyz_u0(:)); 
          m2 = min(cod_xyz_u0(:)); 
          if m1 > maxm
              maxm = m1; 
          end
          if m2 < minm
              minm = m2; 
          end
      end

  end
  yl = [minm-3 maxm+3];
  ylim(yl); 

  if ~isempty(handles.motion_frames)
      for j = 1:length(handles.motion_frames)
          plot([handles.motion_frames(j) handles.motion_frames(j)], [yl(1) yl(2)], 'k'); 
      end
  end  

  xlabel('sec'); 
  ylabel('mm');
  hold off


  % Z
  maxm = -10000; 
  minm = 100000;
  axes(handles.z_cod);
  cla reset;
  hold on
  for k = 1:handles.num_roi

      cod_xyz_u0 = handles.roi_cod(k).z; 

      t_plot = 0:(size(cod_xyz_u0, 2));
      t_plot = ts .* t_plot(1:(end-1));
      
      if handles.roi_state(k) > 0.5
          plot(t_plot', cod_xyz_u0, handles.cc{k});
          m1 = max(cod_xyz_u0(:)); 
          m2 = min(cod_xyz_u0(:)); 
          if m1 > maxm
              maxm = m1; 
          end
          if m2 < minm
              minm = m2; 
          end
      end
      
  end
  yl = [minm-3 maxm+3];
  ylim(yl); 

  if ~isempty(handles.motion_frames)
      for j = 1:length(handles.motion_frames)
          plot([handles.motion_frames(j) handles.motion_frames(j)], [yl(1) yl(2)], 'k'); 
      end
  end  

  xlabel('sec'); 
  ylabel('mm');
  hold off

  guidata(gcf, handles); 



function motion_frames = detect_motion(handles)

  roi_choose = find(handles.roi_state > 0.5, 1, 'last'); 

  XY(1,:) = handles.roi_cod(roi_choose).x;
  XY(2,:) = handles.roi_cod(roi_choose).y;
  XY(3,:) = handles.roi_cod(roi_choose).z; 

  num_frame_pts = round(handles.min_frame_length/str2double(get(handles.t_sampling, 'String')));
  num_base_pts = round(handles.base_length/str2double(get(handles.t_sampling, 'String')));



  % select which roi's to use
  %for kn = 1:handles.num_roi
    %if handles.roi_state(kn) > 0.5
      %XY = dd(:,kn:handles.num_roi:end); 
      xtemp = XY(1,:); 
      ytemp = XY(2,:); 
      ztemp = XY(3,:); 
        
      mean_baseline_x = mean(xtemp(1:num_base_pts)); 
      mean_baseline_y = mean(ytemp(1:num_base_pts)); 
      mean_baseline_z = mean(ztemp(1:num_base_pts)); 

      std_baseline_x = std(xtemp(1:num_base_pts)); 
      std_baseline_y = std(ytemp(1:num_base_pts)); 
      std_baseline_z = std(ztemp(1:num_base_pts));

      thr = 2.5 * max([std_baseline_x, std_baseline_y, std_baseline_z]);

      xtemp = xtemp - mean_baseline_x;
      ytemp = ytemp - mean_baseline_y; 
      ztemp = ztemp - mean_baseline_z; 

      if size(xtemp, 2) > size(xtemp, 1)
        xtemp = xtemp';  
        ytemp = ytemp'; 
        ztemp = ztemp';
      end
        
      % do motion detect algorithm for each roi
      order = 4; 
      framel = 11; 

      %xtemp2 = sgolayfilt(xtemp, order, framel);  
      %ytemp2 = sgolayfilt(ytemp, order, framel); 
      %ztemp2 = sgolayfilt(ztemp, order, framel);

      xtemp2 = xtemp; 
      ytemp2 = ytemp; 
      ztemp2 = ztemp;  
          
      inds_x = abs(xtemp2((num_frame_pts+2):end) - xtemp2(1:(end-num_frame_pts-1))) > thr;
      inds_y = abs(ytemp2((num_frame_pts+2):end) - ytemp2(1:(end-num_frame_pts-1))) > thr;
      inds_z = abs(ztemp2((num_frame_pts+2):end) - ztemp2(1:(end-num_frame_pts-1))) > thr;

      inds_all = inds_x + inds_y + inds_z;
      inds_all(inds_all > 1.5) = 1;
          
      % find groups of 3-5 nearly consecutive numbers,  indicates consistent jump
      inds_new = [];
      for ki = 1:(length(inds_all)-num_frame_pts)
        temp_sum = sum(inds_all(ki:(ki+num_frame_pts-1)));
        if temp_sum > (num_frame_pts-1.5)
          inds_new = [inds_new, ki];
        end
      end

      inds_new = inds_new + num_frame_pts;
      inds_new = [1, inds_new];
      inds_new2 = 1;

      for kk = 2:(length(inds_new))

        if inds_new(kk-1) > (inds_new(kk) - 2)
          %%%%
        else
          inds_new2 = [inds_new2, inds_new(kk)];
        end
      end
      inds_new2 =  [inds_new2,  length(xtemp2)];  
      inds_all = inds_new2;

      frame_length = str2double(get(handles.t_sampling, 'String')) .* (inds_all(2:end) - inds_all(1:(end-1))); 

      num_motion_frames = length(inds_all) - 1;

      %inds_all = inds_all + 2; 
      motion_frames_init = inds_all .* str2double(get(handles.t_sampling, 'String')); 
        
    %end
  %end

  motion_frames = motion_frames_init


function mip_256x830 = interp_mip(mip_x)
  % interpolates 2D mip to 256 x 830
  % original mip must be in transaxial (row) x axial (column) image matrix format

  % try a few common cases for faster performance
  if size(mip_x, 1) == 512 && size(mip_x, 2) == 1660 % high resolution
    mip_256x830 = (mip_x(1:2:end, :) + mip_x(2:2:end, :)) ./ 2; 
    mip_256x830 = (mip_256x830(:, 1:2:end) + mip_256x830(:, 2:2:end)) ./ 2; 

  elseif size(mip_x, 1) == 512 && size(mip_x, 2) == 830 % high resolution transaxial
    mip_256x830 = (mip_x(1:2:end, :) + mip_x(2:2:end, :)) ./ 2; 

  else
    igrid_oldx = repmat(linspace(0, 1, size(mip_x, 2)), size(mip_x, 1), 1); 
    igrid_oldy = repmat((linspace(0, 1, size(mip_x, 1)))', 1, size(mip_x, 2)); 
    igrid_newx = repmat(linspace(0, 1, 830), 256, 1); 
    igrid_newy = repmat((linspace(0, 1, 256))', 1, 830); 

    mip_256x830 = interp2(igrid_oldx, igrid_oldy, mip_x, igrid_newx, igrid_newy); 

  end



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


function check_state = get_roi_status(handles)

  check_state = zeros(8,1); 

  check_state(1) = get(handles.roi_check_1, 'Value');
  check_state(2) = get(handles.roi_check_2, 'Value');
  check_state(3) = get(handles.roi_check_3, 'Value');
  check_state(4) = get(handles.roi_check_4, 'Value');
  check_state(5) = get(handles.roi_check_5, 'Value');
  check_state(6) = get(handles.roi_check_6, 'Value');
  check_state(7) = get(handles.roi_check_7, 'Value');
  check_state(8) = get(handles.roi_check_8, 'Value');


function myfun1(p)

  handles = guidata(gcf);
  handles.r1_pos = p; 
  handles.r2_pos = [p(1) handles.r2_pos(2) p(3) handles.r2_pos(4)]; 
  setPosition(handles.r2, handles.r2_pos); 
  guidata(gcf, handles); 


function myfun2(p)

  handles = guidata(gcf);
  handles.r2_pos = p; 
  handles.r1_pos = [p(1) handles.r1_pos(2) p(3) handles.r1_pos(4)]; 
  setPosition(handles.r1, handles.r1_pos); 
  guidata(gcf, handles); 


% --- Executes on button press in reset.
function handles = reset_Callback(hObject, eventdata, handles)

  disp('Reset'); 

  axes(handles.x_cod);
  cla reset;

  axes(handles.y_cod);
  cla reset;

  axes(handles.z_cod);
  cla reset;

  axes(handles.mip_image);
  cla reset;
  set(gca, 'box', 'on', 'xticklabel', {''}, 'yticklabel',  {''}, 'ticklength',  [0 0]);


  set(handles.plot_cod, 'Enable', 'off'); 
  set(handles.detect_motion, 'Enable', 'off'); 


  % system variables
  handles.num_crystals_all = 564480; 
  handles.num_threads = 8; 
  handles.run_recon = '0'; 
  handles.write_lm = '0'; 
  handles.write_sino = '0';
  handles.write_sino_block = '0';  
  handles.write_histo_img = '0'; 
  handles.num_iter = 5; 
  handles.cc = {'r', 'g', 'b', 'c', 'm', 'y', 'k', 'r', 'g', 'b', 'c', 'm', 'k'}; 

  % detect motion variables
  handles.min_frame_length = 10;  
  handles.thr_mm = 1.0; 
  handles.base_length = 30; 
  handles.motion_frames = []; 
  handles.num_roi = 0;
  handles.roi_state = zeros(8,1); 
  handles.motion_frames = []; 
  handles.roi_cod = []; 
  t_sampling_init = 1;  
  set(handles.t_sampling, 'String', num2str(t_sampling_init)); 

  set(handles.outfolder_name, 'String', ''); 
  set(handles.dynamic_frames, 'String', ''); 

  guidata(hObject, handles); 


  set(handles.add_roi, 'enable', 'off'); 

  set(handles.roi_check_all, 'Value', 0);
  set(handles.roi_check_all, 'Enable', 'off');

  set(handles.roi_check_1, 'Value', 0);
  set(handles.roi_check_1, 'Enable', 'off');

  set(handles.roi_check_2, 'Value', 0);
  set(handles.roi_check_2, 'Enable', 'off');

  set(handles.roi_check_3, 'Value', 0);
  set(handles.roi_check_3, 'Enable', 'off');

  set(handles.roi_check_4, 'Value', 0);
  set(handles.roi_check_4, 'Enable', 'off');

  set(handles.roi_check_5, 'Value', 0);
  set(handles.roi_check_5, 'Enable', 'off');

  set(handles.roi_check_6, 'Value', 0);
  set(handles.roi_check_6, 'Enable', 'off');

  set(handles.roi_check_7, 'Value', 0);
  set(handles.roi_check_7, 'Enable', 'off');

  set(handles.roi_check_8, 'Value', 0);
  set(handles.roi_check_8, 'Enable', 'off');

  guidata(hObject, handles); 

