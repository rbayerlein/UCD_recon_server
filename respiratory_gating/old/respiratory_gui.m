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

% Last Modified by GUIDE v2.5 20-Sep-2019 12:55:33

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
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


set(handles.outfolder_name, 'String', ''); 
set(handles.dynamic_frames, 'String', ''); 

%handles.path_choose = reset; 
%handles.fname_choose_base = reset; 




% --- Executes on button press in choose_listmode.
function choose_listmode_Callback(hObject, eventdata, handles)
% hObject    handle to choose_listmode (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[fname_choose, pathname_choose] = uigetfile('.raw')

ind_ext = strfind(fname_choose, '.raw') - 2; 
fname_choose_base = fname_choose(1:ind_ext)


handles.fname_choose_base = fname_choose_base; 
handles.path_choose = pathname_choose; 

% get os type
if strcmp(pathname_choose(1), '/') > 0.5
  handles.os_type = 'unix'; 
  disp('unix')
else
  handles.os_type = 'windows'; 
  disp('windows')
end


handles.num_threads = 8; 

handles.run_recon = '0'; 
handles.write_lm = '0'; 
handles.write_sino = '0';
handles.write_sino_block = '1';  





guidata(hObject, handles);







function outfolder_name_Callback(hObject, eventdata, handles)
% hObject    handle to outfolder_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outfolder_name as text
%        str2double(get(hObject,'String')) returns contents of outfolder_name as a double



% --- Executes during object creation, after setting all properties.
function outfolder_name_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outfolder_name (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in start_pushbutton.
function start_pushbutton_Callback(hObject, eventdata, handles)
% hObject    handle to start_pushbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

guidata(hObject, handles);

dyn_frames = get(handles.dynamic_frames, 'String'); 
[val, dyn_frames_nospace] = check_dyn_frames(dyn_frames); 

if val < 0.5
    errordlg('Invalid characters in dynamic frames input, please check'); 
    return
end 


num_frames = getNumFrames(dyn_frames_nospace)


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

  sh_str = []; 
  for fnum = 1:8
   
    sh_str_temp = ''; 
   
    infile_fname = [handles.fname_choose_base, num2str(fnum), '.raw']; 
    infile_fullpath = fullfile(handles.path_choose, infile_fname); 
    %infile_fullpath = ['"',infile_fullpath,'"']; 
   
    
    fname_config = ['/run/media/meduser/data/explorer-master_24092019/read_lm/Reconstruction_Parameters_',num2str(fnum)];
    fid_config = fopen(fname_config, 'w'); 
   
    fprintf(fid_config, '%s\n%s\n%s\n%s\n%s\n%s\n%s', outfolder, infile_fullpath, handles.run_recon, dyn_frames_nospace, handles.write_lm, handles.write_sino, handles.write_sino_block); 
   
    fclose(fid_config); 
   
    if fnum == -1 || fnum == -1
      sh_str_temp = ['/run/media/meduser/data/explorer-master_24092019/read_lm/process_lm_sino_explorer_singlefile  ',fname_config, ' ; ']; 
    else
      sh_str_temp = ['/run/media/meduser/data/explorer-master_24092019/read_lm/process_lm_sino_explorer_singlefile  ',fname_config, ' & '];
    end
    
    sh_str = [sh_str, sh_str_temp]; 
 
    
  end
  
  sh_str = [sh_str, ' wait; ']; 
  
  sh_str = [sh_str, '/run/media/meduser/data/explorer-master_24092019/read_lm/sum_sino ']; 
  for kk = 1:8
  	ff = [outfolder, 'sinogram_block_f0_prompts.',num2str(kk),'.raw ']; 
  	sh_str = [sh_str, ff]; 
  end
  
  sh_str = [sh_str, ' & ']; 
  
  sh_str = [sh_str, '/run/media/meduser/data/explorer-master_24092019/read_lm/sum_sino ']; 
  for kk = 1:8
  	ff = [outfolder, 'sinogram_block_f0_randoms.',num2str(kk),'.raw ']; 
  	sh_str = [sh_str, ff]; 
  end
  
  sh_str = [sh_str, ' & wait']; 
  
  

  system(sh_str);

end



%sum_sino(outfolder, 0); 






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

% str_e = str_n;
% numstr = frame_str(str_b:str_e); 
% fnum = str2double(numstr); 
% nframes = nframes + fnum;
% nframes = round(nframes); 
