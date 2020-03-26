function varargout = EOSoutline(varargin)
% EOSOUTLINE MATLAB code for EOSoutline.fig
%      EOSOUTLINE, by itself, creates a new EOSOUTLINE or raises the existing
%      singleton*.
%
%      H = EOSOUTLINE returns the handle to a new EOSOUTLINE or the handle to
%      the existing singleton*.
%
%      EOSOUTLINE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EOSOUTLINE.M with the given input arguments.
%
%      EOSOUTLINE('Property','Value',...) creates a new EOSOUTLINE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EOSoutline_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EOSoutline_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EOSoutline

% Last Modified by GUIDE v2.5 26-Mar-2020 14:26:17

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EOSoutline_OpeningFcn, ...
                   'gui_OutputFcn',  @EOSoutline_OutputFcn, ...
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


% --- Executes just before EOSoutline is made visible.
function EOSoutline_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EOSoutline (see VARARGIN)

 % get mat file detailing scan directory
[fn, fp] = uigetfile('*.mat', 'Select Scan Directory File'); 
load([fp, '\', fn]); 
handles.base_fp = base_fp;
handles.img_fp = img_fp; 
handles.patient_list = patient_list;

handles.current_patient = 1;

% Choose default command line output for EOSoutline
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

shownewpatient(hObject, eventdata, handles)

% UIWAIT makes EOSoutline wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EOSoutline_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in inputCor.
function inputCor_Callback(hObject, eventdata, handles)
% hObject    handle to inputCor (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axesCor); 
OL = roipoly; 

% --- Executes on button press in inputSag.
function inputSag_Callback(hObject, eventdata, handles)
% hObject    handle to inputSag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in show3D.
function show3D_Callback(hObject, eventdata, handles)
% hObject    handle to show3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in nextButton.
function nextButton_Callback(hObject, eventdata, handles)
% hObject    handle to nextButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current_patient = handles.current_patient+1;
if current_patient <= length(handles.patient_list)
    handles.current_patient = current_patient;
    guidata(hObject, handles);
    shownewpatient(hObject, eventdata, handles);
end

% --- Executes on button press in prevButton.
function prevButton_Callback(hObject, eventdata, handles)
% hObject    handle to prevButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
current_patient = handles.current_patient-1;
if current_patient >= 1
    handles.current_patient = current_patient;
    guidata(hObject, handles);
    shownewpatient(hObject, eventdata, handles);
end

function shownewpatient(hObject, eventdata, handles)
% get DICOM info 
ifoCor = dicomread([handles.base_fp,...
    num2str(handles.patient_list(handles.current_patient)),...
    handles.img_fp,...
    'cor']);
ifoSag = dicomread([handles.base_fp,...
    num2str(handles.patient_list(handles.current_patient)),...
    handles.img_fp,...
    'sag']);
handles.ifoCor = ifoCor; handles.ifoSag = ifoSag; 
% get DICOM img
imgCor = dicomread([handles.base_fp,...
    num2str(handles.patient_list(handles.current_patient)),...
    handles.img_fp,...
    'cor']);
imgSag = dicomread([handles.base_fp,...
    num2str(handles.patient_list(handles.current_patient)),...
    handles.img_fp,...
    'sag']);
% convert int to double [0,1]

% show/store 
handles.imgCor = imgCor; handles.imgSag = imgSag; 
axes(handles.axesCor); imshow(imgCor, []); 
axes(handles.axesSag); imshow(imgSag, []);

% Update handles structure
guidata(hObject, handles);
