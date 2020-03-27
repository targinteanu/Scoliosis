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

% Last Modified by GUIDE v2.5 27-Mar-2020 17:37:18

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
hold off; imshow(handles.imgCor); hold on;
OL = roipoly; 
[handles.imgCorFilt, handles.splCorObj, handles.splCorSmp] = process(handles.imgCor, OL); 
guidata(hObject, handles);

% --- Executes on button press in inputSag.
function inputSag_Callback(hObject, eventdata, handles)
% hObject    handle to inputSag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axesSag); 
hold off; imshow(handles.imgSag); hold on;
OL = roipoly; 
[handles.imgSagFilt, handles.splSagObj, handles.splSagSmp] = process(handles.imgSag, OL); 
guidata(hObject, handles);

% --- Executes on button press in show3D.
function show3D_Callback(hObject, eventdata, handles)
% hObject    handle to show3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% set dimensions to mm
splCorScl = splineRescale(handles.splCorObj, ...
    handles.ifoCor.PixelSpacing(1), handles.ifoCor.PixelSpacing(2));
splSagScl = splineRescale(handles.splSagObj, ...
    handles.ifoSag.PixelSpacing(1), handles.ifoSag.PixelSpacing(2));
zCor = handles.splCorSmp(:,1)*handles.ifoCor.PixelSpacing(1); 
zSag = handles.splSagSmp(:,1)*handles.ifoSag.PixelSpacing(1);
zScl = max(zCor(1),zSag(1)):.5:min(zCor(end),zSag(end)); 
xScl = ppval(zScl, splSagScl.p);
yScl = ppval(zScl, splCorScl.p);

% store 
handles.splSclObj = {splSagScl, splCorScl};
handles.splSclRng = [min(zScl); max(zScl)];
handles.splSclSmp = [xScl, yScl, zScl];
guidata(hObject, handles);

% display 
axes(handles.axes3); plot3(xScl, yScl, zScl, 'b'); grid on; 

function splObj = splineRescale(splObj, zscl, xscl)
splPP = splObj.p.coefs;
zmult = (1/zscl).^(0:(size(splPP, 2)-1));
splPP = splPP.*zmult; 
splPP = splPP * xscl; 
splObj.p.coefs = splPP; 
splObj.p.breaks = splObj.p.breaks*zscl;

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
ifoCor = dicominfo([handles.base_fp,...
    num2str(handles.patient_list(handles.current_patient)),...
    handles.img_fp,...
    'cor']);
ifoSag = dicominfo([handles.base_fp,...
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
% convert uint to double [0,1]
imgCor = double(imgCor) / double(2^ifoCor.BitsStored - 1);
imgSag = double(imgSag) / double(2^ifoSag.BitsStored - 1);
% show/store 
handles.imgCor = imgCor; handles.imgSag = imgSag; 
axes(handles.axesCor); hold off; imshow(imgCor); 
axes(handles.axesSag); hold off; imshow(imgSag);

% Update handles structure
guidata(hObject, handles);

function [imgFiltered, splineObj, splineSample] = process(img, outln)
% img is a double img [0, 1]; outln is a BW outline

% filtering 
im = imreconstruct(double(~outln), img); 
imgFiltered = img - im; 
imgFiltered = imgFiltered .* outln;

% spline fitting 
[r, c] = find(imgFiltered); 
[splineObj, gof] = fit(r, c, 'smoothingspline', 'Weights', imgFiltered(find(imgFiltered(:))));
z = min(r):max(r); 
splineSample = [z; ppval(z, splineObj.p)]';

% display results on current axes 
visboundaries(outln); 
plot(splineSample(:,2), splineSample(:,1), 'b');