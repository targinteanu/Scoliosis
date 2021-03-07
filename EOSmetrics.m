function varargout = EOSmetrics(varargin)
% EOSMETRICS MATLAB code for EOSmetrics.fig
%      EOSMETRICS, by itself, creates a new EOSMETRICS or raises the existing
%      singleton*.
%
%      H = EOSMETRICS returns the handle to a new EOSMETRICS or the handle to
%      the existing singleton*.
%
%      EOSMETRICS('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in EOSMETRICS.M with the given input arguments.
%
%      EOSMETRICS('Property','Value',...) creates a new EOSMETRICS or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before EOSmetrics_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to EOSmetrics_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help EOSmetrics

% Last Modified by GUIDE v2.5 04-Mar-2021 23:34:47

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @EOSmetrics_OpeningFcn, ...
                   'gui_OutputFcn',  @EOSmetrics_OutputFcn, ...
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


% --- Executes just before EOSmetrics is made visible.
function EOSmetrics_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to EOSmetrics (see VARARGIN)

 % get mat file detailing scan directory
[fn, fp] = uigetfile('*.mat', 'Select Scan Directory File'); 
load([fp, '\', fn]); 
handles.base_fp = base_fp;
handles.img_fp = img_fp; 
handles.patient_list = patient_list;

handles.current_patient = 1;

update3Dview(eventdata, handles);

handles.metricSelecter.String = {...
    'Curve Only', ...
    'Numeric X Derivative', ...
    'Numeric Y Derivative', ...
    'Numeric X 2nd Derivative', ...
    'Numeric Y 2nd Derivative', ...
    'Numeric Second Derivative', ...
    'Numeric Curvature', ...
    'Numeric Torsion', ...
    'Lewiner 1st Derivative', ...
    'Lewiner 2nd Derivative', ...
    'Lewiner 3rd Derivative', ...
    'Lewiner Curvature', ...
    'Lewiner Torsion'};

handles.metricFuncs = {...
    @(HO,ED,H) showpatient(HO,ED,H), ...
    @(HO,ED,H) showderivD(HO,ED,H,1), ...
    @(HO,ED,H) showderivD(HO,ED,H,2), ...
    @(HO,ED,H) showderiv2D(HO,ED,H,1), ...
    @(HO,ED,H) showderiv2D(HO,ED,H,2), ...
    @(HO,ED,H) showderiv2(HO,ED,H), ...
    @(HO,ED,H) showNumericCurvature(HO,ED,H), ...
    @(HO,ED,H) showNumericTorsion(HO,ED,H), ...
    @(HO,ED,H) showLewinerQuantity(HO,ED,H,1), ...
    @(HO,ED,H) showLewinerQuantity(HO,ED,H,2), ...
    @(HO,ED,H) showLewinerQuantity(HO,ED,H,3), ...
    @(HO,ED,H) showLewinerQuantity(HO,ED,H,5), ...
    @(HO,ED,H) showLewinerQuantity(HO,ED,H,4)};

% Choose default command line output for EOSmetrics
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

shownewpatient(hObject, eventdata, handles)

% UIWAIT makes EOSmetrics wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = EOSmetrics_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in pushbutton3.
function pushbutton3_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton3 (see GCBO)
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


function showLinearHeatmap(hObject, eventdata, handles, XYZ, heatvals)
axes(handles.axes3);
x = XYZ(:,1)'; y = XYZ(:,2)'; z = XYZ(:,3)'; h = heatvals';
cla;
surface([x;x], [y;y], -[z;z], [h;h], ...
    'facecol', 'no', 'edgecol', 'interp', 'LineWidth', 5);
grid on;
colorbar;


function showLewinerQuantity(hObject, eventdata, handles, varToShow)
% varToShow key: 
%   1 - first derivative 
%   2 - second derivative 
%   3 - third derivative 
%   4 - torsion
%   5 - curvature 
q = inputdlg('Number of Fit Points (q):', 'q selection');
q = str2num(q{1});
p = handles.splfilt;
    vertebrae = (1+q):(size(p,1)-q); 
    tau = zeros(size(vertebrae)); % local torsion at each point on the spine 
    kappa = zeros(size(vertebrae)); % local curvature at each point on the spine 
    d2 = zeros(size(vertebrae)); d1 = d2; d3 = d2;
    for vertebra = vertebrae
        [tau(vertebra-q), d,dd,ddd, kappa(vertebra-q)] = ...
            lewinerTorsion(p, vertebra, q);
        d2(vertebra-q) = norm(dd); d1(vertebra-q) = norm(d); d3(vertebra-q) = norm(ddd);
    end
vars = {d1, d2, d3, tau, kappa}; var = vars{varToShow};
showLinearHeatmap(hObject, eventdata, handles, p(vertebrae,:), var');


function showNumericCurvature(hObject, eventdata, handles)
R = handles.splfilt; ddR = diff(diff(R));
Ri = .5 * ( R(2:end,:) + R(1:(end-1),:) );
dR = diff(Ri);
K = arrayfun(@(i) norm(cross(dR(i,:), ddR(i,:)))./( norm(dR(i,:)).^3 ), ...
    1:size(dR,1));
showLinearHeatmap(hObject, eventdata, handles, R(2:(end-1),:), K');

function showNumericTorsion(hObject, eventdata, handles)
R = handles.splfilt; dddR = diff(diff(diff(R)));
dR = diff(R(2:(end-1),:));
Ri = .5 * ( R(2:end,:) + R(1:(end-1),:) );
ddR = diff(diff(Ri));
T = arrayfun(@(i) ( dR(i,:) * cross( ddR(i,:), dddR(i,:))' ) ./ ...
                  ( norm( cross( dR(i,:), ddR(i,:) ) ).^2 ) , ...
             1:size(dR,1));
showLinearHeatmap(hObject, eventdata, handles, Ri(2:(end-1),:), T');


function showderivD(hObject, eventdata, handles, D)
R = handles.splfilt;
dR = diff(R);
deriv1 = dR(:,D)./dR(:,3);
Ri = .5 * ( R(2:end,:) + R(1:(end-1),:) );
showLinearHeatmap(hObject, eventdata, handles, Ri, deriv1)

function showderiv2D(hObject, eventdata, handles, D)
R = handles.splfilt;
dR = diff(R);
deriv1 = dR(:,D)./dR(:,3); ddR = diff(deriv1);
Ri = .5 * ( R(2:end,:) + R(1:(end-1),:) ); dRi = diff(Ri);
deriv2 = ddR./dRi(:,3);
showLinearHeatmap(hObject, eventdata, handles, R(2:(end-1),:), deriv2)

function showderiv2(hObject, eventdata, handles)
R = handles.splfilt;
dR = diff(R); ddR = diff(dR);
Ri = .5 * ( R(2:end,:) + R(1:(end-1),:) ); dRi = diff(Ri);
ds2 = sum(dRi.^2, 2);
deriv2 = ddR./ds2; deriv2 = sum(deriv2.^2, 2).^.5;
showLinearHeatmap(hObject, eventdata, handles, R(2:(end-1),:), deriv2)


function showpatient(hObject, eventdata, handles)
ifoCor = handles.ifoCor; ifoSag = handles.ifoSag;
imgCor = handles.imgCor; imgSag = handles.imgSag; 
axes(handles.axesCor); hold off; imshow(imgCor); 
axes(handles.axesSag); hold off; imshow(imgSag);

splSclObj = handles.splSclObj;
splSclRng = handles.splSclRng;
splSclSmp = handles.splSclSmp;
splSclBnd = handles.splSclBnd;

imgSagFilt = handles.imgSagFilt; 
splSagObj = handles.splSagObj; 
splSagSmp = handles.splSagSmp;
SagOL = handles.SagOL;
splSagObjBound = handles.splSagObjBound; 
splSagSmpBound = handles.splSagSmpBound;
imgCorFilt = handles.imgCorFilt; 
splCorObj = handles.splCorObj; 
splCorSmp = handles.splCorSmp;
CorOL = handles.CorOL;
splCorObjBound = handles.splCorObjBound; 
splCorSmpBound = handles.splCorSmpBound;

    axes(handles.axesCor); hold on; 
    plot(splCorSmp(:,2), splCorSmp(:,1), ':b');
    plot(splCorSmpBound(2), splCorSmpBound(1), '*b', 'LineWidth', 1.25);
    
    axes(handles.axesSag); hold on; 
    plot(splSagSmp(:,2), splSagSmp(:,1), ':b');
    plot(splSagSmpBound(2), splSagSmpBound(1), '*b', 'LineWidth', 1.25);
    
        splfilt = handles.splfilt;

        splfiltPix = splfilt/handles.ifoCor.PixelSpacing(1);
                
        axes(handles.axesSag); hold on;
        plot(splfiltPix(:,1), splfiltPix(:,3), 'm', 'LineWidth', 2);
        axes(handles.axesCor); hold on;
        plot(splfiltPix(:,2), splfiltPix(:,3), 'm', 'LineWidth', 2);
        
        axes(handles.axes3); cla; 
        plot3(splSclSmp(:,1), splSclSmp(:,2), -splSclSmp(:,3), ':b'); hold on;
        plot3(splfilt(:,1), splfilt(:,2), -splfilt(:,3), 'm', 'LineWidth', 2);
        grid on; update3Dview(eventdata, handles); 
        xlim([0, size(imgSag,2)]*ifoSag.PixelSpacing(2));
        ylim([0, size(imgCor,2)]*ifoCor.PixelSpacing(2));
        zlim(-[size(imgCor,1), 0]*ifoCor.PixelSpacing(1));

function shownewpatient(hObject, eventdata, handles)
set(handles.PatientNum, 'String', ...
    num2str(handles.patient_list(handles.current_patient)));
handles.metricSelecter.Value = 1;

% get DICOM info 
ifoCor = dicominfo([handles.base_fp,...
    num2str(handles.patient_list(handles.current_patient)),...
    handles.img_fp,...
    'cor']);
ifoSag = dicominfo([handles.base_fp,...
    num2str(handles.patient_list(handles.current_patient)),...
    handles.img_fp,...
    'sag']);
% PixelSpacing 
if ~exist('ifoCor.PixelSpacing', 'var')
    ifoCor.PixelSpacing = ifoCor.ImagerPixelSpacing;
end
if ~exist('ifoSag.PixelSpacing', 'var')
    ifoSag.PixelSpacing = ifoSag.ImagerPixelSpacing;
end
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

% use saved file if available
fn = [handles.base_fp,...
    num2str(handles.patient_list(handles.current_patient)),...
    handles.img_fp,...
    'patient',num2str(handles.patient_list(handles.current_patient)),...
    ' EOSoutline data.mat'];

fn2 = [handles.base_fp,...
    num2str(handles.patient_list(handles.current_patient)),...
    handles.img_fp,...
    'patient',num2str(handles.patient_list(handles.current_patient)),...
    ' filtered data.mat'];

if exist(fn, 'file')
    disp(['loading ',fn])
    load(fn);
    
handles.splSclObj = splSclObj;
handles.splSclRng = splSclRng;
handles.splSclSmp = splSclSmp;
handles.splSclBnd = splSclBnd;

handles.imgSagFilt = imgSagFilt; 
handles.splSagObj = splSagObj; 
handles.splSagSmp = splSagSmp;
handles.SagOL = SagOL;
handles.splSagObjBound = splSagObjBound; 
handles.splSagSmpBound = splSagSmpBound;
handles.imgCorFilt = imgCorFilt; 
handles.splCorObj = splCorObj; 
handles.splCorSmp = splCorSmp;
handles.CorOL = CorOL;
handles.splCorObjBound = splCorObjBound; 
handles.splCorSmpBound = splCorSmpBound;

    if exist(fn2, 'file')
        disp(['loading ',fn2])
        load(fn2);
        
        handles.splfilt = splfilt;
                
    end
    
else
clear handles.splSclObj handles.splSclRng handles.splSclSmp
clear handles.imgSagFilt handles.splSagObj handles.splSagSmp handles.SagOL
clear handles.imgCorFilt handles.splCorObj handles.splCorSmp handles.CorOL
end

showpatient(hObject, eventdata, handles);

% Update handles structure
guidata(hObject, handles);


% --- Executes on slider movement.
function slider1_Callback(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
update3Dview(eventdata, handles)


% --- Executes during object creation, after setting all properties.
function slider1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider2_Callback(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider
update3Dview(eventdata, handles)


% --- Executes during object creation, after setting all properties.
function slider2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


function update3Dview(eventdata, handles)
axes(handles.axes3);
view(handles.slider1.Value*360, handles.slider2.Value*360);


% --- Executes on selection change in metricSelecter.
function metricSelecter_Callback(hObject, eventdata, handles)
% hObject    handle to metricSelecter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns metricSelecter contents as cell array
%        contents{get(hObject,'Value')} returns selected item from metricSelecter
fcn = handles.metricFuncs{hObject.Value};
fcn(hObject, eventdata, handles);


% --- Executes during object creation, after setting all properties.
function metricSelecter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to metricSelecter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
