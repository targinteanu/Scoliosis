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

% Last Modified by GUIDE v2.5 08-Sep-2020 02:47:47

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

qstOL = questdlg('Outline Spine?', 'Outline', 'Yes', 'No', 'Yes');
qstOL = strcmp(qstOL, 'Yes');
if qstOL
    uiwait(msgbox(['Closely outline as much of the spine as possible. ' ...
        'Include at least C8-S1. Include vertebral bodies only.'], 'Instructions'));
    [OL, vertx, verty] = roipoly;
    handles.CorOL = OL; handles.CorVertX = vertx; handles.CorVertY = verty; 
end

qstBnd = questdlg('Separate C8 Bound?', 'Bound', 'Yes', 'No', 'Yes');
qstBnd = strcmp(qstBnd, 'Yes');
if qstBnd
    uiwait(msgbox(['Separate C8 from T1 with a straight line. ' ...
        'Click to set the two endpoints of the line. Double click when complete.'], ...
        'Instructions'));
    [boundX, boundY] = getpts;
    handles.CorBound = [boundX, boundY];
end

qstFH = questdlg('Locate Femoral Heads?', 'Femoral Heads', 'Yes', 'No', 'Yes');
qstFH = strcmp(qstFH, 'Yes');
if qstFH
    uiwait(msgbox(['Indicate a diameter line of the left femoral head. ' ...
        'Click to set the two endpoints of the line. Double click when complete. ' ...
        'Only the first two points will be considered. '], ...
        'Instructions'));
    [FH_L_X, FH_L_Y] = getpts; FH_L = [FH_L_X'; FH_L_Y']; 
    FH_L = FH_L(1:2,:);
    FH_L = mean(FH_L, 2);
    
    uiwait(msgbox(['Indicate a diameter line of the right femoral head. ' ...
        'Click to set the two endpoints of the line. Double click when complete. '...
        'Only the first two points will be considered. '], ...
        'Instructions'));
    [FH_R_X, FH_R_Y] = getpts; FH_R = [FH_R_X'; FH_R_Y']; 
    FH_R = FH_R(1:2,:);
    FH_R = mean(FH_R, 2);
    
    handles.femheadsCor = [FH_L, FH_R];
end

if qstOL | qstBnd
[handles.imgCorFilt, handles.splCorObj, handles.splCorSmp, ...
    handles.splCorObjBound, handles.splCorSmpBound] = ...
    processOL(handles.imgCor, OL, boundX, boundY, ...
    handles.ifoCor.PixelSpacing(1), handles.ifoCor.PixelSpacing(2)); 
end

if qstOL | qstBnd | qstFH
guidata(hObject, handles);
end

% --- Executes on button press in inputSag.
function inputSag_Callback(hObject, eventdata, handles)
% hObject    handle to inputSag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axesSag); 
hold off; imshow(handles.imgSag); hold on;

qstOL = questdlg('Outline Spine?', 'Outline', 'Yes', 'No', 'Yes');
qstOL = strcmp(qstOL, 'Yes');
if qstOL
    uiwait(msgbox(['Closely outline as much of the spine as possible. ' ...
        'Include at least C8-S1. Include vertebral bodies only.'], 'Instructions'));
    [OL, vertx, verty] = roipoly;
    handles.SagOL = OL; handles.SagVertX = vertx; handles.SagVertY = verty; 
end

qstBnd = questdlg('Separate L5 Bound?', 'Bound', 'Yes', 'No', 'Yes');
qstBnd = strcmp(qstBnd, 'Yes');
if qstBnd
    uiwait(msgbox(['Separate L5 from S1 with a straight line. ' ...
        'Click to set the two endpoints of the line. Double click when complete.'], ...
        'Instructions'));
    [boundX, boundY] = getpts;
    handles.SagBound = [boundX, boundY];
end

qstFH = questdlg('Locate Femoral Heads?', 'Femoral Heads', 'Yes', 'No', 'Yes');
qstFH = strcmp(qstFH, 'Yes');
if qstFH
    uiwait(msgbox(['Indicate a diameter line of the left femoral head. ' ...
        'Click to set the two endpoints of the line. Double click when complete. ' ...
        'Only the first two points will be considered. '], ...
        'Instructions'));
    [FH_L_X, FH_L_Y] = getpts; FH_L = [FH_L_X'; FH_L_Y']; 
    FH_L = FH_L(1:2,:);
    FH_L = mean(FH_L, 2);
    
    uiwait(msgbox(['Indicate a diameter line of the right femoral head. ' ...
        'Click to set the two endpoints of the line. Double click when complete. '...
        'Only the first two points will be considered. '], ...
        'Instructions'));
    [FH_R_X, FH_R_Y] = getpts; FH_R = [FH_R_X'; FH_R_Y']; 
    FH_R = FH_R(1:2,:);
    FH_R = mean(FH_R, 2);
    
    handles.femheadsSag = [FH_L, FH_R];
end

if qstOL | qstBnd
[handles.imgSagFilt, handles.splSagObj, handles.splSagSmp, ...
    handles.splSagObjBound, handles.splSagSmpBound] = ...
    processOL(handles.imgSag, OL, boundX, boundY, ...
    handles.ifoSag.PixelSpacing(1), handles.ifoSag.PixelSpacing(2)); 
end

if qstOL | qstBnd | qstFH
guidata(hObject, handles);
end

function [imgFiltered, splineObj, splineSample, splineObjBound, splineSampleBound] ...
    = processOL(img, outln, boundX, boundY, zscl, xscl)
% img is a double img [0, 1]; outln is a BW outline

% filtering 
im = imreconstruct(double(~outln), img); 
imgFiltered = img - im; 
imgFiltered = imgFiltered .* outln;

% spline fitting 
[r, c] = find(imgFiltered);
[splineObj, gof] = fit(r*zscl, c*xscl, 'smoothingspline', 'Weights', imgFiltered(find(imgFiltered(:))));
boundX = boundX*xscl; boundY = boundY*zscl;
[xi, yi] = lineSeg_PP_intersect(splineObj.p, [boundX(1), boundY(1)], [boundX(2), boundY(2)]); %xi, yi in mm
splineObjBound = [yi; xi];
z = min(r):max(r); 
z = z*zscl;
splineSample = [z/zscl; ppval(z, splineObj.p)/xscl]';
splineSampleBound = [yi/zscl; xi/xscl];
% splineSample is in pixels, but splineObj is in mm

% display results on current axes 
visboundaries(outln); 
plot(splineSample(:,2), splineSample(:,1), 'b');
plot(splineSampleBound(2), splineSampleBound(1), '*r');

function [xi, yi] = lineSeg_PP_intersect(pp, xy1, xy2)
x1 = xy1(1); y1 = xy1(2); x2 = xy2(1); y2 = xy2(2);
m = (y2-y1)./(x2-x1); b = y1-m*x1;
ypp = linspace(min(y1, y2), max(y1, y2));
xpp = ppval(pp, ypp); 
[~, i] = min(ypp - m.*xpp - b); xi = xpp(i); yi = ypp(i);


% --- Executes on button press in show3D.
function show3D_Callback(hObject, eventdata, handles)
% hObject    handle to show3D (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

splCorScl = handles.splCorObj; splSagScl = handles.splSagObj;
zCor = handles.splCorSmp(:,1)*handles.ifoCor.PixelSpacing(1); % mm
zSag = handles.splSagSmp(:,1)*handles.ifoSag.PixelSpacing(1); % mm
zScl = max(zCor(1),zSag(1)):.5:min(zCor(end),zSag(end)); 
xScl = ppval(zScl, splSagScl.p);
yScl = ppval(zScl, splCorScl.p);

botBoundScl = handles.splSagObjBound; 
topBoundScl = handles.splCorObjBound; 
[~,iBot] = min(sum(([zScl; xScl] - botBoundScl).^2));
[~,iTop] = min(sum(([zScl; yScl] - topBoundScl).^2));

fhXZsag = handles.femheadsSag./flipud(handles.ifoSag.PixelSpacing); % [x1, x2; z1, z2]
fhYZcor = handles.femheadsCor./flipud(handles.ifoCor.PixelSpacing); % [y1, y2; z1, z2]
fhXYZ = [fhXZsag(1,:), fhYZcor(1,:), mean([fhXZsag(2,:); fhYZcor(2,:)])];

% store 
handles.splSclObj = {splSagScl, splCorScl};
handles.splSclRng = [min(zScl); max(zScl)];
handles.splSclSmp = [xScl', yScl', zScl'];
handles.splSclBnd = [iBot, iTop];
handles.femheadsScl = fhXYZ; % [x1, x2; y1, y2; z1, z2]
guidata(hObject, handles);

% display 
axes(handles.axes3); hold off;
plot3(xScl, yScl, zScl, 'b'); grid on; hold on;
plot3([xScl(iBot),xScl(iTop)], [yScl(iBot),yScl(iTop)], [zScl(iBot),zScl(iTop)], '*r');
plot3(fhXYZ(1,:), fhXYZ(2,:), fhXYZ(3,:), 'or');
xlim([1, size(handles.imgSag,2)] * handles.ifoSag.PixelSpacing(2));
ylim([1, size(handles.imgCor,2)] * handles.ifoCor.PixelSpacing(2));
xlabel('x(mm)'); ylabel('y(mm)'); zlabel('z(mm)');

% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fn = [handles.base_fp,...
    num2str(handles.patient_list(handles.current_patient)),...
    handles.img_fp,...
    'patient',num2str(handles.patient_list(handles.current_patient)),...
    ' EOSoutline data.mat'];

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

femheadsCor = handles.femheadsCor;
femheadsSag = handles.femheadsSag; 
femheadsScl = handles.femheadsScl;

save(fn, 'splSclObj', 'splSclRng', 'splSclSmp', 'splSclBnd', ...
    'imgSagFilt', 'splSagObj', 'splSagSmp', 'SagOL', ...
    'imgCorFilt', 'splCorObj', 'splCorSmp', 'CorOL', ...
    'femheadsCor', 'femheadsSag', 'femheadsScl', ...
    'splCorObjBound', 'splCorSmpBound', 'splSagObjBound', 'splSagSmpBound');
disp(['saving ',fn])


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
%show3D_Callback(hObject, eventdata, handles)

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
%show3D_Callback(hObject, eventdata, handles)

function shownewpatient(hObject, eventdata, handles)
set(handles.PatientNum, 'String', ...
    num2str(handles.patient_list(handles.current_patient)));

set(handles.Writhe, 'String', 'Wr = ');

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
axes(handles.axesCor); hold off; imshow(imgCor); 
axes(handles.axesSag); hold off; imshow(imgSag);

% initialize vars that may not have been set
femheadsScl = [0,0;0,0;0,0];
femheadsCor = [0,0;0,0];
femheadsSag = [0,0;0,0];

splSclObj = []; % fix!!!
splSclRng = [0;0];
splSclSmp = [0,0,0];
splSclBnd = [0,0];

imgSagFilt = [];
splSagObj = [];
splSagSmp = [0,0];
SagOL = [];
splSagObjBound = [];
splSagSmpBound = [0,0];
imgCorFilt = [];
splCorObj = [];
splCorSmp = [0,0];
CorOL = [];
splCorObjBound = [];
splCorSmpBound = [0,0];

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
      
handles.femheadsScl = femheadsScl;
handles.femheadsCor = femheadsCor; 
handles.femheadsSag = femheadsSag; 
    
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

    axes(handles.axesCor); hold on; 
    visboundaries(CorOL); 
    plot(splCorSmp(:,2), splCorSmp(:,1), ':b');
    plot(splCorSmpBound(2), splCorSmpBound(1), '*b', 'LineWidth', 1.25);
    plot(femheadsCor(1,:), femheadsCor(2,:), 'ob', 'LineWidth', 1.25);
    
    axes(handles.axesSag); hold on; 
    visboundaries(SagOL); 
    plot(splSagSmp(:,2), splSagSmp(:,1), ':b');
    plot(splSagSmpBound(2), splSagSmpBound(1), '*b', 'LineWidth', 1.25);
    plot(femheadsSag(1,:), femheadsSag(2,:), 'ob', 'LineWidth', 1.25);
    
    if exist(fn2, 'file')
        disp(['loading ',fn2])
        load(fn2);
        
        set(handles.Writhe, 'String', ['Wr = ' num2str(Wr)]);
        
        splfiltPix = splfilt/handles.ifoCor.PixelSpacing(1);
                
        axes(handles.axesSag); hold on;
        plot(splfiltPix(:,1), splfiltPix(:,3), 'm', 'LineWidth', 2);
        axes(handles.axesCor); hold on;
        plot(splfiltPix(:,2), splfiltPix(:,3), 'm', 'LineWidth', 2);
    end
    
else
clear handles.splSclObj handles.splSclRng handles.splSclSmp 
clear handles.imgSagFilt handles.splSagObj handles.splSagSmp handles.SagOL
clear handles.imgCorFilt handles.splCorObj handles.splCorSmp handles.CorOL
clear handles.femheadsSag handles.femheadsCor handles.femheadsScl
end

% Update handles structure
guidata(hObject, handles);
