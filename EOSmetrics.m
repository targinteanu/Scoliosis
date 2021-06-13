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

% Last Modified by GUIDE v2.5 08-Mar-2021 03:31:51

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

handles.defaultFilter = {'.0033', '.0034', '1', '50'};

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
    'Lewiner Torsion', ...
    'Plumbline X', ...
    'Plumbline Y', ...
    'Plumbline 3D', ...
    'Writhe', ...
    'Cobb Angle', ...
    'Sag/Cor Vertical Alignment', ...
    'Re-Filter'};

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
    @(HO,ED,H) showLewinerQuantity(HO,ED,H,4), ...
    @(HO,ED,H) showPlumblineDistance(HO,ED,H,1), ...
    @(HO,ED,H) showPlumblineDistance(HO,ED,H,2), ...
    @(HO,ED,H) showPlumblineDistance(HO,ED,H,3), ...
    @(HO,ED,H) showWrithe(HO,ED,H), ...
    @(HO,ED,H) cobbAngleMinMax(HO,ED,H), ...
    @(HO,ED,H) SCVA(HO,ED,H), ...
    @(HO,ED,H) refilter(HO,ED,H)};

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


% --- Executes on button press in localMinMax.
function localMinMax_Callback(hObject, eventdata, handles)
% hObject    handle to localMinMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

ifoCor = handles.ifoCor; ifoSag = handles.ifoSag;
SplSclHt = handles.SplSclHt;
x = SplSclHt(:,1); y = SplSclHt(:,2); z = SplSclHt(:,3); h = SplSclHt(:,4);

[pk_min,lc_min,pw_min,pp_min] = findpeaks(-h);%, 'MinPeakWidth', length(h)/20); 
[pk_max,lc_max,pw_max,pp_max] = findpeaks(h);%, 'MinPeakWidth', length(h)/20);
figure; plot(h); hold on; grid on;
errorbar(lc_min, -pk_min, pp_min,pp_min, pw_min,pw_min, 'ob'); 
errorbar(lc_max, pk_max, pp_max,pp_max, pw_max,pw_max, 'or');
pp_min = 100*pp_min/max(pp_min); pp_max = 100*pp_max/max(pp_max);

handles.lc_min = lc_min; handles.pw_min = pw_min; handles.pp_min = pp_min; handles.pk_min = pk_min;
handles.lc_max = lc_max; handles.pw_max = pw_max; handles.pp_max = pp_max; handles.pk_max = pk_max;

axes(handles.axes3); hold on;
plot3(x(lc_min), -y(lc_min), -z(lc_min), 'ob', 'LineWidth', 1.25);
plot3(x(lc_max), -y(lc_max), -z(lc_max), 'or', 'LineWidth', 1.25);

axes(handles.axesSag); hold on; % x
errorbar(x(lc_min)/ifoSag.PixelSpacing(2), z(lc_min)/ifoSag.PixelSpacing(1), ...
    pw_min,pw_min, pp_min,pp_min, 'ob', 'LineWidth', 1.25);
errorbar(x(lc_max)/ifoSag.PixelSpacing(2), z(lc_max)/ifoSag.PixelSpacing(1), ...
    pw_max,pw_max, pp_max,pp_max, 'or', 'LineWidth', 1.25);

axes(handles.axesCor); hold on; % y
errorbar(y(lc_min)/ifoCor.PixelSpacing(2), z(lc_min)/ifoCor.PixelSpacing(1), ...
    pw_min,pw_min, pp_min,pp_min, 'ob', 'LineWidth', 1.25);
errorbar(y(lc_max)/ifoCor.PixelSpacing(2), z(lc_max)/ifoCor.PixelSpacing(1), ...
    pw_max,pw_max, pp_max,pp_max, 'or', 'LineWidth', 1.25);

guidata(hObject, handles);


function SCVA(hObject, eventdata, handles)

ifoCor = handles.ifoCor; ifoSag = handles.ifoSag;

XYZ = handles.splfilt; 
p1 = XYZ(1,:); p2 = XYZ(end,:); 
p3 = [p1(1), p1(2), p2(3)];
p4 = .5*(p2-p3) + p3;
p123 = [p1; p3; p2];

% x - sag
SVA = -(p1(1)-p2(1));
axes(handles.axesSag);
xy = p123(:,[1,3])./flipud(ifoSag.PixelSpacing)';
P4 = p4([1,3])./flipud(ifoSag.PixelSpacing)';
plot(xy(:,1), xy(:,2), 'g');
text(P4(1), P4(2), ...
    [num2str(SVA) ' mm'], ...
    'VerticalAlignment', 'top', 'Color', 'g');

% y - cor
CVA = (p1(2)-p2(2));
axes(handles.axesCor);
xy = p123(:,[2,3])./flipud(ifoCor.PixelSpacing)';
P4 = p4([2,3])./flipud(ifoCor.PixelSpacing)';
plot(xy(:,1), xy(:,2), 'g');
text(P4(1), P4(2), ...
    [num2str(CVA) ' mm'], ...
    'VerticalAlignment', 'top', 'Color', 'g');

% 3D
VA3D = norm(p3-p2);
axes(handles.axes3);
plot3(p123(:,1), -p123(:,2), -p123(:,3), 'g');
text(p4(1), -p4(2), -p4(3), ...
    [num2str(VA3D), ' mm'], ...
    'VerticalAlignment', 'top', 'Color', 'g');


function [theta, n1, n2] = numericCobbAngle(R, t1, t2)
dR = diff(R); 
Ri = .5 * ( R(2:end,:) + R(1:(end-1),:) ); dRi = diff(Ri);
ds = sum(dR.^2, 2).^.5;
dsi = sum(dRi.^2, 2).^.5;
dR = dR./ds;
ddR = diff(dR);
ddR = ddR./dsi;
xy1 = R(t1,:); xy2 = R(t2,:);
t1 = t1-1; t2 = t2-1;
t1 = min(t1, size(ddR,1)); t1 = max(t1,1);
t2 = min(t2, size(ddR,1)); t2 = max(t2,1);
n1 = ddR(t1,:); n1 = n1/norm(n1);
n2 = ddR(t2,:); n2 = n2/norm(n2);
if length(xy1) < 3
    [~,sgn] = intersect2d(xy1', n1', xy2', n2'); sgn=sign(sgn); 
    n1=sgn(1)*n1; n2=sgn(2)*n2;
end
theta = acos(n1 * n2') * 180/pi;
%theta = min(theta, 180-theta); % come up with a better way of picking the right angle 

function [xy3,c] = intersect2d(xy1, a, xy2, b)
c = ([a, b])^-1 * (xy1 - xy2); c(1) = -c(1);
xy3 = c(2)*b + xy2;

function dispCobbAngle(hObject, eventdata, handles, R, t1, t2)

ifoCor = handles.ifoCor; ifoSag = handles.ifoSag;

% x - sag
xy = R(:,[1,3]); 
[thta, a, b] = numericCobbAngle(xy, t1, t2);
xy1 = xy(t1,:); xy2 = xy(t2,:); 
%xy3 = intersect2d(xy1', a', xy2', b');
xy3 = xy1 + 100*a; xy4 = xy2 + 100*b;
if xy3(1) > xy1(1)
    al = 'left';
else
    al = 'right';
end
xytxt = xy3./flipud(ifoSag.PixelSpacing)'; 
xytxt(1) = max(xytxt(1), 1); xytxt(1) = min(xytxt(1), ifoSag.Width);
axes(handles.axesSag); hold on; 
plot([xy1(1), xy3(1), xy4(1), xy2(1)]/ifoSag.PixelSpacing(2), ...
    [xy1(2), xy3(2), xy4(2), xy2(2)]/ifoSag.PixelSpacing(1), 'g');
text(xytxt(1), xytxt(2), ...
    [num2str(thta) '\circ'], ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', al, 'Color', 'g');

% y - cor
xy = R(:,[2,3]); 
[thta, a, b] = numericCobbAngle(xy, t1, t2);
xy1 = xy(t1,:); xy2 = xy(t2,:); 
%xy3 = intersect2d(xy1', a', xy2', b');
xy3 = xy1 + 100*a; xy4 = xy2 + 100*b;
if xy3(1) > xy1(1)
    al = 'left';
else
    al = 'right';
end
xytxt = xy3./flipud(ifoCor.PixelSpacing)'; 
xytxt(1) = max(xytxt(1), 1); xytxt(1) = min(xytxt(1), ifoCor.Width);
axes(handles.axesCor); hold on; 
plot([xy1(1), xy3(1), xy4(1), xy2(1)]/ifoCor.PixelSpacing(2), ...
    [xy1(2), xy3(2), xy4(2), xy2(2)]/ifoCor.PixelSpacing(1), 'g');
text(xytxt(1), xytxt(2), ...
    [num2str(thta) '\circ'], ...
    'VerticalAlignment', 'middle', 'HorizontalAlignment', al, 'Color', 'g');

% 3d 
[theta, n1, n2] = numericCobbAngle(R, t1, t2);
xyz1 = R(t1,:); xyz2 = R(t2,:); 
xyz1 = xyz1 + [zeros(size(n1)); n1*100]; xyz2 = xyz2 + [zeros(size(n2)); n2*100];
xyz3 = .5*( xyz1(2,:) + xyz2(2,:) );
axes(handles.axes3); hold on;
plot3(xyz1(:,1), -xyz1(:,2), -xyz1(:,3), 'g');
plot3(xyz2(:,1), -xyz2(:,2), -xyz2(:,3), 'g');
text(xyz3(:,1), -xyz3(:,2), -xyz3(:,3), [num2str(theta) '\circ'], 'Color', 'g');

function idx = minSgn(vals, sgn)
origIdx = 1:length(vals);
vals = vals*sgn; 
sel = vals > 0;
origIdx = origIdx(sel); vals = vals(sel);
[~,i] = min(vals); idx = origIdx(i);

function cobbAngleMinMax(hObject, eventdata, handles)
SplSclHt = handles.SplSclHt;
x = SplSclHt(:,1); y = SplSclHt(:,2); z = SplSclHt(:,3); R = [x,y,z];
neutIdx = [1; handles.lc_min; length(z)]; apIdx = handles.lc_max;
boundIdx = zeros(length(apIdx),2);
for i = 1:length(apIdx)
    d = z(neutIdx)-z(apIdx(i));
    boundIdx(i,1) = neutIdx(minSgn(d, -1));
    boundIdx(i,2) = neutIdx(minSgn(d, 1));
end
dispCobbAngles(hObject, eventdata, handles, R, boundIdx);

function dispCobbAngles(hObject, eventdata, handles, R, boundIdx)
for i = 1:size(boundIdx,1)
    t1 = boundIdx(i,1); t2 = boundIdx(i,2);
    dispCobbAngle(hObject, eventdata, handles, R, t1, t2);
end


% --- Executes on button press in saveButton.
function saveButton_Callback(hObject, eventdata, handles)
% hObject    handle to saveButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
handles.splfilt = handles.splfiltNew;
handles.splfiltNew = [];
guidata(hObject, handles);
showpatient(hObject, eventdata, handles);


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
x = XYZ(:,1); y = XYZ(:,2); z = XYZ(:,3); h = heatvals;
cla;
surface([x,x]', -[y,y]', -[z,z]', [h,h]', ...
    'facecol', 'no', 'edgecol', 'interp', 'LineWidth', 5);
grid on;
colorbar;
handles.SplSclHt = [XYZ,h];
guidata(hObject, handles);
    hold on; 
    fhXYZ = handles.femheadsScl;
    plot3(fhXYZ(1,:), -fhXYZ(2,:), -fhXYZ(3,:), 'or');


function showPlumblineDistance(hObject, eventdata, handles, varToShow)
projuv = @(u,v) ((u*v')/(v*v'))*v;
R = handles.splfilt;
ifoSag = handles.ifoSag; ifoCor = handles.ifoCor;

if (nargin < 4)|(varToShow==3)
    R0 = R;
else
    R0 = R(:,[varToShow,3]);
end

plumb = R0([1,end],:); plumb3 = R([1,end],:);
vPlumb = diff(plumb); 
R0 = R0 - plumb(1,:);

axes(handles.axes3); 
hold on; plot3(plumb3(:,1), -plumb3(:,2), -plumb3(:,3), 'c');
axes(handles.axesSag); 
hold on; plot(plumb3(:,1)/ifoSag.PixelSpacing(2), ...
    plumb3(:,3)/ifoSag.PixelSpacing(1), 'c'); 
axes(handles.axesCor); 
hold on; plot(plumb3(:,2)/ifoCor.PixelSpacing(2), ...
    plumb3(:,3)/ifoCor.PixelSpacing(1), 'c'); 

dist = arrayfun(@(i) norm( R0(i,:) - projuv(R0(i,:), vPlumb) ), 1:size(R0,1))';
showLinearHeatmap(hObject, eventdata, handles, R, dist);


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
R = handles.splfilt; 
dR = diff(R); ds = sum(dR.^2, 2).^.5; dR = dR./ds;
Ri = .5 * ( R(2:end,:) + R(1:(end-1),:) );
dRi = diff(Ri); dsi = sum(dRi.^2, 2).^.5; dRi = dRi./dsi;
ddR = diff(dR)./dsi;
K = arrayfun(@(i) norm(cross(dRi(i,:), ddR(i,:)))./( norm(dRi(i,:)).^3 ), ...
    1:size(dRi,1));
showLinearHeatmap(hObject, eventdata, handles, R(2:(end-1),:), K');

function showNumericTorsion(hObject, eventdata, handles)
% change to include ds
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


function showWrithe(hObject, eventdata, handles)
cm = handles.splfilt;
range = 2:size(cm,1);
Wr = zeros(length(range));
for i = range
    for j = range
        dr1 = cm(i,:) - cm((i-1),:);
        dr2 = cm(j,:) - cm((j-1),:);
        dr12 = cm(i,:) - cm(j,:);
        if norm(dr12)
            Wr(i,j) = ( cross(dr1, dr2)*dr12' )/( norm(dr12)^3 );
        end
    end
end
Wr = Wr/(4*pi); Writhe = sum(Wr(:));
figure; imagesc(Wr); title(['Writhe = ' num2str(Writhe)]); colorbar;


function refilter(hObject, eventdata, handles)
filtfilt3d = @(filtobj, sig3d) [filtfilt(filtobj, sig3d(:,1)), filtfilt(filtobj, sig3d(:,2)), sig3d(:,3)];
sig = handles.splSclSmp; splSclBnd = handles.splSclBnd;
fs = 1/mean(diff(sig(:,3)));
showpatient(hObject, eventdata, handles);

if ~isempty(handles.splfiltNew)
    sigf = handles.splfiltNew;
    axes(handles.axes3); hold on; 
    plot3(sigf(:,1), -sigf(:,2), -sigf(:,3), 'c', 'LineWidth', 1.5);
end

fTypeOpts = {'cheby2', 'cheby1', 'ellip', 'butter'};
[sel, ok] = listdlg('ListString', fTypeOpts, 'SelectionMode', 'single');
if ok
    
fType = fTypeOpts{sel};

promptFields  = {'Passband Frequency (1/mm)', ...
    'Stopband Frequency (1/mm)', ...
    'Passband Ripple (dB)', ...
    'Stopband Attenuation (dB)'};
filtspecs = inputdlg(promptFields, 'Build New Filter', 1, handles.defaultFilter);
handles.defaultFilter = filtspecs;
fPass = str2num(filtspecs{1}); fStop = str2num(filtspecs{2}); 
aPass = str2num(filtspecs{3}); aStop = str2num(filtspecs{4}); 

filtobj = designfilt('lowpassiir', 'SampleRate',fs, ...
    'DesignMethod',fType, ...
    'PassbandFrequency',fPass, 'StopbandFrequency',fStop, ...
    'PassbandRipple',aPass, 'StopbandAttenuation',aStop);

plotFilterAndSignal(filtobj, sig(:,1:2), fs);
showpatient(hObject, eventdata, handles);
sigf = filtfilt3d(filtobj, sig); sigf = sigf(splSclBnd(2):splSclBnd(1),:);
handles.splfiltNew = sigf;
axes(handles.axes3); hold on; 
plot3(sigf(:,1), -sigf(:,2), -sigf(:,3), 'c', 'LineWidth', 1.5);

guidata(hObject, handles);
end


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
        
        axes(handles.axes3); cla; colorbar off;
        plot3(splSclSmp(:,1), -splSclSmp(:,2), -splSclSmp(:,3), ':b'); hold on;
        plot3(splfilt(:,1), -splfilt(:,2), -splfilt(:,3), 'm', 'LineWidth', 2);
        grid on; update3Dview(eventdata, handles); 
        xlim([0, size(imgSag,2)]*ifoSag.PixelSpacing(2));
        ylim(-[size(imgCor,2), 0]*ifoCor.PixelSpacing(2));
        zlim(-[size(imgCor,1), 0]*ifoCor.PixelSpacing(1));
        
            hold on;
            fhXYZ = handles.femheadsScl;
            plot3(fhXYZ(1,:), -fhXYZ(2,:), -fhXYZ(3,:), 'or');
            
            femheadsSag = handles.femheadsSag; femheadsCor = handles.femheadsCor;
            axes(handles.axesSag); 
            plot(femheadsSag(1,:), femheadsSag(2,:), 'ob', 'LineWidth', 1.25);
            axes(handles.axesCor); 
            plot(femheadsCor(1,:), femheadsCor(2,:), 'ob', 'LineWidth', 1.25);

function shownewpatient(hObject, eventdata, handles)
set(handles.PatientNum, 'String', ...
    num2str(handles.patient_list(handles.current_patient)));
handles.metricSelecter.Value = 1;

% initialize vars that may not have been set
femheadsScl = [0,0;0,0;0,0];
femheadsCor = [0,0;0,0];
femheadsSag = [0,0;0,0];
handles.splfiltNew = [];

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

    handles.femheadsScl = femheadsScl; 
    handles.femheadsSag = femheadsSag; 
    handles.femheadsCor = femheadsCor;

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
view(-handles.slider1.Value*90, handles.slider2.Value*90);


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


function plotFilterAndSignal(filtObj, x, fs)
X = abs(fftshift(fft(x)));
w = linspace(-fs/2, fs/2, size(X,1));

wmax = w(end-1); 
wordr = log10(wmax);
wunitlog = 3*floor(wordr/3); wunit = 10^wunitlog;

w = w/wunit;
pwr = X.^2; pwr_dB = 10*log10(pwr./sum(pwr));

fvtool(filtObj); hold on;
plot(w, pwr_dB, 'g');
xlabel(['Frequency, ' num2str(wunit) '/mm']);
