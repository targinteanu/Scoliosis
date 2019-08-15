function varargout = gui1(varargin)
% gui1 MATLAB code for gui1.fig
%      gui1, by itself, creates a new gui1 or raises the existing
%      singleton*.
%
%      H = gui1 returns the handle to a new gui1 or the handle to
%      the existing singleton*.
%
%      gui1('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in gui1.M with the given input arguments.
%
%      gui1('Property','Value',...) creates a new gui1 or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui1_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui1_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui1

% Last Modified by GUIDE v2.5 21-Oct-2018 18:58:39

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui1_OpeningFcn, ...
                   'gui_OutputFcn',  @gui1_OutputFcn, ...
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

% --- Executes just before gui1 is made visible.
function gui1_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui1 (see VARARGIN)

handles.cmap = lines; % colormap
handles.currentvertebra = 1; % (none)
handles.currentshape = 1;
vertebrae = cell(1, 26);

% change after testing 
fp = uigetdir; 
%fp = 'C:\Users\Toren\Desktop\scoliosis\patient4\4\C T L-Spines 18-12-29\C T L-SPINES - 3';
handles.filepath = fp;
[vol, imwidth, imheight, imdepth] = dicomreadvol(fp);
handles.imwidth = imwidth; handles.imheight = imheight; handles.imdepth = imdepth;

for n = 1:length(vertebrae)
    vertebrae{n} = false(size(vol));
end
handles.vertebrae = vertebrae;

currentz = floor(imdepth/2); 
imshow(vol(:, :, currentz));
% uncomment after testing 
%%{
filename = [fp, '\background_import.mat'];
if exist(filename, 'file')
    disp(['using ' filename])
    load(filename);
else
    [c, r] = getpts; c = floor(c); r = floor(r);
    marker2d = zeros(imheight, imwidth); 
    pts = (c-1)*imheight + r;
    marker2d(pts) = 1; 
    %imshow(marker2d);
    marker3d = zeros(size(vol)); marker3d(:,:,currentz) = marker2d;
    bg = imreconstruct(marker3d, vol); 
    save(filename, 'bg');
end
%}

% threshold: .85 for patient 4, .96 for patient 1
thresh = .96;

%load('test_background.mat')
vol = vol-bg; 
vol = histeq(vol); 
handles.vol = vol; [min(vol(:)) max(vol(:))]
volthreshed = vol > thresh; handles.volthreshed = volthreshed;

% code from seg(...)
CC = bwconncomp(handles.volthreshed, 6);
objs = CC.PixelIdxList;
minsize = 10000; 
szs = zeros(size(objs));
for n = 1:length(szs)
    szs(n) = length(objs{n});
end
objs = objs(szs > minsize); % objs{n} is an index matrix, not volume 
shapes = cell(size(objs)); Z = false(size(handles.volthreshed));
for n = 1:length(objs)
    I = Z; I(objs{n}) = true;
    shapes{n} = I; % shapes{n} is a volume matrix
end

% uncomment after testing 
% code from cleanshapes(...)
%%{
tf = zeros(size(shapes));
for n = 1:length(tf)
    %{
    cm = CenterOfMass3(shapes{n});
    k1 = min(imdepth, floor(cm(3)) + 5);
    imtoshow = .1*shapes{n}(:,:,k1);
    for k = 1:9
       imtoshow(shapes{n}(:,:,max(1,k1-k))) = .1 + k/10; 
    end
    imshow(imtoshow);
    %}
    imshow(show3dBW(shapes{n}, 10));
    x = inputdlg('accept (1) or reject (0)', 'accept/reject object', 1, {'0'});
    tf(n) = str2num(x{:});
end
tf = ~~tf; shapes = shapes(tf);
%}
handles.originalshapes = shapes;

filename = [fp, '\vertebrae_import.mat'];
if exist(filename, 'file')
    disp(['using ' filename])
    load(filename);
    for n = 1:length(Vertebrae)
        sz = size(Vertebrae{n}.Volume);
        if (length(sz) > 2)
            %vertebrae{n+1} = unpackage(Vertebrae{n}, false);
            vertebrae{n} = unpackage(Vertebrae{n}, false);
            for m = 1:length(shapes)
                %shapes{m}(vertebrae{n+1}) = false;
                shapes{m}(vertebrae{n}) = false;
            end
        end
    end
    handles.vertebrae = vertebrae;
end

handles.shapes = shapes;

handles.erosionelement = strel('sphere', 1);

% code from createdisplay(...)
dispvol = zeros(size(volthreshed));
for n = 1:length(shapes)
    dispvol = dispvol + n*shapes{n};
end
for m = 1:length(vertebrae)
    dispvol = dispvol + (m+n)*vertebrae{n};
end
dispvol = dispvol + 1;
dispvol = uint8(dispvol);
handles.dispvol = dispvol; % each binary shape is a unique integer
%handles.disprange = m+n+1;
handles.disprange = n+1;

handles.currentz = 1; 
showimg(hObject, eventdata, handles);

% Choose default command line output for gui1
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% This sets up the initial plot - only do when we are invisible
% so window can get raised using gui1.

% UIWAIT makes gui1 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui1_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in forwardbutton.
function forwardbutton_Callback(hObject, eventdata, handles)
% hObject    handle to forwardbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
z = handles.currentz; 

if z < handles.imdepth
    handles.currentz = z + 1;
    currentslizetxt.String = num2str(handles.currentz);
end

guidata(hObject, handles);
handles.currentshape = findcurrentshape(hObject, eventdata, handles);
guidata(hObject, handles);
showimg(hObject, eventdata, handles);

function backbutton_Callback(hObject, eventdata, handles)
% hObject    handle to backbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
z = handles.currentz; 

if z > 1
    handles.currentz = z - 1;
    currentslizetxt.String = num2str(handles.currentz);
end

guidata(hObject, handles);
handles.currentshape = findcurrentshape(hObject, eventdata, handles);
guidata(hObject, handles);
showimg(hObject, eventdata, handles);

% resetting shapes;
%{
handles.shapes = handles.originalshapes; 
for n = 1:length(handles.vertebrae)
    handles.shapes{1}(handles.vertebrae{n}) = false;
end
guidata(hObject, handles);
%}

% --------------------------------------------------------------------
function FileMenu_Callback(hObject, eventdata, handles)
% hObject    handle to FileMenu (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --------------------------------------------------------------------
function OpenMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to OpenMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
file = uigetfile('*.fig');
if ~isequal(file, 0)
    open(file);
end

% --------------------------------------------------------------------
function PrintMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to PrintMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
printdlg(handles.figure1)

% --------------------------------------------------------------------
function CloseMenuItem_Callback(hObject, eventdata, handles)
% hObject    handle to CloseMenuItem (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
selection = questdlg(['Close ' get(handles.figure1,'Name') '?'],...
                     ['Close ' get(handles.figure1,'Name') '...'],...
                     'Yes','No','Yes');
if strcmp(selection,'No')
    return;
end

delete(handles.figure1)


% --- Executes on selection change in popupmenu1.
function popupmenu1_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = get(hObject,'String') returns popupmenu1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu1
popup_sel_index = get(handles.popupmenu1, 'Value');

% change current vertebra
handles.currentvertebra = popup_sel_index;
guidata(hObject, handles);

% update display
[dispvol, disprange] = createdisplay(hObject, eventdata, handles);
handles.dispvol = dispvol; 
handles.disprange = disprange; 
guidata(hObject, handles);
showimg(hObject, eventdata, handles);
%dispvol = handles.dispvol;
%save('test', 'dispvol')

handles.currentshape = findcurrentshape(hObject, eventdata, handles);
guidata(hObject, handles);

function curshape = findcurrentshape(hObject, eventdata, handles)
x = zeros(size(handles.shapes));
for n = 1:length(x)
    x(n) = sum(sum(handles.shapes{n}(:,:,handles.currentz)));
end
[xmax, idx] = max(x);

if xmax
    curshape = idx;
else
    curshape = handles.currentshape;
end

% --- Executes during object creation, after setting all properties.
function popupmenu1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
     set(hObject,'BackgroundColor','white');
end

set(hObject, 'String', {'none', ...
    'C1', 'C2', 'C3', 'C4', 'C5', 'C6', 'C7', ...
    'T1', 'T2', 'T3', 'T4', 'T5', 'T6', 'T7', 'T8', 'T9', 'T10', 'T11', 'T12', ...
    'L1', 'L2', 'L3', 'L4', 'L5', ...
    'S1', 'S2', 'S3', 'S4', 'S5'});


function ptselectbtn_Callback(hObject, eventdata, handles)
[c, r] = getpts; c = floor(c); r = floor(r);
marker2d = false(handles.imheight, handles.imwidth); 
pts = (c-1)*handles.imheight + r;
marker2d(pts) = true; 
%imshow(marker2d);
slice = imreconstruct(marker2d, ...
    handles.shapes{handles.currentshape}(:,:,handles.currentz));
    % previously handles.volthreshed(:,:,handles.currentz)
handles.vertebrae{handles.currentvertebra}(:,:,handles.currentz) = ...
    handles.vertebrae{handles.currentvertebra}(:,:,handles.currentz) | slice;

handles.shapes{handles.currentshape}(handles.vertebrae{handles.currentvertebra}) = false;

guidata(hObject, handles);
showimg(hObject, eventdata, handles);

function roiselectbtn_Callback(hObject, eventdata, handles)
roi = roipoly; 
slice = roi & ...
    handles.shapes{handles.currentshape}(:,:,handles.currentz);
    % previously handles.volthreshed(:,:,handles.currentz);
handles.vertebrae{handles.currentvertebra}(:,:,handles.currentz) = ...
    handles.vertebrae{handles.currentvertebra}(:,:,handles.currentz) | slice;

handles.shapes{handles.currentshape}(handles.vertebrae{handles.currentvertebra}) = false;

guidata(hObject, handles);
showimg(hObject, eventdata, handles);

function selectallbtn_Callback(hObject, eventdata, handles)
handles.vertebrae{handles.currentvertebra} = ...
    handles.vertebrae{handles.currentvertebra} | ...
    handles.shapes{handles.currentshape};

handles.shapes{handles.currentshape}(handles.vertebrae{handles.currentvertebra}) = false;

guidata(hObject, handles);
showimg(hObject, eventdata, handles);

function clearbutton_Callback(hObject, eventdata, handles)
slice = handles.vertebrae{handles.currentvertebra}(:,:,handles.currentz);
handles.shapes{handles.currentshape}(:,:,handles.currentz) = ...
    handles.shapes{handles.currentshape}(:,:,handles.currentz) | slice;

slice = false(handles.imheight, handles.imwidth);
handles.vertebrae{handles.currentvertebra}(:,:,handles.currentz) = slice;

guidata(hObject, handles);
showimg(hObject, eventdata, handles);

function savebutton_Callback(hObject, eventdata, handles)
filename = inputdlg('save as:', 'save', 1, {'vertebrae'});
filename = [handles.filepath, '\', filename{:}];
%Vertebrae = handles.vertebrae(2:end); 
Vertebrae = handles.vertebrae;
for n = 1:length(Vertebrae)
    Vertebrae{n} = package(Vertebrae{n});
end
save(filename, 'Vertebrae');

function segbutton_Callback(hObject, eventdata, handles)
%seg(hObject, eventdata, handles);

%function seg(hObject, eventdata, handles)
newvol = false(size(handles.volthreshed));
for n = 1:length(handles.shapes)
    newvol = newvol | handles.shapes{n};
end

CC = bwconncomp(newvol, 6);
objs = CC.PixelIdxList;
minsize = 1000; 
szs = zeros(size(objs));
for n = 1:length(szs)
    szs(n) = length(objs{n});
end
objs = objs(szs > minsize); % objs{n} is an index matrix, not volume 
shapes = cell(size(objs)); Z = false(size(handles.volthreshed));
for n = 1:length(objs)
    I = Z; I(objs{n}) = true;
    shapes{n} = I; % shapes{n} is a volume matrix
end
handles.shapes = shapes;
guidata(hObject, handles);
[dispvol, disprange] = createdisplay(hObject, eventdata, handles);
handles.dispvol = dispvol; 
handles.disprange = disprange; 
guidata(hObject, handles);
showimg(hObject, eventdata, handles);

function cleanshapes(hObject, eventdata, handles)
shapes = handles.shapes;
tf = zeros(size(shapes));
for n = 1:length(tf)
    cm = CenterOfMass3(shapes{n});
    k1 = min(imdepth, floor(cm(3)) + 5);
    imtoshow = .1*shapes{n}(:,:,k1);
    for k = 1:9
       imtoshow(shapes{n}(:,:,max(1,k1-k))) = .1 + k/10; 
    end
    imshow(imtoshow);
    x = inputdlg('accept (1) or reject (0)', 'accept/reject object', 1, {'0'});
    tf(n) = str2num(x{:});
end
tf = ~~tf; shapes = shapes(tf);
handles.shapes = shapes; 
guidata(hObject, handles);
createdisplay(hObject, eventdata, handles);

function [dispvol, disprange] = createdisplay(hObject, eventdata, handles)
shapes = handles.shapes;
vertebrae = handles.vertebrae;
dispvol = zeros(size(handles.volthreshed));
for n = 1:length(shapes)
    dispvol = dispvol + n*shapes{n};
end
%for m = 1:length(vertebrae)
for m = 2:length(vertebrae)
    dispvol = dispvol + (mod(m,2) + n + 1)*vertebrae{m};
end
dispvol = dispvol + 1;
dispvol = uint8(dispvol);
%handles.dispvol = dispvol; % each binary shape is a unique integer
%handles.disprange = m+n+1;
%disprange = m+n+1;
disprange = n+3;
%handles.disprange
%guidata(hObject, handles);
%showimg(hObject, eventdata, handles);

function showimg(hObject, eventdata, handles)
GB = handles.dispvol(:,:,handles.currentz);
%handles.disprange
%dispvol = handles.dispvol;
%save('test', 'GB', 'dispvol')
imtoshow = zeros([size(GB), 3]);
imtoshow(:,:,1) = GB; imtoshow(:,:,2) = GB; imtoshow(:,:,3) = GB;
R = handles.vertebrae{handles.currentvertebra}(:,:,handles.currentz);
%imtoshow(:,:,1) = imtoshow(:,:,1) - R;

GB(R) = 0; imtoshow(:,:,1) = GB;

imtoshow = imtoshow/handles.disprange;
imshow(imtoshow);
%colormap(handles.cmap);
