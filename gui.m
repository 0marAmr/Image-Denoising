function varargout = gui(varargin)
% GUI MATLAB code for gui.fig
%      GUI, by itself, creates a new GUI or raises the existing
%      singleton*.
%
%      H = GUI returns the handle to a new GUI or the handle to
%      the existing singleton*.
%
%      GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in GUI.M with the given input arguments.
%
%      GUI('Property','Value',...) creates a new GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help gui

% Last Modified by GUIDE v2.5 18-Dec-2021 22:04:56

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @gui_OpeningFcn, ...
                   'gui_OutputFcn',  @gui_OutputFcn, ...
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


% --- Executes just before gui is made visible.
function gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to gui (see VARARGIN)

% Choose default command line output for gui
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes during object creation, after setting all properties.
function after_CreateFcn(hObject, eventdata, handles)
function dt_Callback(hObject, eventdata, handles)
function dt_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
function iteration_Callback(hObject, eventdata, handles)
function iteration_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in apply.
function apply_Callback(hObject, eventdata, handles)


dt = get(handles.dt,'string');
dt = str2double(dt);
iteration = get(handles.iteration,'string');
iteration = str2double(iteration);

%apply the denosing
choice = get(handles.popup3,'value');
flag = get(handles.speckle,'Value');
switch choice
   
    case 1
	J = getappdata(0,'J');
        [frth]=fpdepyou(J,iteration,dt);
        [k] = speckle_sol(frth);
        

        if flag
            axes(handles.after);
            imshow(k);
        else     
            axes(handles.after);
            imshow(frth);
        end
        
    case 2
        J = getappdata(0,'J');
        [frth]=fpdepyou(J,iteration,dt);
        [k] = speckle_sol(frth);
        if flag
            row = get(handles.row_num,'string');
            row = str2double(row);
            row = uint8(row);
            axes(handles.after);
            [~,x] = size(k);
            x1 = linspace(0,x,x);
            y1 = k(row,:);
            plot(x1,y1);
        else
            row = get(handles.row_num,'string');
            row = str2double(row);
            row = uint8(row);
            axes(handles.after);
            [~,x] = size(frth);
            x1 = linspace(0,x,x);
            y1 = frth(row,:);
            plot(x1,y1);
        end
        
end








function popup_noise_Callback(hObject, eventdata, handles)  
function popup_noise_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in nosie.
function nosie_Callback(hObject, eventdata, handles)
cla

%choose photo
choice = get(handles.photo,'value');
switch choice
   
    case 1
        C = imread('photo1.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
    case 2
        C = imread('photo2.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
    case 3
        C = imread('photo3.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
    case 4
        C = imread('photo4.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
    case 5
        C = imread('photo5.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
end



%nosie on before
choice = get(handles.popup_noise,'value');

switch choice
   
    case 1
	C = getappdata(0,'C');
    J = imnoise(C,'Gaussian');  
	setappdata(0,'J',J);
    case 2
	C = getappdata(0,'C');
    J = imnoise(C,'salt & pepper');
	setappdata(0,'J',J);
    case 3        
	C = getappdata(0,'C');
    J = imnoise(C,'poisson');
	setappdata(0,'J',J);
    case 4        
	C = getappdata(0,'C');
    J = imnoise(C,'speckle');
	setappdata(0,'J',J);
end

% (photo or row) original and before
choice = get(handles.popup3,'value');

switch choice
   
    case 1
        axes(handles.original);
        imshow(C);
        axes(handles.before);
        imshow(J); 

    case 2
        row = get(handles.row_num,'string');
        row = str2double(row);
        row = uint8(row);
        axes(handles.original);
        [~,x] = size(C);
        x1 = linspace(0,x,x);
        y1 = C(row,:);
        plot(x1,y1);
        axes(handles.before);
        [~,x] = size(C);
        x2 = linspace(0,x,x);
        y2 = J(row,:);
        plot(x2,y2);
        
end

set(handles.apply,'Enable','on');


function popup3_Callback(hObject, eventdata, handles)
cla
cla
cla
%choose photo
choice = get(handles.photo,'value');
switch choice
   
    case 1
        C = imread('photo1.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
    case 2
        C = imread('photo2.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
    case 3
        C = imread('photo3.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
    case 4
        C = imread('photo4.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
    case 5
        C = imread('photo5.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
end


%(photo or row) original
choice = get(handles.popup3,'value');
switch choice
   
    case 1
        axes(handles.original);
        imshow(C);

    case 2
        row = get(handles.row_num,'string');
        row = str2double(row);
        row = uint8(row);
        axes(handles.original);
        [~,x] = size(C);
        x1 = linspace(0,x,x);
        y1 = C(row,:);
        plot(x1,y1);
        
end


function popup3_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function row_num_Callback(hObject, eventdata, handles)
function row_num_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in photo.
function photo_Callback(hObject, eventdata, handles)
cla
cla

choice = get(handles.photo,'value');
switch choice
   
    case 1
        C = imread('photo1.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
    case 2
        C = imread('photo2.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
    case 3
        C = imread('photo3.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
    case 4
        C = imread('photo4.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
    case 5
        C = imread('photo5.jpg');
        C = rgb2gray(C);
        setappdata(0,'C',C);
end


choice = get(handles.popup3,'value');
switch choice
   
    case 1
        axes(handles.original);
        imshow(C);

    case 2
        row = get(handles.row_num,'string');
        row = str2double(row);
        row = uint8(row);
        axes(handles.original);
        [~,x] = size(C);
        x1 = linspace(0,x,x);
        y1 = C(row,:);
        plot(x1,y1);
        
end


% hObject    handle to photo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns photo contents as cell array
%        contents{get(hObject,'Value')} returns selected item from photo


% --- Executes during object creation, after setting all properties.
function photo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to photo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in speckle.
function speckle_Callback(hObject, eventdata, handles)
% hObject    handle to speckle (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of speckle
