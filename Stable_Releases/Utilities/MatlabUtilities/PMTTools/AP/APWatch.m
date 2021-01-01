function varargout = APWatch(varargin)
% APGLOBAL M-file for apglobal.fig
%      APGLOBAL, by itself, creates a new APGLOBAL or raises the existing
%      singleton*.
%
%      H = APGLOBAL returns the handle to a new APGLOBAL or the handle to
%      the existing singleton*.
%
%      APGLOBAL('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in APGLOBAL.M with the given input arguments.
%
%      APGLOBAL('Property','Value',...) creates a new APGLOBAL or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before APWatch_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to APWatch_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help apglobal

% Last Modified by GUIDE v2.5 19-Feb-2010 15:58:10

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @APWatch_OpeningFcn, ...
                   'gui_OutputFcn',  @APWatch_OutputFcn, ...
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



% --- Executes just before apglobal is made visible.
function APWatch_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to apglobal (see VARARGIN)

% Choose default command line output for apglobal
handles.output = hObject;

handles.serials = getPMTserials;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes apglobal wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = APWatch_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


function aphistory_Callback(hObject, eventdata, handles)

% Clear axes
axes(handles.mainviewer);
cla(handles.mainviewer,'reset');

try
    APHist(handles.PMTSerial)
    set(handles.status,'String','')
catch
    set(handles.status,'String','No Datasets Found')
end

guidata(hObject, handles);

function pmtselect_Callback(hObject, eventdata, handles)

handles.PMTSelection = get(hObject,'Value');
PMTAllSerials = get(hObject,'String');
handles.PMTSerial = char(PMTAllSerials{handles.PMTSelection});

guidata(hObject, handles);

function pmtselect_CreateFcn(hObject, eventdata, handles)

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

serials = getPMTserials;

set(hObject,'String',serials);

guidata(hObject, handles);


function apglobal_Callback(hObject, eventdata, handles)

% Clear axes
axes(handles.mainviewer);
cla(handles.mainviewer,'reset');

set(handles.status,'String','Loading...'); drawnow
serials = getPMTserials;

APGlobalStatus(serials);

set(handles.status,'String','')

