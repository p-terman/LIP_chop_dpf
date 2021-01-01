function varargout = PMTGainCalibration(varargin)
% PMTGAINCALIBRATION M-file for PMTGainCalibration.fig
%      PMTGAINCALIBRATION, by itself, creates a new PMTGAINCALIBRATION or raises the existing
%      singleton*.
%
%      H = PMTGAINCALIBRATION returns the handle to a new PMTGAINCALIBRATION or the handle to
%      the existing singleton*.
%
%      PMTGAINCALIBRATION('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in PMTGAINCALIBRATION.M with the given input arguments.
%
%      PMTGAINCALIBRATION('Property','Value',...) creates a new PMTGAINCALIBRATION or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before PMTGainCalibration_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to PMTGainCalibration_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help PMTGainCalibration

% Last Modified by GUIDE v2.5 06-Oct-2009 16:28:29

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @PMTGainCalibration_OpeningFcn, ...
                   'gui_OutputFcn',  @PMTGainCalibration_OutputFcn, ...
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

%% Initialization
% --- Executes just before PMTGainCalibration is made visible.
function PMTGainCalibration_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to PMTGainCalibration (see VARARGIN)

% Choose default command line output for PMTGainCalibration
handles.output = hObject;

% set(hObject,'toolbar','figure');
set(handles.LoadingStatus,'String','')
set(handles.filepath,'String','/Volumes/paRAID02/analysis/matlab/LUX01/PMT_RQs')
handles.pathbase = get(handles.filepath,'String');
handles.filename = get(handles.filenamefield,'String');
handles.leftedge = get(handles.LeftEdgeValue,'String');
handles.rightedge = get(handles.RightEdgeValue,'String');
handles.currentCh = 1;

% Get this from the LUG
handles.maxCh = 4;

% Update handles structure
guidata(hObject, handles);


%% Main function
function varargout = PMTGainCalibration_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
guidata(hObject, handles); %updates the handles


%% Loading
function filepath_Callback(hObject, eventdata, handles)
handles.pathbase = get(hObject,'String');
guidata(hObject, handles); %updates the handles

function filenamefield_Callback(hObject, eventdata, handles)
handles.filename = get(hObject,'String');
guidata(hObject, handles); %updates the handles

function AnalyzeDataButton_Callback(hObject, eventdata, handles)

set(handles.LoadingStatus,'String','Loading...'); drawnow;
handles.calibration = LUXLoadPMTCal(handles.filename,handles.pathbase);
set(handles.LoadingStatus,'String','Done Loading'); drawnow;

[xx nn] = PlotSpheSpectrum(hObject,handles,handles.currentCh);
FitSpheSpectrum(hObject,handles,xx,nn);

guidata(hObject, handles); %updates the handles

%% Ch Navigation
function next_Callback(hObject, eventdata, handles)
if handles.currentCh < handles.maxCh
    handles.currentCh = handles.currentCh + 1;
    [xx nn] = PlotSpheSpectrum(hObject,handles,handles.currentCh);
    FitSpheSpectrum(hObject,handles,xx,nn)
    set(handles.ChTextEdit,'String',handles.currentCh); drawnow;
    guidata(hObject, handles); %updates the handles
end

function previous_Callback(hObject, eventdata, handles)
if handles.currentCh > 1
    handles.currentCh = handles.currentCh - 1;
    [xx nn] = PlotSpheSpectrum(hObject,handles,handles.currentCh);
    FitSpheSpectrum(hObject,handles,xx,nn)
    set(handles.ChTextEdit,'String',handles.currentCh); drawnow;
    guidata(hObject, handles); %updates the handles
end

function ChTextEdit_Callback(hObject, eventdata, handles)
handles.ChTemp = str2double(get(hObject,'String'));
guidata(hObject, handles); %updates the handles

function GoToCh_Callback(hObject, eventdata, handles)
if handles.ChTemp >= 1 && handles.ChTemp <= handles.maxCh
    handles.currentCh = handles.ChTemp;
    [xx nn] = PlotSpheSpectrum(hObject,handles,handles.currentCh);
    FitSpheSpectrum(hObject,handles,xx,nn);
    set(handles.ChTextEdit,'String',handles.currentCh); drawnow;
    guidata(hObject, handles); %updates the handles
end

%% Plot and Fit Sphe Spectrum
function [xx nn] = PlotSpheSpectrum(hObject,handles,Ch)
axes(handles.SpheSpectrumPlot);
cla(handles.SpheSpectrumPlot,'reset')

areas = handles.calibration.areas_mVns(Ch,:);
areas_sort = sort(areas);
[nn xx] = hist(areas_sort(1:end-10),400);
plotlogstairs(xx,nn);
axis([-20 max(xx)*1.5 1 max(nn)*3])
xlabel('Pulse area [mVns]')
title(sprintf('Sphe Spectrum for Ch. %d',Ch))

function FitSpheSpectrum(hObject,handles,xx,nn)

% Fit noise
axes(handles.SpheSpectrumPlot); hold on

binwidth = xx(2)-xx(1);
noise_R = binwidth*4;

leg{1} = 'Raw Data';
[h_noise w_noise] = max(nn);
fit_noise_range = inrange(xx,xx(w_noise)-noise_R,xx(w_noise)+noise_R);
[A_noise Y_noise] = FH_Gauss(xx(fit_noise_range),nn(fit_noise_range));
xx_smooth = linspace(xx(1),xx(end),500);
plot(xx_smooth,gaussian(A_noise,xx_smooth)','k');
leg{2} = sprintf('Noise Fit, \\sigma  =%1.1f mVns',A_noise(3));
nn_noise = gaussian(A_noise,xx)';
nn2 = nn-nn_noise;
nn2(inrange(xx,xx(1),xx(w_noise)+noise_R)) = 0;

% Fit SPHE
SPHE_R = 1.5;
[h w] = max(nn2);
fit_sphe_range = inrange(xx,xx(w)/SPHE_R,xx(w)*SPHE_R);
[A Y] = FH_Gauss(xx(fit_sphe_range),nn(fit_sphe_range));
plot(xx,gaussian(A,xx)','r','linewidth',2);
% plot(xx(fit_sphe_range),nn(fit_sphe_range),'g')
gain = mVns_to_gain(A(2));
gain_error = mVns_to_gain(A(3)/sqrt(sum(nn(fit_sphe_range))));
resolution = A(3)/A(2);
error_mean = A(3)/sqrt(sum(nn(fit_sphe_range)));
error_sigma = A(3)/sqrt(2*(sum(nn(fit_sphe_range))-1));
error_resolution = sqrt((error_mean/A(2))^2 + (error_sigma/A(3))^2);

leg{3} = 'Sphe Peak';

set(handles.SpheAreaValue,'String',sprintf('%2.2f',A(2)));
set(handles.GainValue,'String',sprintf('%2.2fe6',mVns_to_gain(A(2))/1e6));
set(handles.ResolutionValue,'String',sprintf('%2.1f',A(3)/A(2)*100));

legend(leg)

%% Manual Fitting
function FitSphe_Callback(hObject, eventdata, handles)

left = str2double(handles.leftedge);
right = str2double(handles.rightedge);

axes(handles.SpheSpectrumPlot);
cut = inrange(handles.xx,left,right);
hold on
plotlogstairs(handles.xx(cut),handles.nn(cut),'r')
[A Y] = FH_Gauss(handles.xx(cut),handles.nn(cut));
plot(handles.xx(cut),Y,'k','linewidth',3)

axis(handles.ax)

set(handles.SpheAreaValue,'String',sprintf('%2.2f',A(2)));
set(handles.GainValue,'String',sprintf('%2.2fe6',mVns_to_gain(A(2))/1e6));
set(handles.ResolutionValue,'String',sprintf('%2.1f',A(3)/A(2)*100));

% Update handles structure
guidata(hObject, handles);

function LeftEdgeValue_Callback(hObject, eventdata, handles)
handles.leftedge = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

function RightEdgeValue_Callback(hObject, eventdata, handles)
handles.rightedge = get(hObject,'String');

% Update handles structure
guidata(hObject, handles);

function FitSpheManuallyCheckBox_Callback(hObject, eventdata, handles)


%% Submit to LUG
function SubmitToLUGButton_Callback(hObject, eventdata, handles)

[login pass] = logindlg;




%% Function Creation Routines

function filepath_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function SpheAreaValue_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function LeftEdgeValue_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function RightEdgeValue_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function filenamefield_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function ChTextEdit_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


