function varargout = LUXPMTHealthMonitor(varargin)
% LUXPMTHEALTHMONITOR M-file for LUXPMTHealthMonitor.fig
%      LUXPMTHEALTHMONITOR, by itself, creates a new LUXPMTHEALTHMONITOR or raises the existing
%      singleton*.
%
%      H = LUXPMTHEALTHMONITOR returns the handle to a new LUXPMTHEALTHMONITOR or the handle to
%      the existing singleton*.
%
%      LUXPMTHEALTHMONITOR('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in LUXPMTHEALTHMONITOR.M with the given input arguments.
%
%      LUXPMTHEALTHMONITOR('Property','Value',...) creates a new LUXPMTHEALTHMONITOR or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before LUXPMTHealthMonitor_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to LUXPMTHealthMonitor_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES
%
% 20091028 - JRV




% Edit the above text to modify the response to help LUXPMTHealthMonitor

% Last Modified by GUIDE v2.5 29-Oct-2009 22:18:24

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @LUXPMTHealthMonitor_OpeningFcn, ...
                   'gui_OutputFcn',  @LUXPMTHealthMonitor_OutputFcn, ...
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


% --- Executes just before LUXPMTHealthMonitor is made visible.
function LUXPMTHealthMonitor_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to LUXPMTHealthMonitor (see VARARGIN)

% Choose default command line output for LUXPMTHealthMonitor
handles.output = hObject;

handles.local_acquisition = 0;


% change lug table for upload

% handles.lug_table = 'ap_software_testing';
handles.lug_table = 'AP';


% Update handles structure
guidata(hObject, handles);

% UIWAIT makes LUXPMTHealthMonitor wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = LUXPMTHealthMonitor_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
guidata(hObject, handles); %updates the handles



function calcAPbutton_Callback(hObject, eventdata, handles)

try
    handles = rmfield(handles, {'filename', 'PMTList', 'plot_ch', 'range', 'APR','num_pC_main','ymean','tus','info_to_LUG','PMTbias','num_ch'});
catch
end

if (iscell(handles.Datasets) && ~strcmp(handles.Datasets{1},'.')...   % error protection in case of no dataset selection
        && ~strcmp(handles.Datasets{1},'..') && ~isempty(get(handles.filepath,'String')))
    
    set(handles.LoadingStatus,'String','Loading...'); drawnow;

    if iscell(handles.Datasets)
        filename = handles.Datasets;
    else
        filename = {handles.Datasets};
    end
    %handles.pathbase = get(handles.filepath,'String');

    
    
    % Clear axes
    axes(handles.plotAP);
    cla(handles.plotAP,'reset');
   
    
   if handles.local_acquisition == 1
       % If the acquisition was done using only 1 channel, and was stored in the local acquisition
       % table in the PMT database, then do this
       % 20091029 CHF
       
       [PMT_serial PMT_bias_V amp] = LUXGetLocalAcquisitionSettings(filename{1});

       Ch = 1;
       
       handles.num_ch = 1;
       handles.PMTbias{Ch} = PMT_bias_V;
       handles.PMTList{Ch} = PMT_serial;
       
       [num_pC_main APR tus ymean tot_sum range APratio] = LUXGetAP(handles.pathbase,filename{1},Ch,handles.PMTbias,amp);
       
   else

       channels = LUXGetChannels(filename{1});

       for ii = 1:numel(channels)
           handles.PMTList{ii} = channels(ii).serial;
       end

       amp = LUX01LoadAmp(filename{1}); %moved here from LUXGetAP, now an input - 20091029 CHF

       handles.num_ch = length(handles.PMTList);
       
       for Ch = 1:handles.num_ch                  % cycle though all channels in dataset, calculating the PMT AP ratios
           handles.PMTbias{Ch} = LUX01LoadPMTHV(filename{1},Ch);
       end
       
       [num_pC_main APR tus ymean tot_sum range APratio] = LUXGetAP(handles.pathbase,filename{1},handles.num_ch,handles.PMTbias,amp);

       
   end
    
   
   
    handles.tus = tus;
    handles.ymean = ymean;
    handles.num_pC_main = num_pC_main;
    handles.APR = APR;
    handles.range = range;
    handles.filename = filename{1};

    
    info_to_LUG.num_pC_main = num_pC_main;
    info_to_LUG.APR = APR;
    info_to_LUG.APratio = APratio;
    info_to_LUG.tus = tus;
    info_to_LUG.ymean = ymean;
    info_to_LUG.range = range;
    info_to_LUG.filename = filename{1};
    info_to_LUG.PMTList = handles.PMTList;
    handles.info_to_LUG = info_to_LUG;
    
    handles.plot_ch = 1;
    LUXPlotAP(tus,ymean{handles.plot_ch},handles.PMTList{handles.plot_ch},num_pC_main{handles.plot_ch},APratio{handles.plot_ch},filename{1},handles.PMTbias{handles.plot_ch}, range{handles.plot_ch},handles.plot_ch);

    set(handles.LoadingStatus,'String',''); drawnow;
else
    set(handles.LoadingStatus,'String','Please select a dataset!'); drawnow;
end


guidata(hObject, handles); %updates the handles


function DatasetList_Callback(hObject, eventdata, handles)



DatasetSelection = get(hObject,'Value');
DatasetList = get(hObject,'String');

handles.Datasets = [];


if ~isempty(DatasetList) && ~isempty(DatasetSelection)
    handles.Datasets{1} = DatasetList{DatasetSelection(1)};
    if numel(DatasetSelection) == 1
        set(handles.description,'String',handles.commentList(DatasetSelection));
    else
        set(handles.description,'String','Please select only one dataset!');
    end
else
    set(handles.LoadingStatus,'String','Please select a dataset!'); drawnow;
end

guidata(hObject, handles); %updates the handles




function DatasetList_CreateFcn(hObject, eventdata, handles)

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




function filepath_Callback(hObject, eventdata, handles)

set(handles.DatasetList,'Value',1);
handles.pathbase = get(hObject,'String');

[datasetList commentList] = LUXGetAPdatasetList(handles.pathbase);
set(handles.DatasetList,'String',datasetList)
set(handles.LoadingStatus,'String','')
handles.pathbase = get(handles.filepath,'String');
handles.commentList = commentList;

guidata(hObject, handles); %updates the handles




function filepath_CreateFcn(hObject, eventdata, handles)

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in submitLUG.
function submitLUG_Callback(hObject, eventdata, handles)
% hObject    handle to submitLUG (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

info_to_LUG = handles.info_to_LUG;
if isfield(handles,'comment')
    info_to_LUG.comment = handles.comment;
end
info_to_LUG.PMTbias = handles.PMTbias;

% get entry_group_id from LUG and increment it
xmlsettings = XMLReader('defaultLUGSettings.xml');
gid_query = sprintf('select entry_group_id from lug.%s',handles.lug_table);
[result_gid query_time_gid url_alive_gid] = LUGQuery(gid_query,xmlsettings);
result_gid = regexprep(result_gid, '\n', ' ');
result_gid = sscanf(result_gid{1},'%f');
gid = max(result_gid);

if isempty(gid)
    gid = 1;
else
    gid = gid+1;
end
info_to_LUG.gid = gid;

[username password] = logindlg;
info_to_LUG.username = username;
info_to_LUG.password = password;

% submit data to lug
if ~isempty(info_to_LUG.username) && ~isempty(info_to_LUG.password)
    for ii = 1:handles.num_ch
        LUXSubmitAP(info_to_LUG,ii,handles.lug_table,handles.local_acquisition)
    end
    set(handles.LoadingStatus,'String','Submitted');
else
    msgbox('Please input both a username and a password!','Error');
end


guidata(hObject, handles);


% --- Executes on button press in dirBrowse.
function dirBrowse_Callback(hObject, eventdata, handles)
% hObject    handle to dirBrowse (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


set(handles.DatasetList,'Value',1);
save_pathbase = get(handles.filepath, 'String');
handles.pathbase = uigetdir; % browse for directory

if handles.pathbase == false; % if no directory is selected keep current directory
    handles.pathbase = save_pathbase;
end

set(handles.filepath,'String',handles.pathbase); drawnow;


 [datasetList commentList ] = LUXGetAPdatasetList(handles.pathbase);
 set(handles.DatasetList,'String',datasetList);
 set(handles.LoadingStatus,'String','');
 handles.pathbase = get(handles.filepath,'String');
 handles.commentList = commentList;
 
guidata(hObject, handles);


% --- Executes on button press in previous_plot.
function previous_plot_Callback(hObject, eventdata, handles)
% hObject    handle to previous_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hold off;
%handles.info_to_LUG.APratio
if handles.plot_ch > 1
    handles.plot_ch = handles.plot_ch - 1;
    LUXPlotAP(handles.tus,handles.ymean{handles.plot_ch},handles.PMTList{handles.plot_ch},handles.num_pC_main{handles.plot_ch},handles.info_to_LUG.APratio{handles.plot_ch},handles.filename,handles.PMTbias{handles.plot_ch}, handles.range{handles.plot_ch},handles.plot_ch);
    guidata(hObject, handles);
else
    handles.plot_ch = 1;
end


% --- Executes on button press in next_plot.
function next_plot_Callback(hObject, eventdata, handles)
% hObject    handle to next_plot (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
hold off;
%handles.info_to_LUG.APratio
if handles.plot_ch < handles.num_ch
    handles.plot_ch = handles.plot_ch + 1;
    LUXPlotAP(handles.tus,handles.ymean{handles.plot_ch},handles.PMTList{handles.plot_ch},handles.num_pC_main{handles.plot_ch},handles.info_to_LUG.APratio{handles.plot_ch},handles.filename,handles.PMTbias{handles.plot_ch}, handles.range{handles.plot_ch},handles.plot_ch);
    guidata(hObject, handles);
else
    handles.plot_ch = handles.num_ch;
end



function datasetComment_Callback(hObject, eventdata, handles)
% hObject    handle to datasetComment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.comment = get(hObject, 'String');
guidata(hObject, handles);


% Hints: get(hObject,'String') returns contents of datasetComment as text
%        str2double(get(hObject,'String')) returns contents of datasetComment as a double


% --- Executes during object creation, after setting all properties.
function datasetComment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to datasetComment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function local_acquisition_Callback(hObject, eventdata, handles)

handles.local_acquisition = get(hObject,'Value');
guidata(hObject, handles); %updates the handles
