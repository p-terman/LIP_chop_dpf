function return_status = LUXFirstPass(filename_prefix, filenumber_list, data_path, output_path, analysis_xml_file, options)
% return_status = LUXFirstPass(filename_prefix, filenumber_list, data_path, output_path, analysis_xml_file, options)
%
% S1S2 analysis of Struck DAQ data
% Reads in DAQ binary output and creates RQs (S1S2 data)
% 'REEF' means 'Really Equal Equal Footing' -- all pulses treated as equal
%
% Inputs:
%  filename_prefix: string with dataset name
%  filenumber_list: files to process (e.g. 1:1000)
%  data_path: string with path to data
%  output_path: path where REEF files are to be stored (REQUIRED)
%  options: structure with following optional fields:
%   .flag_overwrite: set to 1 to overwrite existing RQs (default: 0)
%   .flag_write_bin: set to 0 to write .mat files instead of .rq1 files
%    (default: 1)
%   .flag_debug: set to 1 to turn on user-interactive mode with verbose
%    output (default: 0)
%   .flag_log_file: set to 1 to turn on log file (saved in RQ folder)
%    (default: 0)
%	.flag_use_lug_settings: set 1 to use lug_settings settings, 0 to use hard coded in analysis.xml file
%	 (default: 1)
%
% Outputs:
%  return_status:
%   returns -1 if existing RQs for the given dataset are
%    located and 'flag_overwrite' is false
%   returns +1 on input arguments error or on failure to find traces file
%   returns  0 otherwise
%
% For older comments see the file LUX01REEFrunner_v2.m
% v3.0 2009-09-04 JJC - start rewriting to use POD mode
% v3.1 2009-10-06 JJC - tag pulses that are bipolar (noise ringing) as whichpeak = 3.
%  also add RQ peak_height_mV
% v3.2 2009-10-23 JJC - add RQ whichpulse - which pulse in the event (sum) the peak belongs to.
% v4.0 2010-01-27 JJC, DCM - rewriting as a series of smaller functions;
%  trying to trim out unnecessary code
% v5.0 2010-04-19 JJC - incorporating new file format with xml header
%
%

% YOU MUST update the version field ('verzsh', bc 'version' is a MATLAB reserved variable)
% if you modify ANYTHING about the output of this function
%


return_status = 1;

%% Function constants

verzsh = 'v5.0';

% REEF Settings
warning_state = warning('off','MATLAB:divideByZero');

%% Create options variables
if nargin > 5
    field_names = fieldnames(options);
    for ii=1:length(field_names)
        eval([field_names{ii} ' = options.' field_names{ii} ';']);
    end
end

if ~exist('flag_overwrite','var')
    flag_overwrite = false;
end
if ~exist('flag_write_bin','var')
    flag_write_bin = true;
end
if ~exist('flag_debug','var')
    flag_debug = false;
end
if ~exist('flag_log_file','var')
    flag_log_file = false;
end
if ~exist('flag_use_lug_settings','var')
    flag_use_lug_settings = true;
end

%% Function inputs
if isempty(filenumber_list)
    filenumber_list = 1:length(dir([data_path, filesep, '*.evt'])) ;
    
    if isempty(filenumber_list)
        fprintf('There are no .evt files in %s\n',[data_path filesep filename_prefix]);
        return
    end
end

% Using output_path input variable now -- REQUIRED
if ~exist('output_path','var') || isempty(output_path)
    disp('output_path must be specified');
    disp('See help LUXFirstPass');
    return
end

% Check if output folder exists
if ~exist(output_path,'dir')
    mkdir(output_path);
end

if flag_write_bin
    filename_output_extension = '.rq1';
else
    filename_output_extension = '.rqm.mat';
end


%% Analysis/Dataset Settings
% s is the structure returned by the SuperLoader, containing both
% acquisition and analysis settings

%s = LUXSuperLoader_v4(filename_prefix, data_path, [], analysis_xml_file, flag_use_lug_settings);



job_start_time = clock;


%%% Loop over all files in this job %%%
options.load_settings=1;
options.query_lug=1;
filename1 = [filename_prefix sprintf('_f%09d.evt',filenumber_list(1))];
s = LUXSuperLoader(filename1,data_path,options);
s.ana_settings = XMLReader(analysis_xml_file);
for filenumber = filenumber_list
    % Set output file name
    filename = [filename_prefix sprintf('_f%09d.evt',filenumber)];
    output_filename = [filename_prefix sprintf('_f%09d.rq1',filenumber)];
    if flag_debug; fprintf('output filename = %s\n', output_filename); end
    % check for existing file if flag_overwrite is false
    if ~flag_overwrite && exist([output_path filesep filename],'file')
        return_status = -1;
        fprintf('Output file already exists and flag_overwrite is false -- skipping file %d.',filenumber);
        continue;
    end
    %[s dual_fid] = LUXSuperLoader_v5(filename, data_path, options);
    %options.load_settings=0; % don't waste time loading these now that you've done it once.
    
    event_struct = LUXEventLoader(filename,data_path,[],s);
    evt_list = LUXNumberEvents(filename,data_path);
    
    %max_nb_s2s = par.itt_max_s2;% how many potential S2-like pulses we want to find
    %max_nb_pulses = par.itt_max;% how many peaks total we will look for
    %Nbefore = par.Nbefore; % how many peaks to look for in 2nd-pass
    %thr = par.thr_pikpeeks; % threshold for peak_peaks (as a percentage)
    %s1h = par.s1h_per_chan; % threshold per channel (phe/10ns) -- for coincidence test only
    %s2window = par.s2window;
    %s1window = par.s1window;
    %x50 = 1.25;
    %pretrigger = xml_settings.sis3301.global.pretrigger;
    
    % Reads the list of events in the given file
    %evt_list = LUXNumberEvents(s.acq,filenumber);
    if flag_debug; fprintf('length of evt_list = %d\n',length(evt_list)); end
    if isempty(evt_list)
        evt_list=0;
    end
    
    if evt_list==0
        nb_evts_in_file=0;
    else
        nb_evts_in_file = length(evt_list);
    end
    
    % Create & initialize variables for use in analysis
    LUXFirstPass_InitVars;
    if flag_debug; fprintf('Variables Initialized\n'); end
    
    if nb_evts_in_file==0 % don't do this if there are no events in the file
        continue;
    end
    % chunk_size = floor(nb_evts_in_file/evt_num_limit);
    % evt_num_limit = ceil(nb_evts_in_file/chunk_size);
    
    

    livetime = LUXGetLivetime(data_path, filenumber);
    if flag_debug; livetime
    
    
end

%for ss = 1:evt_num_limit + 1
% set the event range to analyze, based on evt_num_limit
evt_start = evt_list(1); % + chunk_size * (ss-1);
%if evt_start >= (evt_list(1)-1+nb_evts_in_file); break; end
%evt_end = (evt_list(1)-1) + chunk_size * ss;
%evt_end = min( evt_end , (evt_list(1)-1)+nb_evts_in_file );
evt_end = evt_list(end);
if flag_debug; fprintf('evt_start = %d\tevt_end = %d\n', evt_start, evt_end); end

% Load Data

%clear event_struct;
%event_struct_temp = LUXEventLoader(s.acq, evt_start:evt_end);
%if flag_debug; fprintf('data loaded\n'); end
%for ii=1:length(evt_list)
%    event_struct(ii) = event_struct_temp(evt_list(ii)); % remove empty events in front of struct
%end
fprintf('Summing Channels\n');
%event_traces = LUXEvent2Trace_old(xml_settings,evt_start:evt_end);
%event_areas.clock_time_sec  = event_areas.clock_time_sec+((double(event_traces.timestamps(end))-double(event_traces.timestamps(1)))*10e-9);

% Do event calibration and sum

if flag_use_lug_settings
    amp_gain = s.lug_settings.amplification.preamp * s.lug_settings.amplification.postamp.out1_area ;
    mVns_per_phe = [s.lug_settings.ch.sphe_area_mVns] ;
else
    amp_gain = s.acq.global.preamp * s.acq.global.postamp ;
    mVns_per_phe = s.ana_settings.mVns_per_phe ;
end

event_struct = LUXCalibratePulses(event_struct, mVns_per_phe, amp_gain);
if flag_debug; fprintf('Pulses Calibrated\n'); end
event_struct = LUXSumPOD(event_struct);
if flag_debug; fprintf('POD Summed\n'); end


% RUN the analysis

% ------------- Loop over events in the file -----------------------------------------

%%%%% RUNNING ANALYSIS ON EACH EVENT WINDOW %%%%%%
    %keyboard
for nn = (evt_start:evt_end)% - (evt_list(1)-1)% event #
    nnn = nn - (evt_list(1)-1);
    if event_struct(nn).empty==0 % if the event isn't empty of pulses...
        % adjust some variables for easier handling...
        clear datasum;
        for ii=1:length(event_struct(nn).chsum.pulse) % structure containing all the pulses in the channel sum, and the time vector
            %datasum(ii).data_phe = -event_struct(nn).chsum.pulse(ii).data_phe ; % flip it to positive now
            %datasum(ii).tt = event_struct(nn).chsum.pulse(ii).tt ;
            datasum(ii) = event_struct(nn).chsum.pulse(ii) ;
            datasum(ii).data_phe = -datasum(ii).data_phe;
            
        end
        
        % fill in some RQs now...
        timestamp(nnn) = event_struct(nn).timestamp;
        
        evt_vec=[];
        for ii=1:length(datasum)
            evt_vec = [evt_vec datasum(ii).data_phe']; % create one vector will all pulses mashed together
            evt_mean(nnn) = mean(evt_vec); % mean of all nonzero data - used to be p_p(1,nn)
            evt_std(nnn) = std(evt_vec); % standard deviation of all nonzero data - used to be p_p(2,nn)
            for var_thresh=1:20
                npts_above_thresh(var_thresh,nnn) = sum(evt_vec>var_thresh); % number of points above variable threshold (1:20 phe) of sum
            end
            
        end
        
        % trying some optimization (untested!) - DCM
        %tempvec = datasum.data_phe';
        %evt_vec = [tempvec]; % create one vector with all pulses mashed together
        %evt_mean(nn) = mean(evt_vec); % mean of all nonzero data - used to be p_p(1,nn)
        %evt_std(nn) = std(evt_vec); % standard deviation of all nonzero data - used to be p_p(2,nn)
        %for var_thresh=1:20
        %	npts_above_thresh(var_thresh,nn) = sum(evt_vec>var_thresh); % number of points above variable threshold (1:20 phe) of sum
        %end
        
        for ch=1:length(event_struct(nn).ch)
            % another of DCM's attempts to make things so-called "better". we'll see.
            % yeah it broke. knew it would.
            %evt_sum_per_ch(ch,nn) = sum([event_struct(nn).ch(ch).pulse(ii).pulse_data_phe]);
            
            evt_sum_per_ch(ch,nnn) = 0;
            for ii=1:length(event_struct(nn).ch(ch).pulse)
                evt_sum_per_ch(ch,nnn) = evt_sum_per_ch(ch,nnn) + sum(event_struct(nn).ch(ch).pulse(ii).pulse_data_phe); % total sum of all pulses in event by channel - used to be a00
            end
        end
        
        % Look for S2-like pulses
        LUXFirstPass_S2Finder ;
        
        % Now look for S1-like pulses
        LUXFirstPass_S1Finder ;
        
        LUXFirstPass_AssignEventRQs ;
        
    else % if event isn't empty of pulses...
        evt_empty(nnn)=1;
    end
end % end nn loop
%end % end evt_num_limit loop
% ------------- End events loop ------------------------------------------------------

%  Save variables

warning(warning_state)
return_status = 0;

%else
%    nb_evts_in_file
%end % if >0 events in file
if flag_debug
    nb_evts_in_file
end

% Now we transfer all the variables we have been calculating
% into the dp{1} structure

LUXFirstPass_CreateRQStruct ;

if ~flag_write_bin
    save(sprintf([output_path '/' filename]) , 'dp')
else
    LUXFirstPass_BinaryRQWriter;
end

elapsed_time = etime(clock,job_start_time);
filenumber_ind = find(filenumber==filenumber_list);

seconds_per_file = elapsed_time/filenumber_ind;

est_total_runtime_secs = seconds_per_file * length(filenumber_list);
est_total_runtime_string = secs2hms(est_total_runtime_secs);

est_total_time_remaining_secs = seconds_per_file * (length(filenumber_list) - filenumber_ind);
est_total_time_remaining_string = secs2hms(est_total_time_remaining_secs);

if ~isempty(filenumber_ind)
    percent_complete = filenumber_ind/length(filenumber_list)*100;
    fprintf('Total runtime so far is %s\n',secs2hms(elapsed_time));
    fprintf('Job is %2.1f%% complete (out of %d files) estimated time remaining is %s\n',percent_complete, length(filenumber_list), est_total_time_remaining_string);
    fprintf('(estimated time remaining is %d seconds)\n',est_total_time_remaining_secs);
end
end





