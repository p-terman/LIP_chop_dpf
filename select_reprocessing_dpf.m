function status = select_reprocessing_dpf(filename_evt)

addpath(genpath(['~/LUXcode/Stable_Releases/v2.0/DataProcessing']));
pathbase= '/home/paul/LUXdata'

filename_prefix=filename_evt(1:19); %looking at rq file 
data_path_evt = [pathbase filesep filename_prefix filesep];
filename_rq = strrep(filename_evt,'evt','rq');
data_path_rq = [pathbase filesep];

lug_iqs_xml_file = '~/MATLAB/lug_iqs_new4.xml'
iq_xml_path=lug_iqs_xml_file;

data_processing_settings_path = '~/MATLAB/dp_settings_new.xml'
data_processing_xml_path =data_processing_settings_path;

%load_multi % we should put some type of input into load multi to be able
%to handle a speific rq file. ?

%dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);
dp = load([data_path_rq '/matfiles/' filename_rq '.mat']);

settings.evt_settings = dp.admin.evt_settings;
settings.daq_settings = dp.admin.daq_settings;
settings.filename_prefix = dp.admin.filename_prefix;
 
event_number = dp.event_number;
livetime = dp.admin.livetime;

%% Bookkeeping

dp_settings_xml = XMLReader_framework(data_processing_xml_path);
lug_iqs_xml = XMLReader_framework(iq_xml_path);
data_processing_settings_path = data_processing_xml_path;
lug_iqs_xml_file = iq_xml_path;

%%


if ~isfield(dp, 'pulse_end_samples')
	fprintf('Should have run Pulse Height Timing and Pulse Quantities Minimum Set !\n')
	system(['PulseTiming_HeightTiming ' filename_evt ' ' data_path_evt ' ' filename_rq ' ' pathbase ' 5 ' data_processing_settings_path ' ' lug_iqs_xml_file ])
	status = PulseQuantities_MinimumSet(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
	dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);
end


%set your constants
s1_length = 50; %samples
minimum_sample_length = 4000; %samples required to chop a pulse
max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;
parts = double(max_num_pulses); % number of parts that we will chop into
detector_sample_length= 32000;% 320ms allowed drift time, which is the whole of the 'real length' of the detector
min_pulse_area = 1000; %phd this is the minimum area required to say that this is the 'first major pulse' in event
                     %this min area will help to screen noise but keep s1 -
                     %remember, looking for large pulses, normal dpf will
                     %get the small s1s


is_long_pulse = dp.pulse_length_samples > minimum_sample_length;
re_process = find(sum(is_long_pulse));

if ~isempty(re_process) % a pulse in event is long
    %first step in the dpf
    system(['/home/paul/LUXcode/Stable_Releases/v2.0/DataProcessing/CppModules/bin/InitializeRQFile_Initialize ' filename_evt ' ' data_path_evt ' ' filename_rq ' ' pathbase ' 1 ' data_processing_settings_path ' ' lug_iqs_xml_file ])
    
    %% chop the evt file 
    event_struct = LUXEventLoader_framework(data_path_evt, filename_evt);
    evt_list = fieldnames(event_struct);
    clear new_evt        
    N_new = length(re_process);
    for i = 1:length(evt_list)
        current_evt = evt_list{i};
      
        for j = 1:N_new % go through each event
            event = re_process(j);

            if j==1
                if strcmp(current_evt, 'ch_map') || strcmp(current_evt, 'filename_prefix') ||...
                    strcmp(current_evt, 'posttrigger') || strcmp(current_evt, 'pretrigger')
                    
                    new_evt_file(j).(current_evt) = event_struct(:,1).(current_evt);
                    
                elseif strcmp(current_evt, 'filename')
                    filename_evt_new = [filename_evt(1:end-4) '_chop.evt'];
                    new_evt_file(j).(current_evt) = filename_evt_new;
                else
                    new_evt_file(j).(current_evt) = event_struct(:, event).(current_evt); 
                end
            elseif j == N_new
                if  strcmp(current_evt, 'filename_prefix') ||...
                    strcmp(current_evt, 'posttrigger') || strcmp(current_evt, 'pretrigger')
                    
                    new_evt_file(j).(current_evt) = event_struct(:,end).(current_evt);
                else
                    new_evt_file(j).(current_evt) = event_struct(:,event).(current_evt); 
                end
                
            else
                % j not 1 or last 
            new_evt_file(j).(current_evt) = event_struct(:, event).(current_evt); 
            end
        end
    end % evt_list

    %save the new file
    fprintf('Saving EVT file... ');
    filename_evt_fullpath_new = [data_path_evt filesep filename_evt_new '.mat'];
    save(filename_evt_fullpath_new, 'new_evt_file');
    chop_filename_evt = filename_evt_new;
    clear new_evt event_struct;
    
    dp_chop.admin = dp.admin;
    dp_chop.file_number =  dp.file_number(re_process);
    dp_chop.event_number = dp.event_number(re_process);
    dp_chop.source_filename = dp.source_filename;
    dp_chop.event_timestamp_samples = dp.event_timestamp_samples(re_process);
    dp_chop.luxstamp_samples = dp.luxstamp_samples(re_process);
    dp_chop.time_since_livetime_start_samples = dp.time_since_livetime_start_samples(re_process);
    dp_chop.time_until_livetime_end_samples = dp.time_until_livetime_end_samples(re_process);
             
    chop_filename_rq = strrep(filename_rq , '.rq' , '_chop.rq');
 
    status = LUXBinaryRQWriter_framework(settings, dp_chop, chop_filename_rq, data_path_rq, re_process, livetime);
    
    status = PulseCalibration_BaselineZen_split_sim(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_settings_path,lug_iqs_xml_file);
    status = PODSummer_LUXSumPOD(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_settings_path,lug_iqs_xml_file);

    %% Section in place of Transparent Rubix Cube
    % look at the finalized non-chop result
    % then find where to chop

    new_pulse_start = zeros(parts, N_new, 'int32');
    new_pulse_end = zeros(parts, N_new,  'int32');
    new_index_kept_sumpods = zeros(parts, N_new, 'uint8');
    %chopped_flag = zeros(1, N_new, 'uint8'); %flag = 0 untouched or reg
    %flag = 1 chopped pulse 

    
    for new_evt = 1: N_new %this is to count the event number in the reduced evt
        evt = re_process(new_evt); % only go through events that need reprocessing
       
        large_pulse = find(dp.pulse_area_phe(:, evt) > min_pulse_area);   %find first large pulse 
        first_large_pulse = large_pulse(1);  
                %this can only work if the first large pulse is a currently
                %saved pulse
     %   chopped_flag(new_evt) = 1;
        for k = 1 : parts 
                
             if k == 1 %imposing an s1 size
                new_pulse_start(k, new_evt) = dp.pulse_start_samples(first_large_pulse, evt);
                new_pulse_end(k, new_evt) = dp.pulse_start_samples(first_large_pulse, evt) + s1_length; 
             else
                new_pulse_start(k, new_evt) = dp.pulse_start_samples(first_large_pulse, evt) + s1_length + ((k-2)/(parts-1)) * (detector_sample_length - s1_length);
                new_pulse_end(k, new_evt) = new_pulse_start(k, new_evt) + 1/(parts-1) * (detector_sample_length - s1_length);
                %this should divide the detector equally, into the number of 'parts'
                %NOTE: there might be rounding errors, giving you
                %slighly higher or lower number of total samples
             end % if k==1
                
             if dp.index_kept_sumpods(first_large_pulse, evt)==1
                new_index_kept_sumpods(k, new_evt)=1;
             end

        end % for k
        
    end %going through all events that have long pulse


    new_num_pulses_found = uint32(sum(new_index_kept_sumpods==1));% this will help visualux
   
    %{
    dp_chop.adc_ppe = dp.adc_ppe(re_process);
    dp_chop.adc_sds = dp.adc_sds(re_process);
    dp_chop.zen_applied = dp.zen_applied(re_process);
    if isfield(dp, 'livetime_latch_samples')
        dp_chop.livetime_latch_samples = dp.livetime_latch_samples;
    end
    if isfield(dp, 'livetime_end_samples')
        dp_chop.livetime_end_samples = dp.livetime_end_samples;
    end
    %}
    dp_chop.num_pulses_found = new_num_pulses_found;
    dp_chop.pulse_start_samples = new_pulse_start;
    dp_chop.pulse_end_samples = new_pulse_end;
    dp_chop.index_kept_sumpods = new_index_kept_sumpods;
    %dp_chop.chopped_flag = chopped_flag;
    %now that we have this info let's save it

 
    status = LUXBinaryRQWriter_framework(settings, dp_chop, chop_filename_rq, data_path_rq, re_process, livetime);
    clear dp_chop; 
%{    
    %% now get the required cvt info
    filename_cvt = strrep(filename_evt,'evt','cvt');
    [cvt_struct settings] = LUXCVTLoader_framework(data_path_evt, filename_cvt);

    cvt_list = fieldnames(cvt_struct);
    clear new_cvt;
    for i = 1:length(cvt_list)
        current_cvt = cvt_list{i};

        for j = 1:ending % go through each event reprocessed event
            current_event = re_process(j);
  
            if j==1
                if strcmp(current_cvt, 'ch_map') || strcmp(current_cvt, 'filename_prefix') ||...
                    strcmp(current_cvt, 'posttrigger') || strcmp(current_cvt, 'pretrigger')
            
                    new_cvt(j).(current_cvt) = cvt_struct(:, 1).(current_cvt);
                else
                    new_cvt(j).(current_cvt) = cvt_struct(:,current_event).(current_cvt); 
                end
            elseif j == ending
                if  strcmp(current_cvt, 'filename_prefix') ||...
                    strcmp(current_cvt, 'posttrigger') || strcmp(current_cvt, 'pretrigger')
                    
                    new_cvt(j).(current_cvt) = cvt_struct(:, end).(current_cvt);
                else
                    new_cvt(j).(current_cvt) = cvt_struct(:, current_event).(current_cvt); 
                end    
            else 
            new_cvt(j).(current_cvt) = cvt_struct(:, current_event).(current_cvt); 
            end% j not 1 or end
        end %event list
    end % cvt_list

    %save the new file
    fprintf('Saving CVT file... ');
    chop_filename_cvt = strrep(filename_cvt , '.cvt' , '_chop.cvt');
    filename_cvt_fullpath_new = [data_path_evt filesep chop_filename_cvt];
 
    status = LUXCVTWriter_framework( new_cvt, settings, livetime, filename_cvt_fullpath_new );
    clear new_cvt cvt_struct;
%}
    system(['/home/paul/LUXcode/Stable_Releases/v2.0/DataProcessing/CppModules/bin/PulseTiming_HeightTiming ' chop_filename_evt ' ' data_path_evt ' ' chop_filename_rq ' ' pathbase ' 5 ' data_processing_settings_path ' ' lug_iqs_xml_file ]);

    status = PulseQuantities_MinimumSet_split_sim(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
    status = PulseQuantities_PhotonCounting_split_sim(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
    status = PulseClassifier_MultiDimensional(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
    status = S1S2Pairing_Naive(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
    status = Event_Classification(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
    status = PositionReconstruction_MercuryI (chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
    status = Corrections_PositionCorrection_split_sim(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
    status = Corrections_ApplyCorrections(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
    status = PulseQuantities_TimeSince(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path); %deals with etrains, and looks at the time since the last 'big' event
    status = AdditionalFileFormat_SaveMatFile(chop_filename_evt,data_path_evt,chop_filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
        
    %now the dp_chop is sreaved as a .mat file
    dp_chop = load([data_path_rq '/matfiles/' chop_filename_rq '.mat']);
    %dp_reg = LUXLoadRQ1s_framework(filename_rq, data_path_rq);
    dp = post_dpf_corrections(filename_evt , dp);
   
    
    dp_combined = dp;
    clear dp;
    fields = fieldnames(dp_combined);
    for ii = 1: length(fields)
        current_field = fields{ii};
        if strcmp(current_field, 'admin') || strcmp(current_field, 'source_filename') ...
                || strcmp(current_field, 'livetime_latch_samples') || strcmp(current_field, 'livetime_end_samples')...
                || ~isfield(dp_chop, current_field)
            continue
        elseif ndims(dp_combined.(current_field))>2
            dp_combined.(current_field)(:, :, re_process) = dp_chop.(current_field);
        else
            dp_combined.(current_field)(:, re_process) = dp_chop.(current_field);
        end
    end
    
    dp_combined.chopped_flag = zeros(1, length(dp_combined.file_number));
    dp_combined.chopped_flag(re_process) = 1;

    %save it
    new_filename_rqmat = strrep(filename_rq , '.rq' , '_comb_chop.rq.mat');

    matfiles_dir = dir([data_path_rq '/matfiles_comb_chop/']);
    if isempty(matfiles_dir)
        mkdir(data_path_rq, '/matfiles_comb_chop/');
    end

    save([data_path_rq '/matfiles_comb_chop/' new_filename_rqmat] , 'dp_combined');

    fprintf('Done!\n') 
end % if re_process is empty
re_process
%exit
end
