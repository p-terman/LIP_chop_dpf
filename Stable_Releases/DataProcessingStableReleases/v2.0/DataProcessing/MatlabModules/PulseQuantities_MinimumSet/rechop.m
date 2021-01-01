function status = rechop(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path, pathbase)
%this is used to check that all of the pulse is being taken from chopped pulses
%only use with chopped pulse dpf

dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);
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
total_pulse_area = sum(dp.pulse_area_phe);
total_evt_area = dp.full_evt_area_phe;
max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;


pulse_area_to_evt_area = total_pulse_area ./ total_evt_area;

re_process = find(pulse_area_to_evt_area < 0.75);

if ~isempty(re_process)
    detector_sample_length= 32000;% 20000 samples= 200ms allowed drift time, which is the whole of the 'real length' of the detector
    parts = double(max_num_pulses);
    for evt = re_process
           for k = 1 : parts
                
                %this was for the 'equal break' version
                new_pulse_start(k) = dp.pulse_start_samples(1, evt) + ((k-1)/parts)*detector_sample_length;
                new_pulse_end(k) = dp.pulse_start_samples(1, evt) + detector_sample_length - ((parts-k)/parts)*detector_sample_length;
                %this should divide the found start and end times equally, into the number of 'parts' 
                new_index_kept_sumpods(k)=1;
           end
 
        dp.pulse_end_samples(:,evt) = new_pulse_end;
        dp.pulse_start_samples(:,evt) = new_pulse_start;    
        dp.index_kept_sumpods(:,evt) = new_index_kept_sumpods;
        dp.num_pulses_found(evt) = uint32(sum(dp.index_kept_sumpods(:,evt)==1));% this will help visualux

    end
end
 
%%
status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, event_number, livetime); % Should add livetime input at the end, remove event_number as settings field
%we have saved the new dp file and will now re-run pulsetiming and
%pulsequantities modules
system(['/global/homes/p/pterman/LUXcode/Stable_Releases/v2.0/DataProcessing/CppModules/bin/PulseTiming_HeightTiming ' filename_evt ' ' data_path_evt ' ' filename_rq ' ' pathbase ' 5 ' data_processing_settings_path ' ' lug_iqs_xml_file ])
status = PulseQuantities_MinimumSet(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);

