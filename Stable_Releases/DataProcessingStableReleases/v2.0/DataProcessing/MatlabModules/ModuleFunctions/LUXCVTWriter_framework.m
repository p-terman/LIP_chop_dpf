function status = LUXCVTWriter_framework( event_struct, settings, livetime, filename )
% function status = LUXCVTWriter_framework( event_struct, settings, livetime, filename )
% Write cvt file from event_struct.
% Inputs:   event_struct: data to be written. Can include
%               event_struct(:).chsum, in which case ch 137 will be written with sumpods.
%           settings: xml settings struct to be written to header.
%           livetime: contains livetime.livetime_latch_samples and livetime.livetime_end_samples in samples.
%           filename: including path. Should have .cvt extension.
%
% 2013-02-11 JJC - created
% 2013-02-17 CHF - Changes to make it work as CVT writer
% 2013-03-22 CHF - Fixed sumpod offset readout problem that was plaging CVT files. 
%                  Cleaned up the code.
%%

status = -1;
fid = fopen(filename, 'wb');

settings_string = MakeXMLString_framework(settings);
settings_string_length = numel(settings_string);

fwrite(fid, hex2dec('01020304'), 'uint32');
fwrite(fid, settings_string_length, 'uint32');
fwrite(fid, settings_string, 'int8');

% write file header
location_experiment_run_version = 0;
fwrite(fid, hex2dec('01020304'), 'uint32');
fwrite(fid, 123, 'uint32');
fwrite(fid, location_experiment_run_version, 'uint32');
fwrite(fid, length(event_struct), 'uint32');
for ii=1:length(event_struct)
    fwrite(fid, event_struct(ii).event_number, 'uint32');
    index_loc(ii) = ftell(fid);
    fwrite(fid, 0, 'uint32');
end

% write livetime header
fwrite(fid, 1, 'uint16');
fwrite(fid, livetime.livetime_latch_samples, 'uint64');
fwrite(fid, livetime.livetime_end_samples, 'uint64');

for ii=1:length(event_struct)
    % write event gid
    current_event_gid_loc = ftell(fid);
    fseek(fid, index_loc(ii), 'bof');
    fwrite(fid, current_event_gid_loc, 'uint32');
    fseek(fid, current_event_gid_loc, 'bof');
    fwrite(fid, 1234, 'uint32');
    fwrite(fid, location_experiment_run_version, 'uint32');
    fwrite(fid, event_struct(ii).event_number, 'uint32');
    
    if isfield(event_struct(ii),'sumpod_start_samples')
        max_ch = 137;
    else
        max_ch = 136;
    end
    
    fwrite(fid, max_ch, 'uint32');
    fwrite(fid, 0, 'uint32');
    
    % Trigger stuff
    %     if isfield(settings.daq_settings.sis3301.global,'read_xlm') && settings.daq_settings.sis3301.global.read_xlm == 1
    %         fwrite(fid, event_struct(ii).trigger_timestamp, 'uint64');
    %         fwrite(fid, event_struct(ii).trigseqnum, 'uint32');
    %         fwrite(fid, event_struct(ii).max_filter_response, 'uint32');
    %         fwrite(fid, event_struct(ii).max_ch_ID, 'int8');
    %         fwrite(fid, length(event_struct(ii).S1_hit_vector), 'uint16');
    %         for dd=1:length(event_struct(ii).S1_hit_vector)
    %             fwrite(fid, event_struct(ii).S1_hit_vector(dd), 'int8');
    %             fwrite(fid, event_struct(ii).S2_hit_vector(dd), 'int8');
    %         end
    %         fwrite(fid, event_struct(ii).xlm_ccheck, 'int8');
    %     end
    
    % write record format/size
    fwrite(fid, [0 0], 'uint32');
    
    % write event timestamp
    fwrite(fid, event_struct(ii).timestamp, 'uint64');
    
    for ch=1:max_ch
        
        fwrite(fid, 14, 'int8'); % bit_precision
        fwrite(fid, 2000/2^14, 'double'); % v_range_V
        fwrite(fid, 0, 'double'); % v_offset
        fwrite(fid, 1/1e8, 'double'); % time_step_s
        fwrite(fid, event_struct(1).pretrigger, 'int32'); % pretrigger_samples
        fwrite(fid, 0, 'uint32'); % temp1
        fwrite(fid, 24, 'uint32'); % n_preceding_samples
        fwrite(fid, 31, 'uint32'); % n_follogin_samples
        
        if ch == 137
            nb_pulses = length(event_struct(ii).sumpod_start_samples);
        else
            if ch <= length(event_struct(ii).ch) && isfield(event_struct(ii).ch(ch),'pod_start_samples')
                nb_pulses = length(event_struct(ii).ch(ch).pod_start_samples);
            else
                nb_pulses = 0;
            end
        end
        
        fwrite(fid, nb_pulses, 'uint32');        

        if nb_pulses > 0
            
            if ch == 137 && ~isempty(event_struct(ii).sumpod_data_phe_per_sample)
                
                fwrite(fid, event_struct(ii).sumpod_start_samples, 'int32');
                fwrite(fid, event_struct(ii).sumpod_length_samples, 'uint32');
                
                fwrite(fid, zeros(1,nb_pulses), 'uint32'); % Add fake baseline of zeros
                fwrite(fid, event_struct(ii).sumpod_data_phe_per_sample, 'single');
            
            else
            
                fwrite(fid, event_struct(ii).ch(ch).pod_start_samples, 'int32');
                fwrite(fid, event_struct(ii).ch(ch).pod_length_samples, 'uint32');

                % Add fake baseline of zeros
                fwrite(fid, zeros(1,nb_pulses), 'uint32');            
                fwrite(fid, event_struct(ii).ch(ch).pod_data_phe_per_sample, 'single');

            end % fi ch == 137
            
        end % fi nb_pulses
        
    end % for ch
    
    
    
end % for event struct

fclose(fid);
status = 1;

