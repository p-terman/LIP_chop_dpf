function status = Event_Classification(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)

% status = Event_Classification(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% DESCRIPTION:
%    This module selects golden single scatter events and flags them as such with an array called golden, whith all entries either set to 0 or 1
%    For easy slection of the participating S1 and S2, the helper RQ selected_s1_s2 was created, which reads 1 if the particular pulse has been part of the consideration of bein a single scatter or mulitple scatter.
%
%   Example usage to select the S2 of golden single scatter use: golden==1 && pulse_classification==2 && selected_s1_s2==1
%   Example usage to select the S2s of a multiple scatter use: multiple==1 && pulse_classification==2 && selected_s1_s2==1
%
% Required RQs:
%
% Versioning:

%
%% THIS STUFF YOU ALWAYS NEED
% Load .rq file

dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);

settings.evt_settings = dp.admin.evt_settings;
settings.daq_settings = dp.admin.daq_settings;
settings.filename_prefix = dp.admin.filename_prefix;

event_number = dp.event_number;
livetime = dp.admin.livetime;

%% Bookkeeping - access parameters

myname = 'Event_Classification';
fprintf('\n\n *** Starting module %s\n',myname);

dp_settings_xml = XMLReader_framework(data_processing_xml_path);
lug_iqs_xml = XMLReader_framework(iq_xml_path);

module_names = {dp_settings_xml.data_processing_settings.module.module_name};
index_temp = strfind(module_names,myname);
index_module = find(not(cellfun('isempty', index_temp)));

mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters;

max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;



%% Initialize variables


N = length(dp.event_number);
pulse_event_size = [max_num_pulses N];

dp.golden = zeros(pulse_event_size, 'uint32');
dp.selected_s1_s2 =   zeros(pulse_event_size, 'uint32');
dp.multiple = zeros(pulse_event_size, 'uint32');

s2threshold = mymodule_settings.s2threshold;
echopercent = mymodule_settings.echopercent; % allowance for S2 echo pulses
echotime = mymodule_settings.echotime; %maximum time delay for S2 echo pulses

%% Loop per event and assign S1 S2 pairs

for ii = 1:N

    % find the indices for S1 and S2 pulses above the decided threshold
    s1_inds = find(dp.pulse_classification(:,ii) == 1);
    s2_inds = find(dp.pulse_classification(:,ii) == 2 & dp.pulse_area_phe(:,ii)>=s2threshold); % all S2 above defined s2 threshold (xml file input)

    % how many pulses of each type are there?
    num_s1 = length(s1_inds);
    num_s2 = length(s2_inds);
    
    % are there one or more pulses of each type?
    has_s1 = num_s1 > 0;
    has_s2 = num_s2 > 0;
    
    if has_s1 && has_s2
        
        % find the first s1
        first_s1_ind = find(dp.pulse_classification(:,ii) == 1,1);
        s2_after_s1 = find(     dp.pulse_classification(:,ii) == 2 & ...                   % s2 & ...
                       dp.aft_t0_samples(:,ii) > dp.aft_t0_samples(first_s1_ind,ii) & ...  % after the first s1 & ...
                       dp.pulse_area_phe(:,ii)>=s2threshold);                              % above given s2threshold
        num_s2_after_s1 = length(s2_after_s1);


        if num_s2_after_s1 > 0

            first_s2_ind = s2_after_s1(1);

            %are there more than 1 S1 before the first S2
            allearly_s1_ind = find(     dp.pulse_classification(:,ii) == 1 & ...                   % s1 & ...
                       dp.aft_t0_samples(:,ii) < dp.aft_t0_samples(first_s2_ind,ii)        );  % before the first s2
            num_s1_before_s2=length(allearly_s1_ind);


            if num_s1_before_s2 == 1 % 1 S1 before S2


                % loop over all s2 after the first s1 and apply S2 echo allowance if within the timw window of a typical S2 echo

                validS2_ind = first_s2_ind;

                for ss = 2:num_s2_after_s1
                
                    if dp.pulse_area_phe(s2_after_s1(ss),ii)>dp.pulse_area_phe(first_s2_ind,ii)*echopercent | (dp.aft_t0_samples(s2_after_s1(ss),ii)-dp.aft_t0_samples(first_s2_ind,ii))>echotime

                        validS2_ind = [validS2_ind s2_after_s1(ss)];

                    end
                end

                num_validS2 = length(validS2_ind);

                % golden single scatter event selection

                if num_validS2 == 1

                    dp.golden(:,ii)=1;
                    dp.selected_s1_s2(first_s1_ind,ii)=1;
                    dp.selected_s1_s2(first_s2_ind,ii)=1;
                end

                %looking for golden multiple scatter

                if num_validS2 > 1

                    %are there s1s after the first S2
                    s1_after_s2 = find(     dp.pulse_classification(:,ii) == 1 & ...                   % s1 & ...
                                       dp.aft_t0_samples(:,ii) > dp.aft_t0_samples(first_s2_ind,ii)        );  % before the first valid s2
                    num_s1_after_s2=length(s1_after_s2);

                    if num_s1_after_s2 == 0

                        dp.multiple(:,ii)=1;
                        dp.selected_s1_s2(first_s1_ind,ii)=1;

                        for ss = 1:num_validS2

                            dp.selected_s1_s2(validS2_ind(ss),ii)=1;

                        end

                    end

                    if num_s1_after_s2 > 0

                        if dp.aft_t0_samples(validS2_ind(num_validS2),ii)<dp.aft_t0_samples(s1_after_s2(1),ii)

                            dp.multiple(:,ii)=1;
                            dp.selected_s1_s2(first_s1_ind,ii)=1;

                            for ss = 1:num_validS2

                                dp.selected_s1_s2(validS2_ind(ss),ii)=1;

                            end
                        end
                    end
                end
            end
        end
    end
end


%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, event_number, livetime); % Should add livetime input at the end, remove event_number as settings field



