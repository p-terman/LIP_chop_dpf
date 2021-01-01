function status = PulseTiming_BasicSet(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% This is a PulseTiming module that calculates pulse timing in Matlab
%
% status = PulseTiming_BasicSet(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
%
%
%
%
%
% Required RQs:
%
%
%
%
%
% Versioning:
%
%    v1.0   20140204 SS - Created
%           20140317 SS - Modified hfts to include flags so they are in time order
%
% RQ versions:
%
%
%
%% Load .rq file

status = [];

dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);

settings.evt_settings = dp.admin.evt_settings;
settings.daq_settings = dp.admin.daq_settings;
settings.filename_prefix = dp.admin.filename_prefix;

event_number = dp.event_number;
livetime = dp.admin.livetime;

%% Bookkeeping

myname = 'PulseTiming_BasicSet';
fprintf('\n\n *** Starting module %s\n',myname);

if isempty(event_number)
    fprintf('\n\n *** Skipping module (no events in file) %s\n',myname);
    return
end

dp_settings_xml = XMLReader_framework(data_processing_xml_path);
lug_iqs_xml = XMLReader_framework(iq_xml_path);

module_names = {dp_settings_xml.data_processing_settings.module.module_name};
index_temp = strfind(module_names,myname);
index_module = find(not(cellfun('isempty', index_temp)));

mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters;

max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;

%% Load .cvt file or calculate it

filename_cvt = strrep(filename_evt,'evt','cvt');

if ~exist([data_path_evt filesep filename_cvt],'file')
    fprintf('Did not find .cvt file. Running Summer Module\n');
    status = PulseCalibration_BaselineZen(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
    status = PODSummer_LUXSumPOD(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
end

fprintf('Loading sum from %s\n',filename_cvt);
[cvt_struct settings] = LUXCVTLoader_framework(data_path_evt,strrep(filename_evt,'evt','cvt'));

% We still need the mV data...
event_struct = LUXEventLoader_framework(data_path_evt, filename_evt);

%% Initialize variables

pmt_chs = 1:122;

N = length(cvt_struct);

pulse_event_size = [max_num_pulses N];
per_channel_event_size = [max_num_pulses length(pmt_chs) N];

% should make these sparse!



% Timing RQs	  

dp.aft_t0_samples                   = zeros(pulse_event_size);
dp.aft_tlx_samples                  = zeros(pulse_event_size);
dp.aft_t1_samples                   = zeros(pulse_event_size);
dp.aft_trx_samples                  = zeros(pulse_event_size);
dp.aft_t2_samples                   = zeros(pulse_event_size);

dp.hft_t0_samples                   = zeros(pulse_event_size);
dp.hft_tl0l_samples                 = zeros(pulse_event_size);
dp.hft_t50l_samples                 = zeros(pulse_event_size);
dp.hft_t1_samples                   = zeros(pulse_event_size);
dp.hft_t50r_samples                 = zeros(pulse_event_size);
dp.hft_t10r_samples                 = zeros(pulse_event_size);
dp.hft_t2_samples                   = zeros(pulse_event_size);




%% Compute RQs

fprintf('\n');

dp.event_timestamp_samples = [event_struct(:).timestamp];

for evt = 1:N
    
    if mod(evt,100) == 1
        fprintf('\n')
    end
    
    fprintf('.');
    
    % Check if fields exist and has data
    if isfield(cvt_struct(evt),'sumpod_data_phe_per_sample') && cvt_struct(evt).empty == 0
        loop_max = min([sum(isfinite(dp.index_kept_sumpods(:,evt))) max_num_pulses]);

        % for every pulse (NOT pod)
        for pp = 1:loop_max
            %dis('pp=%d,evt=%d ... index_kept=%d\n',pp,evt,dp.index_kept_sumpods(pp,evt))
            % -1 and +1 signs are needed to get the right time cut
            %pulse_cut = inrange(cvt_struct(evt).sumpod_time_samples,dp.pulse_start_samples(pp,evt)-1,dp.pulse_end_samples(pp,evt)+1);
            
            % use explicit range cuts, not "inrange" function
            pulse_cut = (cvt_struct(evt).sumpod_time_samples >= dp.pulse_start_samples(pp,evt)) ...
                      & (cvt_struct(evt).sumpod_time_samples <= dp.pulse_end_samples(pp,evt));
            pulse_data_phe = cvt_struct(evt).sumpod_data_phe_per_sample(pulse_cut);
           
            


           % Pulse area fractional timing
    if numel(pulse_data_phe) > 0
           fullBoxCsum(1) = pulse_data_phe(1); % Calculate cumulative area
           pulse_length = dp.pulse_end_samples(pp,evt) - dp.pulse_start_samples(pp,evt)+1; %pulse_length gives number of final entry in array
           for tt=1:pulse_length-1
              if tt<numel(pulse_data_phe) %
		          fullBoxCsum(tt+1) = fullBoxCsum(tt) + pulse_data_phe(tt+1);
              else 
		          fullBoxCsum(tt+1) = fullBoxCsum(tt);
              end
           end
           
           %Initiliaze timings & flags
           t0  = 1;
           tlx = 1;
           t1  = 1;
           trx = 1;
           t2  = 1;
        
	       fflag_t0  = 0;
           fflag_tlx = 0;
           fflag_t1  = 0;
           fflag_trx = 0;
           fflag_t2  = 0;
    
     
				
           for tt = 1:pulse_length
	          if fullBoxCsum(tt) >= mymodule_settings.edgeFraction*fullBoxCsum(pulse_length-1) && ~fflag_t0
		          t0 = tt;  
                  fflag_t0 = 1;
	          elseif fullBoxCsum(tt) >= mymodule_settings.txFraction*fullBoxCsum(pulse_length-1) && ~fflag_tlx 
		         tlx = tt;
                 fflag_tlx = 1;	    	       			
              elseif fullBoxCsum(tt) >= 0.50*fullBoxCsum(pulse_length-1) && ~fflag_t1
	 	         t1 = tt;
                 fflag_t1 = 1;	
              elseif fullBoxCsum(tt) >= (1-mymodule_settings.txFraction)*fullBoxCsum(pulse_length-1) && ~fflag_trx
		         trx = tt;
                 fflag_trx = 1;		
              elseif fullBoxCsum(tt) >= (1-mymodule_settings.edgeFraction)*fullBoxCsum(pulse_length-1) && ~fflag_t2 
                 t2 = tt;
                 fflag_t2 = 1;   
              end
           end
           
           % Record values
           dp.aft_t0_samples(pp,evt)  = t0  + dp.pulse_start_samples(pp,evt)-1;
           dp.aft_tlx_samples(pp,evt) = tlx + dp.pulse_start_samples(pp,evt)-1;
           dp.aft_t1_samples(pp,evt)  = t1  + dp.pulse_start_samples(pp,evt)-1;
           dp.aft_trx_samples(pp,evt) = trx + dp.pulse_start_samples(pp,evt)-1;
           dp.aft_t2_samples(pp,evt)  = t2  + dp.pulse_start_samples(pp,evt)-1;
           
           
           % Following added as a fix for when conditions for trx and t2 arent satisfied 
           % before end of pulse (generally 3 sample pulses)
           if dp.aft_trx_samples(pp,evt) < dp.aft_t1_samples(pp,evt)
              dp.aft_trx_samples(pp,evt) = dp.pulse_end_samples(pp,evt);
           end
           if dp.aft_t2_samples(pp,evt) < dp.aft_t1_samples(pp,evt)
              dp.aft_t2_samples(pp,evt) = dp.pulse_end_samples(pp,evt);
           end

     
           
           % Height Fractional Timing
           % First identify the highest point in the pulse
           maxValue = 0;
           maxIndex = 0;
           for ii = 1:pulse_length %changed from pulse_length-1 so looks to the very end of pulse
               if ii <= numel(pulse_data_phe)
	             if pulse_data_phe(ii) > maxValue
	                maxValue = pulse_data_phe(ii);
                    maxIndex = ii;
                 end
               else
                  break
               end
           end
          
          
        
           % Hft initiliazation	    
           t0   = 1;
           t10l = 1;
           t50l = 1;
           t1   = maxIndex;
           t50r = 1;
           t10r = 1;
           t2   = 1;


           hflag_t0   = 0;
           hflag_t10l = 0;
           hflag_t50l = 0;
           hflag_t1   = 0;
           hflag_t50r = 0;
           hflag_t10r = 0;
           hflag_t2   = 0;

           % Look below peak height
           for ii=maxIndex:-1:1
              if ii <= numel(pulse_data_phe)
                 if pulse_data_phe(ii) >= maxValue*0.5 && hflag_t10l ~= 1 && hflag_t0 ~= 1   % 50% rising edge
 	                t50l = ii;
                    hflag_t50l = 1;
                 end
                 if pulse_data_phe(ii) >= maxValue*0.1 && hflag_t0 ~= 1  % 10% rising edge
                    t10l = ii;
                    hflag_t10l = 1;
                 end 
                 if pulse_data_phe(ii) <= 0
                    t0 = ii;                            % 0% rising edge
                    hflag_t0 = 1;
                 end
              end
           end
           
           % Look above peak height
           for ii=maxIndex:pulse_length
              if ii <= numel(pulse_data_phe)
	             if pulse_data_phe(ii) >= maxValue*0.5 && hflag_t10r ~= 1 && hflag_t2 ~= 1  % 50% falling edge
	                t50r = ii;
                    hflag_t50r = 1;
                 end
                 if pulse_data_phe(ii) >= maxValue*0.1 && hflag_t2 ~= 1  % 10% falling edge
	                t10r = ii;
                    hflag_t10r = 1;
                 end          
                 if pulse_data_phe(ii) <= 0  || ii == pulse_length
	                t2 = ii;
                    hflag_t2 = 1; 
                 end
                
              else
                  t2 = pulse_length;
              end 
           end
        
           % Record hfts
	       dp.hft_t0_samples(pp,evt)   = t0   + dp.pulse_start_samples(pp,evt)-1;
           dp.hft_t10l_samples(pp,evt) = t10l + dp.pulse_start_samples(pp,evt)-1;
	       dp.hft_t50l_samples(pp,evt) = t50l + dp.pulse_start_samples(pp,evt)-1;
	       dp.hft_t1_samples(pp,evt)   = t1   + dp.pulse_start_samples(pp,evt)-1;
           dp.hft_t50r_samples(pp,evt) = t50r + dp.pulse_start_samples(pp,evt)-1;
           dp.hft_t10r_samples(pp,evt) = t10r + dp.pulse_start_samples(pp,evt)-1;
           dp.hft_t2_samples(pp,evt)   = t2   + dp.pulse_start_samples(pp,evt)-1;
           
            
         

      end      
    
        end % for pulse pp
  
    end % fi sumpod
    
end % for event evt


 
fprintf('Done!\n');


%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, event_number, livetime); % Should add livetime input at the end, remove event_number as settings field
