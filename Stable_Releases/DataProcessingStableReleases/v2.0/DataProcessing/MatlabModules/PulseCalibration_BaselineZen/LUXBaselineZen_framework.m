function event_struct = LUXBaselineZen_framework(event_struct,cutoff_samples,flatten_data)
% function event_struct = LUXBaselineZen_framework(event_struct,cutoff_samples,flatten_data)
%
% this function addresses a sagging baseline pathology that is observed primarily
% towards the tail end of S2 pulses. The pathology appears to originate from
% closely spaced (in time) pods during a continuous but sparse pulse (e.g. S2).
% Basically, in those situations, the ADC count baseline returned by the daq is off
% by a significant amount.
%
% Remedy: this function checks the time gap dT between pod n and pod n+1
% if dT<cutoff_samples, it pegs the baseline of pod n+1 to the baseline recorded for
% pod n. It MODIFIES two fields, pod_baseline_mV and pod_data_mV, shifting them to a
% new, balanced, eye-pleasing baseline.
%
% Caution: because fields in event_struct are modified by this function, it is
% imperative to always call it in the same order:
%
% >> LUXEventLoader_framework;
% >> LUXBaselineZen_framework;
% >> LUXCalibratePulses_framework;
% >> LUXSumPOD_framework;
%
% if this order is not respected, baseline shifts are integrated into (corrupt) the
% chsum field returned by LUXSumPOD
%
% Relax: no data on disk is modified, ever. ADC count data is not modified, EVER.
% No other fields in event_struct are modified.
%
% Note: I would have expected cutoff_samples ~ 32 to be appropriate, but empirically
% this is not sufficient. ~50 samples may be sufficient? Or should it be set to the
% characteristic width of S2 (~200)?
%
%
% 121224 pfs - created.
% 121226 pfs - added temporary fix for occasional zero ADC counts pathology
%              defined by ~any(event_struct(evt).ch(ch).pod(pp).pod_data)
%              see line 76
% 130212 chf - converted to LUXBaselineZen_framework. Changed event_structure field
%              names
% 130307 chf - now uses new event_struct format with all pods in the same
%              vector. To get a particular pod_data_mV, we are using a cut
%              in time to isolate it:
%               this_pod_cut (for pp),
%               next_pod_cut (for pp + 1)
% 130312 pfs - added if ~empty check on pods so we don't barf on sparse
%               events
% 130328 pfs - added zen_applied rq 
%%

if nargin<2
    cutoff_samples = 50;
end
if nargin<3
    flatten_data=false;
end
flatten_threshold = 1; % mV. should be set per channel and optimized

for evt=1:length(event_struct)
    
    event_struct(evt).adc_ppe = false(1,122); % initialize
    event_struct(evt).adc_sds = false(1,122); % initialize
	event_struct(evt).zen_applied = false(1,122); % initialize
    if ~event_struct(evt).empty
    for ch=1:122
        if ~event_struct(evt).ch(ch).empty
        n_pp = length(event_struct(evt).ch(ch).pod_start_samples);
        
        for pp=1:n_pp
            
            this_pod_cut = inrange(event_struct(evt).ch(ch).pod_time_samples,...
                event_struct(evt).ch(ch).pod_start_samples(pp)-1,...
                event_struct(evt).ch(ch).pod_start_samples(pp)+event_struct(evt).ch(ch).pod_length_samples(pp)+1);
            
            
            if pp < n_pp % can't operate on last pod due to pp+1 indexing
                
                dT = event_struct(evt).ch(ch).pod_start_samples(pp+1) ...
                    - (event_struct(evt).ch(ch).pod_start_samples(pp) + event_struct(evt).ch(ch).pod_length_samples(pp));
                
                next_pod_cut = inrange(event_struct(evt).ch(ch).pod_time_samples,...
                    event_struct(evt).ch(ch).pod_start_samples(pp+1)-1,...
                    event_struct(evt).ch(ch).pod_start_samples(pp+1)+event_struct(evt).ch(ch).pod_length_samples(pp+1)+1);
                
                
                %check what we are selecting
                if 0
                    figure(32); clf; hold on
                    plot(event_struct(evt).ch(ch).pod_time_samples,event_struct(evt).ch(ch).pod_data_mV,'k.-');
                    plot(event_struct(evt).ch(ch).pod_time_samples(this_pod_cut),event_struct(evt).ch(ch).pod_data_mV(this_pod_cut),'ro-')
                    plot(event_struct(evt).ch(ch).pod_time_samples(next_pod_cut),event_struct(evt).ch(ch).pod_data_mV(next_pod_cut),'bo-')
                    keyboard
                    
                end
                
                if dT < cutoff_samples % apply a healthy dose of zen
                    
                    event_struct(evt).ch(ch).pod_data_mV(next_pod_cut) = ...
                        event_struct(evt).ch(ch).pod_data_mV(next_pod_cut) + ...
                        event_struct(evt).ch(ch).pod_baseline_mV(pp+1) - ...
                        event_struct(evt).ch(ch).pod_baseline_mV(pp) ...
                        ;
                    % now that pod's baseline is sorted,
                    % we have to update the pod_baseline_mV field for successive pods
                    event_struct(evt).ch(ch).pod_baseline_mV(pp+1) = ...
                        event_struct(evt).ch(ch).pod_baseline_mV(pp) ...
                        ;
                	event_struct(evt).zen_applied(ch) = true; % with this logic, if any pod in the event gets zen'd, we record a true. we don't break it down to the pod level...
                end
            end
            
            % now that the baselines all hover around zero with no significant offset,
            % we can apply baseline flattening
            if flatten_data
                event_struct(evt).ch(ch).pod_data_mV( ...
                    this_pod_cut ...
                    & event_struct(evt).ch(ch).pod_data_mV < flatten_threshold ...
                    & event_struct(evt).ch(ch).pod_data_mV > -flatten_threshold ...
                    ) = 0; % flattened
            end
            
            if sum(event_struct(evt).ch(ch).pod_data(this_pod_cut)==0) ...
                    == length(event_struct(evt).ch(ch).pod_data(this_pod_cut))  % this condition, defined
                % as having all ADC counts = 0, corrupts the event sum.
                % While we look for a hardware fix, we can "fix" it in software
                % by setting the offending channel's pod_data_mV = 0 (instead of
                % the maximal, ADC-saturated value they would otherwise obtain). For the
                % analysis, its as if that PMT were off
                event_struct(evt).ch(ch).pod_data_mV(this_pod_cut) ...
                    = zeros(size(event_struct(evt).ch(ch).pod_data(this_pod_cut)));
                event_struct(evt).adc_ppe(ch) = true; % ppe => pod plateau error flag
            end
            if max(abs(diff(double(event_struct(evt).ch(ch).pod_data(this_pod_cut))))) > 16000 % physical pulses should not exceed this (i hope, needs study)
                event_struct(evt).ch(ch).pod_data_mV(this_pod_cut) ...
                    = zeros(size(event_struct(evt).ch(ch).pod_data(this_pod_cut)));
                event_struct(evt).adc_sds(ch) = true; % sds => sudden death syndrome error flag
            end
        end
        end %fi
    end  %rof
    end %fi
 end %rof

%figure(8);for ch=1:122; try;plot(double(event_struct(ii_evt).ch(ch).pod(1).pod_data),'o-');hold off;title(ch);pause;end; end
%hold on;plot(abs(diff(double(event_struct(ii_evt).ch(ch).pod(1).pod_data))),'ro-')
