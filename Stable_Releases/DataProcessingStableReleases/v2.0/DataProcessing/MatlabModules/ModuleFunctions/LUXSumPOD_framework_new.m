function cvt_struct = LUXSumPOD_framework_new(cvt_struct,ch_map)
%
% function cvt_struct = LUXSumPOD_framework_new(cvt_struct,ch_map)
%
% INPUTS:
%         cvt_struct - Output from LUXCVTLoader
%                        Requires LUXCalibratePulses to be run prior to this function.
%                        i.e. cvt_struct.ch.pod_data_phe must exist
%             ch_map - [optional] the vector with the channels to use
%
% OUTPUTS:
%         cvt_struct - Same as input, but adds:
%           sumpod_data_phe_per_sample  - Summed POD across all channels, in phe
%           sumpod_time_samples         - Time vector for above
%
%
% 20120216 CHF - Created
% 20120222 CHF - Fixed POD chunking algorithm
% 20120306 CHF - Fixed issue with timing logic
% 20121108 CHF - Fixed issue with small pulses appearing
% 20121213 CHF - Improved performance by x3 by taking sparse matrix out of
%                loop
% 20121213 CHF - Major improvement in algorithm - using imaginary numbers
%                to differentiate from data zeros in partitioning pulses
% 20130212 CHF - Large number of changes.
%                * removed chsum from structure
%                * renamed pod to podsum
%                * renamed data_phe to sumpod_data_phe_per_sample, tt to podsum_time_samples
%                * changed pod_channels, pod_numbers to sumpod_start_pod and sumpod_end_pod
% 20130217 CHF - Major changes - now using new, simplified event_structure
%                format, which makes summing pods much easier, and does not
%                need to keep track of which pod was where.
% 20130322 CHF - Compressed code to ~4 actual lines for summing pods.
%                Had a mindfuck moment and figured out how to do this very
%                elegantly.
% 20130430 pfs - implemented baseline thresholding, then realized it had to be part of
%                LUXCalibratePulses. oops, heh heh
% 20140501 CHF - Fixed bug where event can be considered a "good event" and still not have any pulses in the ch_map! 
% 20140815 CHF - Fixed missing 'end' in if statement
% 20140826 CHF - Added optional argument ch_map to ignore certain PMTs
% 20140901 AL  - Changed module name, to keep new and old version in DP2.0. This new version will
% 		 only be used by the Kr S1 module for now
%%

% Shortcut
if ~exist('ch_map','var')
    ch_map = 1:122;
end

% Only get index of events that are not flagged as empty
good_evts = find(cell2mat({cvt_struct.empty}) == 0);

% For each event with something in it
for evt = good_evts
   if mean([cvt_struct(evt).ch(ch_map).empty]) < 1
    if sum(~[cvt_struct(evt).ch(ch_map).empty]) > 0
    data = cat(1,cvt_struct(evt).ch(ch_map).pod_data_phe_per_sample);
	% DO NOT USE THIS (it would need to be implemented per channel) -pfs
	%data(data<cvt_struct(1).thr & data>-cvt_struct(1).thr) = 0; % symmetric threshold always applied
    times = cat(2,cvt_struct(evt).ch(ch_map).pod_time_samples)';
    
    inds = times - min(times) + 1;
    
    times_with_data = unique(times);
    inds_with_data = times_with_data - min(times) + 1;
    
    sumpod_data_phe_per_sample = accumarray(inds,data);
    sumpod_time_samples = min(times):max(times);

    cvt_struct(evt).sumpod_data_phe_per_sample = sumpod_data_phe_per_sample(inds_with_data);
    cvt_struct(evt).sumpod_time_samples = sumpod_time_samples(inds_with_data);

	%if 1 % additional code to make thresholded sumpod -pfs
	%	data_thr = data;
	%	data_thr(data_thr<cvt_struct(1).thr & data_thr>-cvt_struct(1).thr) = 0; % symmetric threshold
	%	sumpod_data_thr_phe_per_sample = accumarray(inds,data_thr);
	%	cvt_struct(evt).sumpod_data_thr_phe_per_sample = sumpod_data_thr_phe_per_sample(inds_with_data);
	%end
    
    % Sanity check
    if 0
        figure(432); clf; 
        h1 = subplot(2,1,1); 
        plot(cvt_struct(evt).sumpod_time_samples,cvt_struct(evt).sumpod_data_phe_per_sample,'k.-');
        
        h2 = subplot(2,1,2);
        plot(times,data,'.-');
        
        linkaxes([h1 h2],'xy')
        keyboard
    end
    
    ind_edges = find(diff(inds_with_data) > 1)';
    sumpod_end_samples =  cvt_struct(evt).sumpod_time_samples([ind_edges end]);

    cvt_struct(evt).sumpod_start_samples = cvt_struct(evt).sumpod_time_samples([1 ind_edges+1]);
    cvt_struct(evt).sumpod_length_samples = sumpod_end_samples - cvt_struct(evt).sumpod_start_samples + 1;
    end
   end
end % evt
