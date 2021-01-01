function event_struct = LUXCalibratePulses_framework(event_struct, mVns_per_phe,amp_gain)
% event_struct = LUXCalibratePulses_framework(event_struct, mVns_per_phe,amp_gain)
% Inputs:
%   event_struct - output of LUXEventLoader_framework
%   mVns_per_phe - calibration for PMTs (vector)
%   amp_gain - total area gain of analog chain (preamp * postamp area)
%
% Output:
%   event_struct - same as input, but now pulse_data_phe is added
%
%
%
% 2010-02-02 JJC created - chopped out of LUXSumPOD
% 2011-05-19 JJC, CHF - ch_map being taken from event_struct(1)
% 2011-06-06 LdV - Cleaned up a bit; instead of a generic vargin, now we use
%   explicit inputs in the function declaration; removed references to ch_map
%   and the settings structure as possible inputs - that is, they were considered
%   valid inputs before, but they replaced in the code, so that the input was
%   never used, so I removed them as possible inputs to make the code less
%   confusing;  if mVns_per_phe and amp_gain are not given in the input, it will
%   run LUXSuperLoader to get them;
% 2012-02-17 JJC modified for event_struct to have .pod instead of .pulse
% 2013-02-12 CHF Changed pod_data_phe to pod_data_phe_per_sample
% 2013-02-16 CHF Now using simplified event_struct fields
%                Removed LUXSuperLoader usage of mVns_per_phe and amp_gain.
%                That's just asking for trouble.
% 2013-04-30 pfs added symmetric, universal (all channels same) baseline thresholding

if isfield(event_struct(1),'ch_map')
    %ch_map = event_struct(1).ch_map;
	ch_map = 1:122; % this should not be hard-coded, but the previous line is not good either!
	fprintf('warning: ch_map is hard-coded (LUXCalibratePulses_framework.m)\n');
else
    fprintf('channel map not defined in event_struct\n');
    fprintf('assuming channel map is 1:number of channels\n');
    ch_map = 1:length(event_struct(1).ch);
end

conversion_factor_mV_to_phe_per_10ns = 10 ./(mVns_per_phe .* amp_gain) ;
for evt=1:length(event_struct)
    if ~event_struct(evt).empty
        for ii = 1:length(ch_map)
            ch = ch_map(ii);
            if event_struct(evt).ch(ch).empty == 0 % if it's not empty
				data = event_struct(evt).ch(ch).pod_data_mV;
				data(data<event_struct(1).thr & data>-event_struct(1).thr) = 0; % symmetric threshold always applied
                event_struct(evt).ch(ch).pod_data_phe_per_sample = ...
                	data ...
                	.* conversion_factor_mV_to_phe_per_10ns(ii);
            end
        end
    end
end




