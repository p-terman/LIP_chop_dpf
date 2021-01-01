function event_struct = LUXCalibratePulses(varargin)
% event_struct = LUXCalibratePulses(event_struct, mVns_per_phe,amp_gain)
% Inputs:
%   event_struct - output of LUXEventLoader
%   mVns_per_phe - calibration for PMTs (vector)
%   amp_gain - total area gain of analog chain (preamp * postamp area)
%   
%   if 1 input, call LUXSuperLoader for settings
%   if 2 inputs, assume second is settings structure with calibrations
%
% Output:
%   event_struct - same as input, but now pulse_data_phe is added
%
%
%
% 2010-02-02 JJC created - chopped out of LUXSumPOD

%if ~exist('mVns_per_phe','var') || ~exist('amp_gain','var') || isempty(mVns_per_phe) || isempty(amp_gain)
 if nargin<2
    event_struct = varargin{1};
    fprintf('No calibrations specified. Looking in LUG...\n');
    %s = LUXSuperLoader_v5(event_struct(end).filename_prefix, [], [], [], 1);
    mVns_per_phe = [s.lug_settings.ch.sphe_area_mVns];
    amp_gain = s.lug_settings.amplification.preamp * s.lug_settings.amplification.postamp.out1_area ;

 elseif nargin<3 % second input is settings structure
    event_struct = varargin{1};
    s = varargin{2};
    mVns_per_phe = [s.lug_settings.ch.sphe_area_mVns];
    amp_gain = s.lug_settings.amplification.preamp * s.lug_settings.amplification.postamp.out1_area ;
 
 elseif nargin==3
    event_struct = varargin{1};
    mVns_per_phe = varargin{2};
    amp_gain = varargin{3};
 end
    
conversion_factor_mV_to_phe_per_10ns = 10./(mVns_per_phe .* amp_gain) ;

for evt=1:length(event_struct)
	if ~event_struct(evt).empty
		for ch=1:min(4,length(event_struct(evt).ch))
			if ~isempty(event_struct(evt).ch(ch).pulse)
				for ps=1:length(event_struct(evt).ch(ch).pulse)
					event_struct(evt).ch(ch).pulse(ps).pulse_data_phe = (event_struct(evt).ch(ch).pulse(ps).pulse_data_mV(2:end)-event_struct(evt).ch(ch).pulse(ps).baseline_mV).*conversion_factor_mV_to_phe_per_10ns(ch);
				end
			end
		end
	end
end




