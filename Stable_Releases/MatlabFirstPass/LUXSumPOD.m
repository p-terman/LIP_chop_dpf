function [event_struct] = LUXSumPOD(event_struct)
% function [event_struct] = LUXSumPOD(event_struct)
% 
% Inputs:		- event_struct (output from LUXEventLoader)
%               
% Output:		- event_struct (same as input, only new channel added 
%                  containing the sum) 
%                  Finds all overlapping pulses and calculates a summed 
%                  pulse to be used by REEFrunner
% 
% v1.0 - JJC created 2009-09-03
% v1.1 - JJC, DCM - changed inputs from settings structure to values
% v1.2 - JJC, DCM - expects that LUXCalibratePulses has been run and .data_phe field exists.
% v1.3 - JJC - now creates sumpp as well as sumtt, to label the summed
% pulse that tt belongs to.

%conversion_factor_mV_to_phe_per_10ns = 10./(mVns_per_phe .* amp_gain) ;

if length(event_struct)<2 % catch for sending in raw data output
	event_struct(1) = event_struct;
	event_struct(1).empty=0;
end



% loop over all events

%event_struct = event_struct(1); % for debugging/development
tic
for evt=1:length(event_struct)
if ~event_struct(evt).empty
	event_struct(evt).sumtt = []; % initialize sum tt
	for ch=1:min(4,length(event_struct(evt).ch))
	if ~isempty(event_struct(evt).ch(ch).pulse)
		for ps=1:length(event_struct(evt).ch(ch).pulse)
			%event_struct(evt).ch(ch).pulse(ps).pulse_data_phe = (event_struct(evt).ch(ch).pulse(ps).pulse_data_mV(2:end)-event_struct(evt).ch(ch).pulse(ps).baseline_mV).*conversion_factor_mV_to_phe_per_10ns(ch); % kludge! - fix first sample problem
			event_struct(evt).ch(ch).pulse(ps).length = event_struct(evt).ch(ch).pulse(ps).length-1; % kludge!
			if isfield(event_struct(evt).ch(ch).pulse(ps),'start')
				event_struct(evt).ch(ch).starts(ps) = event_struct(evt).ch(ch).pulse(ps).start ;
				event_struct(evt).ch(ch).pulse(ps).tt = (event_struct(evt).ch(ch).pulse(ps).start):1:((event_struct(evt).ch(ch).pulse(ps).start)+int32(event_struct(evt).ch(ch).pulse(ps).length)-1);
			else
				event_struct(evt).ch(ch).starts(ps) = event_struct(evt).ch(ch).pulse(ps).timestamp ; % if raw data, then start is timestamp
				event_struct(evt).ch(ch).pulse(ps).tt = (event_struct(evt).ch(ch).pulse(ps).timestamp):1:((event_struct(evt).ch(ch).pulse(ps).timestamp)+int32(event_struct(evt).ch(ch).pulse(ps).length)-1);
			end
			%event_struct(evt).ch(ch).pulse(ps).tt = (event_struct(evt).ch(ch).pulse(ps).start):1:((event_struct(evt).ch(ch).pulse(ps).start)+int32(event_struct(evt).ch(ch).pulse(ps).length)-1);
			%event_struct(evt).ch(ch).ends(ps) = event_struct(evt).ch(ch).pulse(ps).tt(end) ;
			event_struct(evt).sumtt = union(event_struct(evt).sumtt, event_struct(evt).ch(ch).pulse(ps).tt);
		end
	end
    end
        event_struct(evt).sumpp = zeros(1,length(event_struct(evt).sumtt));
        
		breaks = find(diff(event_struct(evt).sumtt)>1);
		breaks_ext = [];
		breaks_ext(1:length(breaks)) = breaks;
		breaks_ext(end+1) = length(event_struct(evt).sumtt);
		
		for ii=1:length(breaks_ext)
		
			if ii==1
				event_struct(evt).chsum.pulse(ii).tt = event_struct(evt).sumtt(1:breaks_ext(ii));
                event_struct(evt).sumpp(1:breaks_ext(ii)) = ones(1,length(event_struct(evt).chsum.pulse(ii).tt));
			else
				event_struct(evt).chsum.pulse(ii).tt = event_struct(evt).sumtt((breaks_ext(ii-1)+1):breaks_ext(ii));
                event_struct(evt).sumpp((breaks_ext(ii-1)+1):breaks_ext(ii)) = ii*ones(1,length(event_struct(evt).chsum.pulse(ii).tt));
			end
			event_struct(evt).chsum.pulse(ii).data_phe = zeros(length(event_struct(evt).chsum.pulse(ii).tt),1);
		end
	%end
	%end

	for sps=1:length(event_struct(evt).chsum.pulse)
		ii=0;
		for ch=1:min(4,length(event_struct(evt).ch))
		
			[starts_intersects spsi psi] = intersect(event_struct(evt).chsum.pulse(sps).tt,event_struct(evt).ch(ch).starts);
			if ~isempty(starts_intersects)
			for psii=psi % for every pulse in the channel that belongs in this summed pulse (sps)
				ii = ii+1;
				%event_struct(evt).chsum.pulse(sps).ch(ch).pulse_numbers(ii)=psii; % record which pulses in each channel belong to this summed pulse.
				event_struct(evt).chsum.pulse(sps).pulse_numbers(ii)=psii; % two vectors with pulse numbers and channels to map
				event_struct(evt).chsum.pulse(sps).pulse_channels(ii)=ch;
				[intersects sps_ind psii_ind] = intersect(event_struct(evt).chsum.pulse(sps).tt,event_struct(evt).ch(ch).pulse(psii).tt);
				event_struct(evt).chsum.pulse(sps).data_phe(sps_ind) = (event_struct(evt).chsum.pulse(sps).data_phe(sps_ind)) ... 
				+(event_struct(evt).ch(ch).pulse(psii).pulse_data_phe) ;
			end
			end
        end
         
	end

	%{
	loop through channels
		loop through pulses
			loop through remaining channels (and sum channel)
				if there is intersection, mark channel and pulse number to be added later
	%}

	%for ch=1:lengh(event_struct(evt))
	
	%	for ps=1:length(event_struct(evt).ch(ch))
		
	%		for chch=ch:length(event_struct(evt))
			
	%			for psps=1:length(event_struct(evt).ch(chch))
				
	%				[intsecs psi pspsi] = intersect(event_struct(evt).ch(ch).pulse(ps).tt, event_struct(evt).ch(chch).pulse(psps).tt);
	%				if ~isempty(intsecs) % if pulses overlap
					
	%					event_struct(evt).chsum.pulse(sps).tt = union(event_struct(evt).ch(ch).pulse(ps).tt, event_struct(evt).ch(chch).pulse(psps).tt);
					
	%				end
				
	%			end
			
	%		end
			
			% do loop for sum channel as well
		
	%	end
	
	%end
end
end
toc

