% LUXFirstPass_InitVars
%
% DO NOT CALL THIS MANUALLY, IT DOESN'T MAKE SENSE
%
% To be called from LUXFirstPass
% Find all S2-like pulses
%
% 2010-01-27 JJC, DCM
%


S2count=0;
S2s=[];

% Pass each pulse through the pikS2peek function (find S2s-ish. really just look for things that i think are S2s...) 
for ii=1:length(datasum) % for every pulse in the event window

	ptrace=zeros(1,s.ana_settings.s2window*2+length(datasum(ii).data_phe)); % put some leading and trailing zeros to deal with filters
	ptrace(s.ana_settings.s2window+1:(length(datasum(ii).data_phe)+s.ana_settings.s2window)) = datasum(ii).data_phe ;
                		
	loopstop=0; % if 1, stop the s2 finding loop. there are no more left, or i've found too many
	loopcount=0; % counter
	S2.t0_old = []; % times of previously found s2 so that pulses found do not overlap
	S2.t2_old = [];
	S2.t1_old = [];
	while loopstop==0 % while keep looking for s2s
		[times areas diffs] = pikS2peeks(ptrace,s.ana_settings.s2window,s.ana_settings.s1window,[],[],s.ana_settings.x50,1); % put trace through S2 finder for 1 pulse at a time
        if flag_debug; fprintf('s2window = %d\ts1window = %d\n',s.ana_settings.s2window,s.ana_settings.s1window); end
		loopcount=loopcount+1;
		times=times(:,1)-s.ana_settings.s2window; % account for filter non-kludge above(take out leading s.ana_settings.s2window zeros)
		times(find(times<1))=1; % pulse can't start before beginning of pulse - duh.
		times(find(times>length(datasum(ii).tt)))=length(datasum(ii).tt); % pulse can't end after end of pulse - duh.
		
        if diffs(1)>s.ana_settings.s2threshold % if there is an S2
            % should this really be zero, or should it be some threshold
            % based on the noise? 2010-02-03 JJC
            %keyboard
			LUXFirstPass_AssignS2RQs ; % script to calculate all RQs for this S2
	
		
		else
			loopstop=1; % there are no more S2s in this pulse
		end % if diffs>0 there was an s2
		
		if loopcount>s.ana_settings.max_nb_s2s % stop looking, i've found enough
			loopstop=1;
		end
		ptrace((times(1)+s.ana_settings.s2window):(times(7)+s.ana_settings.s2window))=0; % zero out pulse to look for another S2
		S2.t0_old = datasum(ii).tt(times(1));
		S2.t2_old = datasum(ii).tt(times(7));
		[peak_height_temp maxind_temp] = max(datasum(ii).data_phe(times(1):times(7))); % peak height
		S2.t1_old = datasum(ii).tt(maxind_temp); % time of max
	end % while loop
end % for every pulse (ii)
% now sort these by the "S2-ness" from diffs
if S2count>0
	[S2s.diffs sortI] = sort(S2s.diffs,'descend');
	S2s.t0=S2s.t0(sortI);
	S2s.ft10l=S2s.ft10l(sortI);
	S2s.ft50l=S2s.ft50l(sortI);
	S2s.ft1=S2s.ft1(sortI);
	S2s.t1=S2s.t1(sortI);
	S2s.t2=S2s.t2(sortI);
	S2s.ft50r=S2s.ft50r(sortI);
	S2s.ft10r=S2s.ft10r(sortI);
	S2s.ftareas=S2s.fareas(sortI);
	S2s.t10l=S2s.t10l(sortI);
	S2s.t50l=S2s.t50l(sortI);
	S2s.t50r=S2s.t50r(sortI);
	S2s.t10r=S2s.t10r(sortI);
	S2s.t_mean=S2s.t_mean(sortI);
	S2s.t_std=S2s.t_std(sortI);
	S2s.peak_area_per_ch=S2s.peak_area_per_ch(sortI,:);
	S2s.head=S2s.head(sortI,:);
	S2s.tail=S2s.tail(sortI,:);
	S2s.peak_height_per_ch=S2s.peak_height_per_ch(sortI,:);
	S2s.peak_height_mV_per_ch=S2s.peak_height_mV_per_ch(sortI,:);
	S2s.baseline_mV=S2s.baseline_mV(sortI,:);
	S2s.sat=S2s.sat(sortI,:);
	S2s.bipolar_disorder=S2s.bipolar_disorder(sortI);
	S2s.whichpulse=S2s.whichpulse(sortI);
end