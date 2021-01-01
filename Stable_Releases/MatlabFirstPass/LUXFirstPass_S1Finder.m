% LUXFirstPass_InitVars
%
% DO NOT CALL THIS MANUALLY, IT DOESN'T MAKE SENSE
%
% To be called from LUXFirstPass
% Find S1-like pulses
%
% 2010-02-02 JJC, DCM
%

S1count=0;
S1s=[];
for ii=1:length(datasum)


	ptrace = datasum(ii).data_phe ;
	loopstop=0;
	loopcount=0;
	S1.t0_old = [];
	S1.t2_old = [];
	S1.t1_old = [];
	while loopstop==0
       % keyboard
		peaks = pikpeekser(ptrace,s.ana_settings.thr_pikpeeks); % put traces through peak finder (no filters)
		loopcount=loopcount+1;
		if S2count>=1 % if s2 is only peak in event then it is probably an s1
			not_largest_s2 = ~inrange(datasum(ii).tt(peaks(4)),S2s.t0(1),S2s.t2(1));
		else
			not_largest_s2 = 1;
		end
		if sum(ptrace(peaks(1):peaks(7))>.5) && loopcount<=3 && not_largest_s2 %~inrange(datasum(ii).tt(peaks(4)),S2s.t0(1),S2s.t2(1))% not found too many, and area is high enough (some threshold...) and not in largest S2
		
			LUXFirstPass_AssignS1RQs ;
			
			
		else % if in a peak
			loopstop=1;
		end
		ptrace(peaks(1):peaks(7)) = 0; % zero out pulse to look for more
		S1.t0_old = datasum(ii).tt(peaks(1));
		S1.t2_old = datasum(ii).tt(peaks(7));
		S1.t1_old = datasum(ii).tt(peaks(4));
        %keyboard
	end % while loop
end % for every pulse (ii)

% sort S1 like pulses (the leftovers) by area. do we want to do this, or sort them by something else?
if S1count>0
[S1s.sumarea sortJ] = sort(S1s.sumarea,'descend');
S1s.t0=S1s.t0(sortJ);
S1s.t10l=S1s.t10l(sortJ);
S1s.t50l=S1s.t50l(sortJ);
S1s.t10r=S1s.t10r(sortJ);
S1s.t50r=S1s.t50r(sortJ);
S1s.t1=S1s.t1(sortJ);
S1s.t2=S1s.t2(sortJ);
S1s.t_mean=S1s.t_mean(sortJ);
S1s.t_std=S1s.t_std(sortJ);
S1s.preS1=S1s.preS1(sortJ);
S1s.postS1=S1s.postS1(sortJ);
S1s.prompt_fraction=S1s.prompt_fraction(sortJ);
S1s.peak_area_per_ch=S1s.peak_area_per_ch(sortJ,:);
S1s.head=S1s.head(sortJ,:);
S1s.tail=S1s.tail(sortJ,:);
S1s.peak_height_per_ch=S1s.peak_height_per_ch(sortJ,:);
S1s.peak_height_mV_per_ch=S1s.peak_height_mV_per_ch(sortJ,:);
S1s.baseline_mV=S1s.baseline_mV(sortJ,:);
S1s.sat=S1s.sat(sortJ,:);
S1s.bipolar_disorder=S1s.bipolar_disorder(sortJ);
S1s.whichpulse=S1s.whichpulse(sortJ);
end