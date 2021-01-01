% LUXFirstPass_AssignS1RQs
%
% DO NOT CALL THIS MANUALLY, IT DOESN'T MAKE SENSE
%
% To be called from LUXFirstPass
% Record S1 RQs to S1 structure
%
% 2010-02-02 JJC, DCM
%

S1count=S1count+1;
S1s.whichpulse(S1count) = ii; % which pulse in the summed event the peak belongs to
S1s.t0(S1count) = datasum(ii).tt(peaks(1));
S1s.t10l(S1count) = datasum(ii).tt(peaks(2));
S1s.t50l(S1count) = datasum(ii).tt(peaks(3));
S1s.t1(S1count) = datasum(ii).tt(peaks(4));
S1s.t50r(S1count) = datasum(ii).tt(peaks(5));
S1s.t10r(S1count) = datasum(ii).tt(peaks(6));
S1s.t2(S1count) = datasum(ii).tt(peaks(7));
S1s.t_mean(S1count) = sum(datasum(ii).data_phe(peaks(1):peaks(7)).*double(datasum(ii).tt(peaks(1):peaks(7))'))/sum(datasum(ii).data_phe(peaks(1):peaks(7)));
S1s.t_std(S1count) = sqrt(sum(datasum(ii).data_phe(peaks(1):peaks(7)).*double(datasum(ii).tt(peaks(1):peaks(7))').*double(datasum(ii).tt(peaks(1):peaks(7))'))/sum(datasum(ii).data_phe(peaks(1):peaks(7))) - S1s.t_mean(S1count)^2);
S1s.sumarea(S1count) = sum(ptrace(peaks(1):peaks(7)));
S1s.peak_height(S1count) = max(ptrace(peaks(1):peaks(7)));

S1s.bipolar_disorder(S1count) = 0 ;
%{
if inrange(max(datasum(ii).data_phe), .7*abs(min(datasum(ii).data_phe)), 1.3*abs(min(datasum(ii).data_phe))) % 2009-10-06 JJC max between 90% and 110% of min - bipolar!
	S1s.bipolar_disorder(S1count) = 1 ;
end
%}
% pulse-shape checking on S1 added 20090804 pfs, modified for POD 20090915 jjc
%% PROMPT FRACTION
if s.ana_settings.fp_prebins<peaks(2)
    if peaks(2)+s.ana_settings.fp_windowbins>length(ptrace)
        s.ana_settings.fp_windowbins=length(ptrace)-peaks(2);
    end
    S1s.preS1(S1count) = sum(ptrace((peaks(2)-s.ana_settings.fp_prebins):(peaks(2)+s.ana_settings.fp_windowbins)));
else
    S1s.preS1(S1count) = -1;
end
if peaks(2)<(length(ptrace)-s.ana_settings.fp_postbins-1)
    S1s.postS1(S1count) = sum(ptrace(max(peaks(2)-s.ana_settings.fp_prebins,1):peaks(2)+s.ana_settings.fp_postbins));
else
    S1s.postS1(S1count) = -1;
end
S1s.prompt_fraction(S1count) = S1s.preS1(S1count) / S1s.postS1(S1count);
%% /PROMPT FRACTION

for chi=1:chs
    S1s.peak_area_per_ch(S1count,chi)=0;
    S1s.head(S1count,chi)=0;
    S1s.tail(S1count,chi)=0;
    S1s.peak_height_per_ch(S1count,chi)=0;
    S1s.peak_height_mV_per_ch(S1count,chi)=0;
    S1s.baseline_mV(S1count,chi)=0;
    S1s.sat(S1count,chi)=0;
end
for chi=1:length(datasum(ii).pulse_channels)
    ch=datasum(ii).pulse_channels(chi);
    pn=datasum(ii).pulse_numbers(chi);
    if ~isempty(S1.t2_old) && S1.t2_old<S1s.t1(S1count)
        plow=max(max(S1s.t0(S1count),event_struct(nn).ch(ch).pulse(pn).tt(1)),S1.t2_old); % start sum at either start of peak, or first sample of pulse, or end of last peak that was zeroed out, whichever is largest
    else
        plow=max(S1s.t0(S1count),event_struct(nn).ch(ch).pulse(pn).tt(1)); % start sum at either start of peak, or first sample of pulse, whichever is largest
    end
    if ~isempty(S1.t0_old) && S1.t0_old>S1s.t1(S1count)
        pup=min(min(S1s.t2(S1count),event_struct(nn).ch(ch).pulse(pn).tt(end)),S1.t0_old); % don't run off the pulse to get the end of the peak (could be in a different, longer pulse)
    else
        pup=min(S1s.t2(S1count),event_struct(nn).ch(ch).pulse(pn).tt(end)); % don't run off the pulse to get the end of the peak (could be in a different, longer pulse)
    end
    plow=find(event_struct(nn).ch(ch).pulse(pn).tt==plow);
    pup=find(event_struct(nn).ch(ch).pulse(pn).tt==pup);
    if ~isempty(plow) && ~isempty(pup)
        if pup>plow
            S1s.peak_area_per_ch(S1count,ch) = sum(-event_struct(nn).ch(ch).pulse(pn).pulse_data_phe(plow:pup));
            S1s.head(S1count,ch) = sum(-event_struct(nn).ch(ch).pulse(pn).pulse_data_phe(1:s.daq_settings.sis3301.global.pulse_detect_pretrigger));
            S1s.tail(S1count,ch) = sum(-event_struct(nn).ch(ch).pulse(pn).pulse_data_phe((end-s.daq_settings.sis3301.global.pulse_end_posttrigger+1):end));
            S1s.peak_height_per_ch(S1count,ch) = max(-event_struct(nn).ch(ch).pulse(pn).pulse_data_phe(plow:pup));
            S1s.peak_height_mV_per_ch(S1count,ch) = max(-event_struct(nn).ch(ch).pulse(pn).pulse_data_mV(plow:pup));
        end
    end
    S1s.baseline_mV(S1count,ch) = event_struct(nn).ch(ch).pulse(pn).baseline_mV;
    if any(event_struct(nn).ch(ch).pulse(pn).pulse_data_mV(2:end)>200) || any(event_struct(nn).ch(ch).pulse(pn).pulse_data_mV(2:end)<-1800)
        S1s.sat(S1count,ch) = 1; % that channel saturated the adc
    else
        S1s.sat(S1count,ch) = 0;
    end
    if max(-(event_struct(nn).ch(ch).pulse(pn).pulse_data_mV(2:end)-event_struct(nn).ch(ch).pulse(pn).baseline_mV)) >= s.daq_settings.sis3301.global.pulse_thresh_detect % this pulse went above threshold in this channel
        
        if inrange(max(-event_struct(nn).ch(ch).pulse(pn).pulse_data_phe),0.7*abs(min(-event_struct(nn).ch(ch).pulse(pn).pulse_data_phe)), 1.3*abs(min(-event_struct(nn).ch(ch).pulse(pn).pulse_data_phe)))
            S1s.bipolar_disorder(S1count) = 1 ;
        end
        
    end
end