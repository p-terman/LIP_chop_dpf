% LUXFirstPass_AssignS2RQs
% 
% DO NOT CALL THIS MANUALLY, IT DOESN'T MAKE SENSE
% 
% To be called from LUXFirstPass
% Record S2 RQs to S2 structure
% 
% 2010-02-02 JJC, DCM
% 


S2count=S2count+1;
S2s.whichpulse(S2count) = ii; % which pulse in the summed event this peak belongs to
S2s.t0(S2count) = datasum(ii).tt(times(1));
S2s.ft10l(S2count) = datasum(ii).tt(times(2));
S2s.ft50l(S2count) = datasum(ii).tt(times(3));
S2s.ft1(S2count) = datasum(ii).tt(times(4));
S2s.ft50r(S2count) = datasum(ii).tt(times(5));
S2s.ft10r(S2count) = datasum(ii).tt(times(6));
S2s.t2(S2count) = datasum(ii).tt(times(7));
S2s.diffs(S2count) = diffs(:,1);
S2s.fareas(S2count) = areas(:,1);

S2s.bipolar_disorder(S2count) = 0 ;
%{
                				if inrange(max(datasum(ii).data_phe), .7*abs(min(datasum(ii).data_phe)), 1.3*abs(min(datasum(ii).data_phe))) % 2009-10-06 JJC max between 90% and 110% of min - bipolar!
                					S2s.bipolar_disorder(S2count) = 1 ;
                				end
%}



% generate some other RQs, including the non-filtered times
[S2s.peak_height(S2count) maxind] = max(datasum(ii).data_phe(times(1):times(7))); % peak height
S2s.t1(S2count) = datasum(ii).tt(maxind); % time of max
tt50l = find(datasum(ii).data_phe(times(1):maxind) < (0.5*S2s.peak_height(S2count)),1,'last');
if ~isempty(tt50l)
    S2s.t50l(S2count)=datasum(ii).tt(tt50l);
else
    S2s.t50l(S2count)=datasum(ii).tt(1);
    tt50l=1;
end
tt50r = find(datasum(ii).data_phe(maxind:end) < (0.5*S2s.peak_height(S2count)),1,'first');
if ~isempty(tt50r)
    S2s.t50r(S2count)=datasum(ii).tt(tt50r);
else
    S2s.t50r(S2count)=datasum(ii).tt(end);
    tt50r=length(datasum(ii).tt);
end
tt10l = find(datasum(ii).data_phe(times(1):tt50l) < (0.1*S2s.peak_height(S2count)),1,'last');
if ~isempty(tt10l)
    S2s.t10l(S2count)=datasum(ii).tt(tt10l);
else
    S2s.t10l(S2count)=datasum(ii).tt(1);
end
tt10r = find(datasum(ii).data_phe(tt50r:end) < (0.1*S2s.peak_height(S2count)),1,'first');
if ~isempty(tt10r)
    S2s.t10r(S2count)=datasum(ii).tt(tt10r);
else
    S2s.t10r(S2count)=datasum(ii).tt(end);
end
S2s.t_mean(S2count) = sum(datasum(ii).data_phe(times(1):times(7)).*double(datasum(ii).tt(times(1):times(7)))')/sum(datasum(ii).data_phe(times(1):times(7)));
S2s.t_std(S2count) = sqrt(sum(datasum(ii).data_phe(times(1):times(7)).*double(datasum(ii).tt(times(1):times(7))').*double(datasum(ii).tt(times(1):times(7))'))/sum(datasum(ii).data_phe(times(1):times(7))) - S2s.t_mean(S2count)^2);
                				
% calculate peak_area_per_ch = pulse area by channel

for chi=1:chs
    S2s.peak_area_per_ch(S2count,chi)=0;
    S2s.head(S2count,chi)=0;
    S2s.tail(S2count,chi)=0;
    S2s.peak_height_per_ch(S2count,chi)=0;
    S2s.peak_height_mV_per_ch(S2count,chi)=0;
    S2s.baseline_mV(S2count,chi)=0;
    S2s.sat(S2count,chi)=0;
end
%for ch=1:length(event_struct(nn).ch)
for chi=1:length(datasum(ii).pulse_channels)
    ch=datasum(ii).pulse_channels(chi);
    pn=datasum(ii).pulse_numbers(chi);
    if ~isempty(S2.t2_old) && S2.t2_old<S2s.t1(S2count)
        plow=max(max(S2s.t0(S2count),event_struct(nn).ch(ch).pulse(pn).tt(1)),S2.t2_old); % start sum at either start of peak, or first sample of pulse, or end of last peak that was zeroed out, whichever is latest
    else
        plow=max(S2s.t0(S2count),event_struct(nn).ch(ch).pulse(pn).tt(1)); % start sum at either start of peak, or first sample of pulse, whichever is latest
    end
    if ~isempty(S2.t0_old) && S2.t0_old>S2s.t1(S2count)
        pup=min(min(S2s.t2(S2count),event_struct(nn).ch(ch).pulse(pn).tt(end)),S2.t0_old); % don't run off the pulse to get the end of the peak (could be in a different, longer pulse)
    else
        pup=min(S2s.t2(S2count),event_struct(nn).ch(ch).pulse(pn).tt(end)); % don't run off the pulse to get the end of the peak (could be in a different, longer pulse)
    end
    plow=find(event_struct(nn).ch(ch).pulse(pn).tt==plow);
    pup=find(event_struct(nn).ch(ch).pulse(pn).tt==pup);
    if ~isempty(plow) && ~isempty(pup)
        if pup>plow
            S2s.peak_area_per_ch(S2count,ch) = sum(-event_struct(nn).ch(ch).pulse(pn).pulse_data_phe(plow:pup));
            S2s.head(S2count,ch) = sum(-event_struct(nn).ch(ch).pulse(pn).pulse_data_phe(1:s.daq_settings.sis3301.global.pulse_detect_pretrigger));
            S2s.tail(S2count,ch) = sum(-event_struct(nn).ch(ch).pulse(pn).pulse_data_phe((end-s.daq_settings.sis3301.global.pulse_end_posttrigger+1):end));
            S2s.peak_height_per_ch(S2count,ch) = max(-event_struct(nn).ch(ch).pulse(pn).pulse_data_phe(plow:pup));
            S2s.peak_height_mV_per_ch(S2count,ch) = max(-event_struct(nn).ch(ch).pulse(pn).pulse_data_mV(plow:pup));
        end
    end
    S2s.baseline_mV(S2count,ch) = event_struct(nn).ch(ch).pulse(pn).baseline_mV;
    if any(event_struct(nn).ch(ch).pulse(pn).pulse_data_mV(2:end)>s.daq_settings.sis3301.global.vrange_top*1000) || any(event_struct(nn).ch(ch).pulse(pn).pulse_data_mV(2:end)<s.daq_settings.sis3301.global.vrange_bot*1000)
        S2s.sat(S2count,ch) = 1; % that channel saturated the adc
    else
        S2s.sat(S2count,ch) = 0;
    end
    
    if max(-(event_struct(nn).ch(ch).pulse(pn).pulse_data_mV(2:end)-event_struct(nn).ch(ch).pulse(pn).baseline_mV)) >= s.daq_settings.sis3301.global.pulse_thresh_detect % this pulse went above threshold in this channel
        
        if inrange(max(-event_struct(nn).ch(ch).pulse(pn).pulse_data_phe),0.7*abs(min(-event_struct(nn).ch(ch).pulse(pn).pulse_data_phe)), 1.3*abs(min(-event_struct(nn).ch(ch).pulse(pn).pulse_data_phe)))
            S2s.bipolar_disorder(S2count) = 1 ;
        end
        
    end
end
