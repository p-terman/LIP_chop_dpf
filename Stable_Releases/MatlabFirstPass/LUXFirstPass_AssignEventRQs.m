% LUXFirstPass_AssignEventRQs
%
% DO NOT CALL THIS MANUALLY, IT DOESN'T MAKE SENSE
%
% To be called from LUXFirstPass
% Record RQs for event
%
% 2010-02-02 JJC, DCM
%

numS2 = min(S2count,s.ana_settings.max_nb_s2s); % limit of s2-like pulses is the total found, or the number requested, whichever is smaller

if S2count>0
    whichpeak(1:numS2,nnn) = 2; % 2 for S2. that makes sense...
    for ii=1:numS2
        if S2s.bipolar_disorder(ii)==1
            whichpeak(ii,nnn) = 3 ; % sorry, not an S2, it's bipolar, which is awesome!
        end
    end
    
    whichpulse(1:numS2,nnn) = S2s.whichpulse(1:numS2);
    t0(1:numS2,nnn) = S2s.t0(1:numS2);
    t10l(1:numS2,nnn) = S2s.t10l(1:numS2);
    t10r(1:numS2,nnn) = S2s.t10r(1:numS2);
    t50l(1:numS2,nnn) = S2s.t50l(1:numS2);
    t50r(1:numS2,nnn) = S2s.t50r(1:numS2);
    t1(1:numS2,nnn) = S2s.t1(1:numS2);
    t2(1:numS2,nnn) = S2s.t2(1:numS2);
    ft10l(1:numS2,nnn) = S2s.ft10l(1:numS2);
    ft10r(1:numS2,nnn) = S2s.ft10r(1:numS2);
    ft50l(1:numS2,nnn) = S2s.ft50l(1:numS2);
    ft50r(1:numS2,nnn) = S2s.ft50r(1:numS2);
    ft1(1:numS2,nnn) = S2s.ft1(1:numS2);
    peak_height(1:numS2,nnn) = S2s.peak_height(1:numS2);
    s2candidate_areadiffs(1:numS2,nnn) = S2s.diffs(1:numS2);
    t_mean(1:numS2,nnn) = S2s.t_mean(1:numS2);
    t_std(1:numS2,nnn) = S2s.t_std(1:numS2);
    peak_area_per_ch(1:numS2,:,nnn) = S2s.peak_area_per_ch(1:numS2,:);
    peak_height_per_ch(1:numS2,:,nnn) = S2s.peak_height_per_ch(1:numS2,:);
    peak_height_mV_per_ch(1:numS2,:,nnn) = S2s.peak_height_mV_per_ch(1:numS2,:);
    baseline_mV(1:numS2,:,nnn) = S2s.baseline_mV(1:numS2,:);
    head(1:numS2,:,nnn) = S2s.head(1:numS2,:);
    tail(1:numS2,:,nnn) = S2s.tail(1:numS2,:);
    sat(:,nnn) = any(S2s.sat(1:end,:)',2);
end %if S2count>0

if S1count>0
    numS1 = min((max_nb_pulses-numS2),S1count); % there may not be enough pulses in the event, but that's okay.
    whichpeak((numS2+1):(numS2+numS1),nnn) = 1; % 1 for S1...
    for ii=1:numS1
        if S1s.bipolar_disorder(ii)==1
            whichpeak(numS2+ii,nnn) = 3 ; % sorry, not an S1, it's bipolar. i'm so sad...
        end
    end
    
    whichpulse((numS2+1):(numS2+numS1),nnn) = S1s.whichpulse(1:numS1);
    t0((numS2+1):(numS2+numS1),nnn) = S1s.t0(1:numS1);
    t10l((numS2+1):(numS2+numS1),nnn) = S1s.t10l(1:numS1);
    t10r((numS2+1):(numS2+numS1),nnn) = S1s.t10r(1:numS1);
    t50l((numS2+1):(numS2+numS1),nnn) = S1s.t50l(1:numS1);
    t50r((numS2+1):(numS2+numS1),nnn) = S1s.t50r(1:numS1);
    t1((numS2+1):(numS2+numS1),nnn) = S1s.t1(1:numS1);
    t2((numS2+1):(numS2+numS1),nnn) = S1s.t2(1:numS1);
    peak_height((numS2+1):(numS2+numS1),nnn) = S1s.peak_height(1:numS1);
    t_mean((numS2+1):(numS2+numS1),nnn) = S1s.t_mean(1:numS1);
    t_std((numS2+1):(numS2+numS1),nnn) = S1s.t_std(1:numS1);
    peak_area_per_ch((numS2+1):(numS2+numS1),:,nnn) = S1s.peak_area_per_ch(1:numS1,:);
    peak_height_per_ch((numS2+1):(numS2+numS1),:,nnn) = S1s.peak_height_per_ch(1:numS1,:);
    peak_height_mV_per_ch((numS2+1):(numS2+numS1),:,nnn) = S1s.peak_height_mV_per_ch(1:numS1,:);
    baseline_mV((numS2+1):(numS2+numS1),:,nnn) = S1s.baseline_mV(1:numS1,:);
    head((numS2+1):(numS2+numS1),:,nnn) = S1s.head(1:numS1,:);
    tail((numS2+1):(numS2+numS1),:,nnn) = S1s.tail(1:numS1,:);
    preS1((numS2+1):(numS2+numS1),nnn) = S1s.preS1(1:numS1);
    postS1((numS2+1):(numS2+numS1),nnn) = S1s.postS1(1:numS1);
    prompt_fraction((numS2+1):(numS2+numS1),nnn) = S1s.prompt_fraction(1:numS1);
    sat(:,nnn) = any(S1s.sat(1:end,:)',2) | sat(:,nnn); % any pulse in the event window saturates
    evt_empty(nnn)=0;
end

