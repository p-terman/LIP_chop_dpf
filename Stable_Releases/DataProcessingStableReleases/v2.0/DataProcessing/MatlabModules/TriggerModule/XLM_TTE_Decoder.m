function [xlm tte] = XLM_TTE_Decoder(event_struct,TTETemplate,EvtTime,options,Debug)
% Developer: Mongkol Moongweluwan
% Version beta
% Released Date: September 3, 2013
%
% Information
% This function is meant to be used in DP framework to translate all
% trigger-related information in .evt into rq.
%
%
% Output: 
% xlm: xlm-related rq. There are, xlm timestamp, trigger number, max filter
% output value and channel, S1 and S2 hit vector.
% tte: tte-related rq. There are, sat value, sat flatness, type, and
% comment.
%
% Version History:
% Version 1.0: Initial released.
%% First XLM business
%Mostly just transtale them to more readable form
if isfield(options,'load_xlm_info')
    if options.load_xlm_info == 1
        try
            xlm.success = 1;
            xlm.xlm_timestamp_samples = event_struct.trigger_timestamp - event_struct.timestamp;
            xlm.trigger_number = event_struct.trigseqnum;
            xlm.max_filter_output_value_mVns = event_struct.max_filter_response*3.82;
            xlm.max_filter_output_channel = 1 + event_struct.max_ch_ID;
            xlm.s1_hit_vector = XLMHitVectorReadout(event_struct.S1_hit_vector);
            xlm.s2_hit_vector = XLMHitVectorReadout(event_struct.S2_hit_vector);
            if Debug.Mode
                Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Successfully calculate xlm rqs.                         ') - Debug.tok;
                pause(Debug.PauseTime)
            end
        catch
            if Debug.Mode
                Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Fail to calculate xlm rqs. xlm entry is corrupted.      ') - Debug.tok;
                pause(Debug.PauseTime)
            end
            xlm.success = 0;
        end
    else
        if Debug.Mode
            Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Fail to calculate xlm rqs. xlm entry is not loaded.         ') - Debug.tok;
            pause(Debug.PauseTime)
        end
        xlm.success = 0;
    end
else
    if Debug.Mode
        Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Fail to calculate xlm rqs. xlm entry is not loaded.             ') - Debug.tok;
        pause(Debug.PauseTime)
    end
    xlm.success = 0;
end
%% Then TTE business
if isfield(options,'load_tte_ch')
if options.load_tte_ch == 1
    tte.success = 1;
    if Debug.Mode
        Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Calculating tte rqs.                                            ') - Debug.tok;
        pause(Debug.PauseTime)
    end
    if length(event_struct.ch) >= 126
        %Use try-catch to catch crashes
        try 
            if Debug.Mode
                Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Analysing tte waveform.                                         ') - Debug.tok;
                pause(Debug.PauseTime)
            end
            %Calculate time separation to flag the case when two or more TTE pulses
            %could be overlapping
            CurrentEvt = find(EvtTime == event_struct.timestamp);
            if CurrentEvt == 1
                TimeLast = 1e6;
                TimeNext = EvtTime(2)-EvtTime(1);
            elseif CurrentEvt == length(EvtTime)
                TimeLast = EvtTime(end)-EvtTime(end - 1);
                TimeNext = 1e6;
            else
                TimeLast = EvtTime(CurrentEvt) - EvtTime(CurrentEvt - 1); 
                TimeNext = EvtTime(CurrentEvt + 1) - EvtTime(CurrentEvt);
            end
    
            %Calculate where to look for the saturation value
            TTEFlatID = find(150 + mean(TTETemplate.Delay([4 7],1)) == event_struct.ch(126).pod_time_samples);
            if ~isempty(TTEFlatID) && ismember(155 + mean(TTETemplate.Delay([4 7],1)), event_struct.ch(126).pod_time_samples)
                tte.tte_sat_value_mV = mean(event_struct.ch(126).pod_data_mV(TTEFlatID:TTEFlatID+5));
                tte.tte_sat_flatness_mV = std(event_struct.ch(126).pod_data_mV(TTEFlatID:TTEFlatID+5));
                %Determine the type based on finding the saturation value
                if ~isempty(find((abs(tte.tte_sat_value_mV - TTETemplate.TTESat(:,1)) < TTETemplate.SatTol(:,1)).*(tte.tte_sat_flatness_mV < TTETemplate.FlatTol(:,1)), 1))
                    switch find((abs(tte.tte_sat_value_mV - TTETemplate.TTESat(:,1)) < TTETemplate.SatTol(:,1)).*(tte.tte_sat_flatness_mV < TTETemplate.FlatTol(:,1)))
                        case 4 %Only DDC
                            tte.tte_type = bin2dec('00001000');
                            tte.tte_comment = NearbyTTECheck(TimeLast,TimeNext);
                        case 5 %Only LED
                            tte.tte_type = bin2dec('00010000');
                            tte.tte_comment = NearbyTTECheck(TimeLast,TimeNext);
                        case 7 %Only Veto
                            tte.tte_type = bin2dec('01000000');
                            tte.tte_comment = NearbyTTECheck(TimeLast,TimeNext);
                        case 9 %DDC and LED
                            tte.tte_type = bin2dec('00011000');
                            tte.tte_comment = NearbyTTECheck(TimeLast,TimeNext);
                        case 10 %DDC and Veto
                            tte.tte_type = bin2dec('01001000');
                            tte.tte_comment = NearbyTTECheck(TimeLast,TimeNext);
                        case 11 %LED and Veto
                            tte.tte_type = bin2dec('01010000');
                            tte.tte_comment = NearbyTTECheck(TimeLast,TimeNext);
                        case 12 %DDC, LED and Veto
                            tte.tte_type = bin2dec('01011000');
                            tte.tte_comment = NearbyTTECheck(TimeLast,TimeNext);
                        otherwise%This suggests the TTE output merged.
                            tte.tte_comment = -1;
                            %'Type is not determined. The TTE pulse may merge with the other nearby TTE.';
                            tte.tte_type = -1;
                        %{
                        %First check how many pulses are merged together.
                        TotPulseMerge = find(TimeNext(evt:end) > 250, 1, 'first');
                        %Start try fitting approach
                        MergeTTE = zeros(350*TotPulseMerge,TotPulseMerge^2);
                        TTEChList = [4 7];
                        switch TotPulseMerge
                            case 1
                                event_struct(evt).TTE_Type = 'Unknown';
                                tte.tte_type = 0;
                                tte.tte_comment = 'Type is not determined. Fitting with template methods fails';
                            case 2
                                for jj = 1:2
                                for kk = 1:2
                                    MergeTTE( TTETemplate.Delay(TTEChList(jj),1)+(1:length(TTETemplate.TTEShape(TTEChList(jj),:))),1) = ...
                                        MergeTTE( TTETemplate.Delay(TTEChList(jj),1)+(1:length(TTETemplate.TTEShape(TTEChList(jj),:))),1) + permute(TTETemplate.TTEShape(TTEChList(jj),:),[2 1]);
                                    MergeTTE( TimeNext(evt)+TTETemplate.Delay(TTEChList(kk),1)+(1:length(TTETemplate.TTEShape(TTEChList(kk),:))),1) = ...
                                        MergeTTE( TimeNext(evt)+TTETemplate.Delay(TTEChList(kk),1)+(1:length(TTETemplate.TTEShape(TTEChList(kk),:))),1) + permute(TTETemplate.TTEShape(TTEChList(kk),:),[2 1]);
                                    % Y = -raw.ch(126).pod(TTEPodID).pod_data_mV((T128(ps) + (0*TTETemplate.Delay(TTEChList(jj),1)) - raw.ch(126).pod(TTEPodID).timestamp + 25):YRangeEnd);
                                    % Idea about getting Y.
                                    % First or all, T128(ps) + TTETemplate.Delay - raw.ch(126).pod(TTEPodID).timestamp + 25 gives 
                                    % the time when the corresponding TTE pod should be above the threshold. Later, we plot
                                    % this relative to T128, hence it is moved
                                    % again by the amount of - TTETemplate.Delay.
                                    YRangeEnd = int32(min(event_struct(evt).ch(126).pod_length_samples(TTEMatchID), + length(MergeTTE) - 24 - event_struct(evt).ch(126).pod_start_samples(TTEMatchID)));
                                    %YRangeEnd = int32(min(raw.ch(126).pod(TTEPodID).length,(T128(ps) - raw.ch(126).pod(TTEPodID).timestamp) + length(MergeTTE)));
                                    Y = event_struct(evt).ch(126).pod_data_mV( (event_struct(evt).ch(126).pod_start_samples(TTEMatchID) + 25):YRangeEnd);
                                    FitRange = 1:length(Y);
                                    X = MergeTTE(24+FitRange,pp);
                                    Resnorm = sum((( X(Y > 0) - Y(Y > 0))).^2)/length(Y(Y > 0));
                                    if Resnorm < 10
                                        event_struct(evt).TTE_Type = TTETemplate.TTEText{TTEChList(jj),1};
                                        event_struct(evt+1).TTE_Type = TTETemplate.TTEText{TTEChList(kk),1};
                                        tte.tte_type = TTEChList(jj);
                                        event_struct(evt+1).TTE_Flag = TTEChList(kk);
                                        tte.tte_comment = 'This pulse merges with the next pulse. Type is determined by fitting with template';
                                        event_struct(evt+1).TTE_Comment = 'This pulse merges with the previous pulse. Type is determined by fitting with template';
                                        evt = evt + 1;
                                    else
                                        event_struct(evt).TTE_Type = 'Unknown';
                                        tte.tte_type = 0;
                                        tte.tte_comment = 'This pulse appears to merge with the next pulse but cannot determine the type.';
                                    end    
                                end
                                end
                            case 3
                                event_struct(evt).TTE_Type = 'Three merged pulses';
                                tte.tte_comment = 'This pulse merges with the next two pulses. Type is not determined';
                                tte.tte_type = 0;
                                event_struct(evt+1).TTE_Type = 'Three merged pulses';
                                event_struct(evt+1).TTE_Comment = 'This pulse merges with the previous and the next pulses. Type is not determined';
                                event_struct(evt+1).TTE_Flag = 0;
                                event_struct(evt+2).TTE_Type = 'Three merged pulses';
                                event_struct(evt+2).TTE_Comment = 'This pulse merges with the previous two pulses. Type is not determined';
                                event_struct(evt+2).TTE_Flag = 0;
                            otherwise
                                event_struct(evt).TTE_Type = 'Unknown';
                                tte.tte_comment = 'Type is not determined. Fitting with template methods fails';
                                tte.tte_type = 0;
                        end
                        %}
                    end
                else
                    tte.tte_sat_value_mV = NaN;
                    tte.tte_sat_flatness_mV = NaN;
                    tte.tte_comment = -2;
                    %'Type is not determined. The TTE pulse is missing. It could be that the veto triggers are too close to generate two spearate gate signal.';
                    tte.tte_type = -2;

                end
            else
                tte.tte_sat_value_mV = NaN;
                tte.tte_sat_flatness_mV = NaN;
                tte.tte_comment = -3;
                %'Can't find TTE pulse that match the trigger'
                tte.tte_type = -3;
            end
        catch exception %Catch in case of crashes
            tte.tte_sat_value_mV = NaN;
            tte.tte_sat_flatness_mV = NaN;
            tte.tte_type = -4;
            tte.tte_comment = -4;
            %'Script crash';
        end
    else
        tte.tte_sat_value_mV = NaN;
        tte.tte_sat_flatness_mV = NaN;
        tte.tte_comment = -5;
        %'Can not find the TTE pulse to analyze. Make sure you load the event_struct with options.load_tte_ch = 1. ';
        tte.tte_type = -5;
    end
else
    if Debug.Mode
        Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Fail to calculate tte rqs. tte is not loaded.           ') - Debug.tok;
        pause(Debug.PauseTime)
    end
    tte.success = 0;
end
else
    if Debug.Mode
        Debug.tok = fprintf(1, '%s\n%s',char(8*ones(1,Debug.tok)),'Fail to calculate tte rqs. tte is not loaded.           ') - Debug.tok;
        pause(Debug.PauseTime)
    end
    tte.success = 0;
end

end

