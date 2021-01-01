
% This one simply check how confident we are in reading out the TTE.
function TextCode = NearbyTTECheck(TimeLast,TimeNext)
    if TimeLast > 350 && TimeNext > 350
        TextCode = 1;
        %'Determined with saturated value. No other nearby TTE pulse within 3.5 us';
    elseif TimeLast > 350 && TimeNext < 350
        TextCode = 2;
        %'Determined with saturated value. The next trigger pulse is within 3.5 us';
    elseif TimeLast < 350 && TimeNext > 350
        TextCode = 3;
        %'Determined with saturated value. The previous trigger pulse is within 3.5 us';
    else
        TextCode = 4;
        %'Determined with saturated value. The previous and the next trigger pulses are within 3.5 us';
    end
end