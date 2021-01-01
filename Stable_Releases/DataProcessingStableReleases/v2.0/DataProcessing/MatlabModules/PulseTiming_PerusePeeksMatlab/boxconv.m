function [delay convtrace] = boxconv(waveform,box_width)
% [delay convtrace] = boxconv(waveform,box_width)
%
% This function performs a waveform convolution with a box of width box_width.
% The output of the convolution (truncated to give the same length as
% original waveform) is given as 'convtrace', and the maximum (delay for box
% enclosing) is given by 'delay'.
%
% Inputs:
%     waveform - the input waveform to convolve
% box_width - the width of the box for convolution
%
% Outputs:
%     delay - the location of the max in convtrace
% convtrace - the convolution output between waveform and a box of width
%             box_width
%
% Versioning:
%   20121111 CHF - Created
%   20130204 CHF - Added error check for screwed up baseline. This causes
%                  the area to always be negative. If delay is < 1, then
%                  just return 1.
%   20130418 CHF - SIGNIFICANT minor bug, thanks Scott, asshole.
%%

% Pad the waveform with zeros to ensure box can slide through it
waveform_long = zeros(1,2*box_width + length(waveform));
waveform_long(box_width:(length(waveform)+box_width-1)) = waveform;

% Make box waveform
box_waveform = zeros(1,length(waveform_long));
box_waveform(2:(box_width+2)) = 1;

% Perform convolution
convtrace_temp = conv(waveform_long,box_waveform);

% Truncate output to give same length as original waveform
startp = box_width;
endp = startp + length(waveform) - 1;
convtrace = convtrace_temp(startp:endp);
[~,delay_temp] = max(convtrace_temp);

% Delay must also be shifted back due to truncation
delay = delay_temp - 2*box_width + 1;

if delay < 1
    delay = 1;
end

