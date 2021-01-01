function [max_area_diff max_long_area max_short_area] = LUXS2FilterMatlab_framework(trace,long_width_samples,short_width_samples)
%
% [max_area_diff max_long_area max_short_area] = LUXS2FilterMatlab_framework(trace,long_width_samples,short_width_samples)
%
% This function will perform two convolutions between the pulse and box
% filters of widths long_width_samples and short_width_samples. 
% The convolutions give the time delay that maximizes the area under the
% filter. Given that time delay, the area under the filters is computed.
%
% The s2filter_max_area_diff is the difference between the area under the
% filters, and can be used for pulse type identifications (S1 vs S2).
%
% Pulses with widths < short_width_samples have s2filter_max_area_diff = 0;
% pulses with widths O(s2filter_max_area_diff) will have a very small
% s2filter_max_area_diff, generally s2filter_max_area_diff < 100.
% Pulses with widths of O(long_width_samples) will generally have a s2filter_max_area_diff
% of 500-1000, or even more.
%
% INPUTS:
%                   trace - the pulse to run the filter on
%        long_width_samples - width of long ("S2") box filter, in samples
%        short_width_samples - width of short ("S1") box filter, in samples
%
%   This S2 filter seems to be optimal for long_width_samples ~ 500, and
%                                          short_width_samples ~  50
%
%   Smaller values will yield large s2filter_max_area_diff for "fat" S1
%   pulses.
%
% OUTPUTS:
%           max_area_diff - the difference between s2filter_max_long_area and s2filter_max_short_area
%             max_long_area - the area under the large ("s2") box filter,
%                           after being delayed to yield the maximum value
%             max_short_area - the area under the small ("s1") box filter,
%                           after being delayed to yield the maximum value
%
% Versioning:
%   20111122 CHF - Created
%   20120826 CHF - Changed 's1' and 's2' in names to 'short' and 'long',
%                  respectively
%   20120826 CHF - Now forcing short box to be inside long box
%                  Added optional plotting for transparent processing
%   20120915 CHF - Changed filter engine to convolution mode. Much faster
%                  and more transparent.
%   20130830 CHF - Fixed 1 sample offset. A delta function input returned
%                  max_area_diff = 1 instead of 0. Problem was related to
%                  incorrect use of 'inrange', which is non-inclusive at
%                  the edges. Also, there was a bug where the actual width
%                  of long and short box filters had an extra box sample added
%                  (e.g. request length 20, had length 21)
%
%% Initialize

M = max([length(trace) long_width_samples*2]);

template_long = zeros(1,M);
template_long(2:(long_width_samples+1)) = 1;

template_short = zeros(1,M);
template_short(2:(short_width_samples+1)) = 1;

tt = 1:length(trace);

%% Run Convolution

[delay_long ftrace_long] = LUXGetFilterOutput(template_long,trace,long_width_samples);

[~,ftrace_short] = LUXGetFilterOutput(template_short,trace,short_width_samples);

% Force short box to be inside long box by zeroing the areas outside this
% range in ftrace_short
% *** Avoid using inrange since it edges are non-inclusive CHF
cut_delay_short_zero = tt >= delay_long & tt <= (delay_long+long_width_samples-short_width_samples);
ftrace_short(~cut_delay_short_zero) = 0;
[~,delay_short] = max(ftrace_short); % replaced maxind here - 20120907 JRV

max_long_area = sum(trace(tt >= delay_long & tt <= delay_long + long_width_samples));
max_short_area = sum(trace(tt >= delay_short & tt <= delay_short + short_width_samples));

max_area_diff = max_long_area - max_short_area;

% Switch to 1 if you want the output plotted - for transparent view
if 0
    figure(432); clf
    plot(tt,trace,'k.-'); hold on
%     plot(tt,ftrace_long,'r--')
    plot((1:M)+delay_long,template_long*max(trace),'ro-');
%     plot(tt,ftrace_short,'b--')
    plot((1:M)+delay_short,template_short*max(trace),'bx-');
    
    l = legend('Data',...%         sprintf('%d \\mus optimal filter output',long_width_samples/100),...
        sprintf('%3.1f \\mus box',long_width_samples/100),...%         sprintf('%d ns optimal filter output',short_width_samples*10),...
        sprintf('%3.1f ns box',short_width_samples*10));
    set(l,'fontsize',12);
    
    xlabel('samples')
    ylabel('phe/sample')
    ylim([-inf max(trace)*1.5])
    ax = axis;
    text(ax(1)+diff(ax(1:2))*0.05,max(trace)*1.3,sprintf('area = %3.1f phe\nmax area diff = %3.1f phe (%2.1f%%)',sum(trace),max_area_diff,max_area_diff./sum(trace)*100),...
        'fontsize',10,'backgroundColor','white','edgeColor','black');
end


function [delay ftrace] = LUXGetFilterOutput(template,trace,width_samples)
% [delay ftraces] = LUXGetFilterOutput(template,traces)
%
% Inputs:
%       template - template for fitting, normalized to 1 in height
%       traces   - data that we'll try to fit template to
%  width_samples - width of the box filter in question
%
% Outputs:
%       delay    - time delay of template
%       ftraces  - filter output
%
% 20111122 CHF - Created based on LUXOptimalFilter.m
%
%% Check inputs

ftrace_temp = conv(trace,template);
startp = width_samples+3;
endp = startp + length(trace) - 1;
ftrace = ftrace_temp(startp:endp);
[~,delay] = max(ftrace);
