function [aft_box_open aft_left_edge aft_center aft_right_edge aft_box_close preBoxArea postBoxArea] = ...
         LUXPerusePeeksMatlab(pulse_data_phe,BoxWidthSamples,edgeFraction)
%
% [aft_box_open aft_left_edge aft_center aft_right_edge aft_box_close preBoxArea postBoxArea] = ...
%         LUXPerusePeeksMatlab(pulse_data_phe,BoxWidthSamples,edgeFraction)
%
% This function is based on the algorithm used in mex function perusePeeks.c
% This is a Matlab-only implementation that should be easily debugged and,
% more importantly, *not crash*.
%
% Inputs:
%   pulse_data_phe - the data in phe/sample
%  BoxWidthSamples - the width of the box, in samples
%     edgeFraction - this defines where the pulse edges are, i.e. at what
%                    fraction of total area do you consider the pulse edges
%                    to be. Examples are 0.01, 0.05, etc. [OPTIONAL, default 0.01]
%
% Outputs:
%   Note that 'aft' means 'area fractonal timing', or the timing where an
%   area fraction in the fractional cumulative sum is met.
%
%   aft_box_open - what sample does the box (with length
%                          BoxWidthSamples) start
%  aft_left_edge - the left edge of the pulse, defined as the first sample
%                          BEFORE fractional cumsum becomes > edgeFraction
%     aft_center - the 0.5 fractional area sample
% aft_right_edge - the right edge of the pulse, defined as the first sample
%                          AFTER fractional cumsum becomes > (1-edgeFraction)
%  aft_box_close - This is just aft_box_open + BoxWidthSamples
%             preBoxArea - The area outside the box on the left
%            postBoxArea - The area outside of the box on the right
%
% If the box (defined by BoxWidthSamples) is longer than the pulse, then 
%  aft_box_open = 1, and
%  aft_box_close = length(pulse_data_phe);
%
% Versioning
%
% 20121111 CHF - Created. Based on Sorensen's perusePeeks.c
% 20130304 JRV - Added a kludge around line 100 to fix index out of bounds
%                problem - CHF or PFS should verify
%


N = length(pulse_data_phe);

if ~exist('edgeFraction','var')
    edgeFraction = 0.01;
end

if N < BoxWidthSamples
    
    % If box completely engulfs pulse, then no need for box convolution
    % Just get the cumsum, edges are beginning and end of pulse
    cumsumtrace_inbox = cumsum(pulse_data_phe);

    aft_box_open = 1;
    aft_box_close = N;
    
    % And there are no wing areas
    preBoxArea = 0;
    postBoxArea = 0;
    
    delay = 0;
    
else
    
    % Do a box convolution
    [delay convtrace] = boxconv(pulse_data_phe,BoxWidthSamples);

    t = 1:length(pulse_data_phe);
    
    aft_box_open = delay;
    aft_box_close = delay + BoxWidthSamples;
    
    cut = inrange(t,aft_box_open,aft_box_close+0.1);

    cuttrace = pulse_data_phe(cut);
    cumsumtrace_inbox = cumsum(cuttrace);

    preBoxArea = sum(pulse_data_phe(1:aft_box_open));
    postBoxArea = sum(pulse_data_phe(aft_box_close:end));
    
end

if aft_box_close > length(pulse_data_phe);
    aft_box_close = length(pulse_data_phe);
end

cumsumtrace_inbox = cumsumtrace_inbox - cumsumtrace_inbox(1);
cumsumtrace_inbox_norm = cumsumtrace_inbox./cumsumtrace_inbox(end);

% Get left edge by doing a find backwards - this guarantees it finds the
% first sample BEFORE cumsum goes to edgeFraction
aft_left_edge = length(cumsumtrace_inbox_norm)-find(cumsumtrace_inbox(end:-1:1) <= edgeFraction,1,'first') + delay + 1;

% Center is just cumsum 0.5. I'm doing this on the right.
aft_center = find(cumsumtrace_inbox_norm >= 0.5,1,'first') + delay;

% Get right edge normally
aft_right_edge = find(cumsumtrace_inbox_norm >= (1-edgeFraction),1,'first') + delay;

% added this to fix crashes due to defining the right edge outside of the
% pulse bounds - 2013-03-04 - JRV
if aft_right_edge > length(pulse_data_phe)
    aft_right_edge = length(pulse_data_phe);
end



%% Just for checking

% if 0
% 
% figure(432); clf
% 
% norm = max(pulse_data_phe)./max(convtrace);
% 
% markloc = -0.1;
% 
% plot(t,pulse_data_phe,'k-'); hold on
% plot((1:sum(cut))+delay-1,cumsumtrace_inbox.*norm,'b--')
% plot(aft_left_edge,markloc,'b+','markers',15)
% plot(aft_center,markloc,'b+','markers',15)
% plot(aft_right_edge,markloc,'b+','markers',15)
% 
% ylim([-0.2 max(pulse_data_phe)*2])
% xlim([delay-50,delay+BoxWidthSamples+50])
% ax = axis;
% 
% plot((1:sum(cut))+delay-1,pulse_data_phe(cut),'r.-')
% plot(delay*[1 1],ax(3:4),'r-')
% plot((delay+BoxWidthSamples)*[1 1],ax(3:4),'r-')
% 
% xlabel('samples')
% 
% end

