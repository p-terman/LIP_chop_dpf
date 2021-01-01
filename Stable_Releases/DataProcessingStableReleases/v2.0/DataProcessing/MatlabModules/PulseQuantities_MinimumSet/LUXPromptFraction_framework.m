function fp = LUXPromptFraction_framework(trace,prebins,twindow,t10l)
%fp = LUXPromptFraction_framework(trace,prebins,twindow,t10l)
% fp is the prompt fraction
% trace is the pulse you wish to calculate (sum from all pmts!)
% prebins and twindow are the parameters for the calculation. prebins
% defaults to 5, twindow defaults to 4.
% t10l is the time for 10% rise of the pulse - when you start the
% integration
%
% 2011-03-22 JJC - created for new first pass (v6)

if ~exist('prebins','var')
    prebins = 5;
end

if ~exist('twindow','var')
    twindow = 4;
end

fp = 0;

num = sum(trace(max((t10l-prebins),1):min((t10l+twindow),length(trace))));
den = sum(trace(1:end));

fp = num/den;