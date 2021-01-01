function out = LUXDoubleExpFcn_framework(param,tt)
% out = LUXDoubleExpFcn_framework(param,tt)
%
% Gives the output for a double S1 (each a double exponential!) given fit parameters
% It merely computes a linear combination of two LUXExpFcn_framework
%
% Inputs:
%   param - Concatenated parameters for each exponential fit
%           So param(1:4) is for first pulse and
%              param(5:8) is for second pulse
%      tt - Time vector over which to compute function
%
% Outputs:
%    out - Function output
%
% Versioning:
%   20140827 CHF&AC - Created

out = LUXExpFcn_framework(param(1:4),tt) + LUXExpFcn_framework(param(5:8),tt);