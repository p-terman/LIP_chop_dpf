function out = LUXDoubleExpFcn_fixedshape_framework(param,input)
% out = LUXDoubleExpFcn_framework(param,tt)
%
% Gives the output for a double S1 (each a double exponential!) given fit parameters
% It merely computes a linear combination of two LUXExpFcn_framework
%
% Inputs:
%   param - Concatenated parameters for each exponential fit
%           So param(1:4) is for first pulse and
%              param(5:8) is for second pulse
%   input.tt - Time vector over which to compute function
%   input.risetime - risetime parameter (in samples)
%   input.falltime - falltime parameter (in samples)
% Outputs:
%    out - Function output
%
% Versioning:
%   20140827 CHF&AC - Created

tt = input.tt;


param2 = zeros(1,8);
param2(1) = param(1);  % amplitude
param2(2) = input.falltime_samples;  % falltime
param2(3) = param(2);  % offset
param2(4) = input.risetime_samples;  % risetime
param2(5) = param(3);  % amplitude
param2(6) = input.falltime_samples;  % falltime
param2(7) = param(4);  % offset
param2(8) = input.risetime_samples;  % risetime

out = LUXExpFcn_framework(param2(1:4),tt) + LUXExpFcn_framework(param2(5:8),tt);