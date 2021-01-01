function Gain = mVns_to_gain(area,Req)
% Gain = mVns_to_gain(xbar,Req)
%
% mVns_to_gain
% This function converts from area (mV*ns units) into Gain
%
% INPUTS:
%         area - value of area to convert
%         Req  - Equivalent termination resistance of circuit in Ohms
%                i.e. 50 Ohm from base in parallel with 50 Ohm from Amplifier
%                Standard value is thus 25 Ohms
%
% OUTPUTS:
%         Gain - Value of area to Gain conversion
%
% cahf 20080814

if nargin < 2
    Req = 25;
end

Gain = area./(Req*1.6e-19*1e12);