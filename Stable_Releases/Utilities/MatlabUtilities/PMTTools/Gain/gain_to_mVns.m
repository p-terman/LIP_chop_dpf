function mVns = gain_to_mVns(Gain,Req)
% mVns = gain_to_mVns(Gain,Req)
%
% gain_to_mVns
% This function converts from Gain into area (mV*ns units)
%
% INPUTS:
%         Gain - Value of area to Gain conversion
%         Req  - Equivalent termination resistance of circuit in Ohms
%                i.e. 50 Ohm from base in parallel with 50 Ohm from Amplifier
%                Standard value is thus 25 Ohms
%
% OUTPUTS:
%         area - value of area to convert
%         
%
% cahf 20080814

if nargin < 2
    Req = 25;
end

mVns = Gain.*(Req*1.6e-19*1e12);