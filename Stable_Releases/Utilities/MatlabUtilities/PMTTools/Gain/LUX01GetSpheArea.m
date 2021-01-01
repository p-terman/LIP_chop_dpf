% sphe_area = LUX01GetSpheArea(ch,bias)
%
% ONLY INTENDED FOR ESTIMATION PURPOSES
% ...when an actual calibration (especially at low gains, below noise level) 
%    is not available
%
% Given a bias and a LUX01 Ch. # (set for Run008-009 configuration), the
% function outputs the estimated sphe area based on a linear fit in loglog space
% gain vs. bias from known calibrations
%
% INPUTS:
%          ch - channel number
%        bias - in Volts
%
% OUTPUTS:
%   sphe_area - in mVns (on a 25 Ohm termination)
%
% 20090709 - CHF

function sphe_area = LUX01GetSpheArea(ch,bias)

%Cold data
% data.bias1 = 1000:100:1200;
% data.bias2 = 1000:100:1200;
% data.bias3 = 1000:100:1200;
% data.bias4 = 1000:100:1200;

% 20091001, Run010
data.bias1 = 1100:100:1400;
data.bias2 = 1100:100:1400;
data.bias3 = 1100:100:1400;
data.bias4 = 1100:100:1400;

data.ch1 = [24.99 46.52 82.83 151.13 ];
data.ch2 = [15.14 28.68 55.14 95.49  ];
data.ch3 = [12.19 20.71 37.29 66.85  ];
data.ch4 = [22.37 42.86 79.51 149.45 ];

data.err1 = ones(1,numel(data.bias1));
data.err2 = ones(1,numel(data.bias2));
data.err3 = ones(1,numel(data.bias3));
data.err4 = ones(1,numel(data.bias4));

[a b] = linfit(log10(data.(['bias' num2str(ch)])),log10(data.(['ch' num2str(ch)])),data.(['err' num2str(ch)]));

sphe_area = 10.^(a + b*log10(bias));
