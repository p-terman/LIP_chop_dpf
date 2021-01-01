function out = LUXExpFcn_framework(param,tt)
% y = LUXExpFcn_framework(param,x)
%
% For LUXExponentialFit_framework
%
% 20110531 - JRV - Created
% 20120124 - PFS - added Gaussian smoothing to leading edge, to simulate the 
%					effect of electronics (etc) shaping. This leads to very good
%					fit quality on alpha S1 events (chi^2/dof). Note that the 
%					smoothing time should NOT be fitted, though it may require adjutment.
% 20120127 - PFS/RJG - modify construction of pulse so that amplitude is decoupled from tau
%					leads to better stability in fit routine


% ensure proper type
param = double(param);
tt = double(tt);

% calculate
ampli = param(1);
tfall = param(2);
offset = param(3);
trise = param(4);

%out = zeros(size(tt));
out = exp(-(tt-offset)./tfall) - exp(-(tt-offset)./trise);
out(tt<offset) = 0;

% The peak time, and peak max
tpeak = -(tfall.*trise).*log(trise./tfall)./(tfall - trise);
peakmax = exp(-tpeak./tfall) - exp(-tpeak./trise);
out = out./peakmax;

if 1 % convolve w Gaussian to create a slower leading edge
    sigma = 10; % ns - smoothing time - tuned by eye -pfs
    dtt = (tt(2)-tt(1))*10; % assume input is in samples, i.e. 1=10 ns, so multiply dx by 10
    A = [1 0 sigma];
    n = ceil(3*sigma/dtt);
    c = gaussian( A , dtt * (-n:n));
    c = c ./ sum(c); % norm area to 1
    out = conv( out , c );
    out = out(n+1:end-n);
end

out = out*ampli;

end