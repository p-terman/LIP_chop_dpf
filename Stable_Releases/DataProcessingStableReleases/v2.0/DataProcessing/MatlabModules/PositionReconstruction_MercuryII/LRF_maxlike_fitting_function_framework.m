function loglike = LRF_maxlike_fitting_function_framework(xy,LRF,signal_evt,pmt_pos_cm)
% This function needs to be used with LUXMercuryPositionReconstruction_framework. It
% will output the log-likelihood value of the poisson-fluctuated hit
% pattern comparison to the LRF functions, given an estimated event
% position
%
% loglike = LRF_maxlike_fitting_function_framework(params,LRF,signal_evt,pmt_pos_cm)
%
% Inputs:
%          xy - [x y] for a given event.
%         LRF - A structure containing the PMT's Light Response
%               Functions. Must include the field .params, which has the
%               fit parameters
%  signal_evt - 122 x 1 vector. The hit pattern for this event.
%  pmt_pos_cm - The pmt centers, in cm, as loaded from pmt_pos_map.mat
%
% Outputs:
%     loglike - The log-likelihood
%
%
% Versioning:
%   20120418 CHF - Created
%   20130106 CHF - Changed how LRF fitting is done - now using 2-parameter
%                  analytical form (instead of 7 coefficient polynomial).
%
%
%%
top = [1:60 121];

dis_r = sqrt( (pmt_pos_cm(:,1) - xy(1)).^2 + (pmt_pos_cm(:,2) - xy(2)).^2 );

params = LRF.params;

QE = 0.517; % make sure this is correct
dis_LRF_vals = params(top,1).*((1 - (1 - params(top,2) - QE).^(1./dis_r(top))) ./ (params(top,2) + QE)).^3;

% dis_LRF_vals = exp( ps(:,1).*dis_r(:).^4 + ps(:,2).*dis_r(:).^3 + ps(:,3).*dis_r(:).^2 + ps(:,4).*dis_r(:) + ps(:,5) );
% dis_LRF_vals = exp( ps(top,1).*dis_r(top).^6 + ps(top,2).*dis_r(top).^5 + ps(top,3).*dis_r(top).^4 +...
%                     ps(top,4).*dis_r(top).^3 + ps(top,5).*dis_r(top).^2 + ps(top,6).*dis_r(top) + ps(top,7) );

N_photons = sum( signal_evt(top)) ./ sum(dis_LRF_vals);

loglike =  -sum( signal_evt(top) .* log(N_photons .* dis_LRF_vals) - (N_photons .* dis_LRF_vals) );

if ~isfinite(loglike)
    loglike = -9999999999;
end


