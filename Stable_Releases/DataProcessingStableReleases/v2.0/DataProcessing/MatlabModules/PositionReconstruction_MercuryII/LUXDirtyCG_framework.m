function [x_cm y_cm] = LUXDirtyCG_framework(signal,thresh)
% This function implements a very fast modified CG for position
% reconstruction. It's very useful as a quick estimate of position (such as
% for the Mercury algorithm).
%
% [x_cm y_cm] = LUXDirtyCG_framework(signal,thresh)
%
% Inputs:
%       signal - A 122 x EVTs matrix of PMT areas
%       thresh - [Optional. Default 0] PMTs below this threshold (as a
%                fraction of maximum PMT signal) will be zeroed out.
%
% Outputs:
%         x_cm - The event position in cm. 0 is center of PMT 121
%         y_cm - The event position in cm. 0 is center of PMT 121
%
% This function will zero out PMTs below threshold (as a fraction of
% maximum PMT signal), and then will perform a CG on the remaining PMTs.
% This way, the CG output is not too biased towards the center. [Hat tip to
% V. Solovov for the recommendation]
%
% Versioning:
%   20120418 CHF - Created based on D. Malling's personal routine
%
%
%% Do weighted CG pos rec

% load PMT positions
load pmt_pos_map

if ~exist('thresh','var')
    thresh = 0;
end

pmt_pos_cm_top = pmt_pos_cm([1:60 121],:);

%%
% separate into top and bottom signals

top_pmts = [1:60 121];
top_sig = signal(top_pmts,:);

top_sig_norm = top_sig ./ repmat( max(top_sig,[],1), [61 1] );
top_sig_norm(top_sig_norm < thresh) = 0;

x_cm = sum( top_sig_norm .* repmat(pmt_pos_cm_top(:,1),[1 size(top_sig,2)]), 1 ) ./ sum( top_sig_norm, 1 );
y_cm = sum( top_sig_norm .* repmat(pmt_pos_cm_top(:,2),[1 size(top_sig,2)]), 1 ) ./ sum( top_sig_norm, 1 );

% cut out events that don't actually have any signal on the given array (<3 PMTs)
x_cm(sum(top_sig>0,1) < 3) = NaN;
y_cm(sum(top_sig>0,1) < 3) = NaN;


