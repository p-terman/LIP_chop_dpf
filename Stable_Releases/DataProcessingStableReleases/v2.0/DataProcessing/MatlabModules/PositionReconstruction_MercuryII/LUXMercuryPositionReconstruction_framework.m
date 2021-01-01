function [x_cm y_cm] = LUXMercuryPositionReconstruction_framework(signal,LRF)
% This function gives the maximum-likelihood estimate of the position in 
% x,y per event, given a matrix of hit patterns and PMT LRFs
%
% [x y] = LUXMercuryPositionReconstruction_framework(signal,LRF_matrix)
%
% Inputs:
%       signal - A 122 x EVTs matrix of PMT areas
%          LRF - [Optional] A structure containing the PMT's Light Response
%                Functions. Must include the field .ps, which is the 6th
%                order coefficients in the fit exp(polyval(ps,x))
%                If not provided, the code will query the LUG to get the
%                latest LRF parametrization from BG data. Failing that, it
%                will use some approximate LRF data.
% Outputs:
%         x_cm - The event position in cm. 0 is center of PMT 121
%         y_cm - The event position in cm. 0 is center of PMT 121
%
%
%
% Example usage:
%
% Let's say we have some S2 data. For now, just look at the second pulse,
% which is likely to be an S2. Let's pick only 1000 events. Make sure to
% squeeze out the extra dimension!
%
% >> signal = squeeze(d.peak_area_phe(2,1:122,1:1000));
%
% Now run the position reconstruction:
%
% >> [x_cm y_cm] = LUXMercuryPositionReconstruction_framework(signal);
%
% The code will indicate the progress as it processes the data.
%
%
% Once it's done, you can visualize the data. Use:
%
% >> figure; LUXPlotTopPMTs;
%
% for plotting the PMT array. Then:
%
% >> hold on; plot(x_cm,y_cm,'k.','markers',1);
%
%
% Versioning:
%   20120418 CHF - Created
%   20120606 JRV - Made options_merc2 a persistent variable; this should
%                  decrease run time significantly on Oscar based on
%                  previous experience.
%
%   20121216 RJG - Put error trapping for fminunc, so function doesn't
%                   just crash
%
%
%% Initialize

% Record file name for progress/error reporting
fn = mfilename;   % 121216 added by RJG

% TODO: This has to be queried from the LUG IQs. Right now it is saved in a
% .mat file. If LUG query fails, this will be loaded instead.
if ~exist('LRF','var')
    load LRF;
    LRF = LRF_matrix{end};
end

M = size(signal,2);

% PMT position map in cm
load pmt_pos_map

persistent options_merc2
if isempty(options_merc2)
    options_merc2=optimset;
    options_merc2.MaxFunEvals = 100;
    options_merc2.MaxIter = 100;
    options_merc2.TolFun = 0.1;
    options_merc2.TolX = 0.1;
    options_merc2.Display = 'off';
end

warning off

%% Compute starting point very quickly with a modified CG
[x0 y0] = LUXDirtyCG_framework(signal,0.3);

a = tic;

x_cm = zeros(1,M);
y_cm = zeros(1,M);

fprintf('*** Starting *** Performing a maximum-likelihood position reconstruction fit with PMT LRFs\n');

% Get new event positions
for evt = 1:M
    
    signal_evt = signal(:,evt);
    
    % Use modified CG result as initial guess
    params(1) = x0(evt);
    params(2) = y0(evt);
    
    try
    % Function that uses a ML estimator to get best fit
    new_params = fminunc(@LRF_maxlike_fitting_function_framework,params,options_merc2,LRF,signal_evt,pmt_pos_cm);
    catch err % 121216 added by RJG
            fprintf('ERROR TRAP (%s) : fminunc() failed for evt=%d. %s\n',fn,evt,err.identifier);
            new_params = [nan nan];
    end
    
    % Update positions
    x_cm(evt) = new_params(1);
    y_cm(evt) = new_params(2);
    
    %     chisq(evt) = LRF_chisq_function([x_cm(evt) y_cm(evt)],LRF,signal_evt);
    
    if mod(evt,ceil(M/10)) == 0
        b = toc(a);
        fprintf('Done with (%d/%d), elapsed time is %3.0f s. Remaining time = %3.2f min\n',evt,M,b,(M-evt)*b/evt / 60);
    end
    
end

fprintf('*** Finished *** Elapsed time was %3.1f mins. Event process rate was %3.1f Hz\n\n',b/60,M/b)

warning on

