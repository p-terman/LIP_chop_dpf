function [param chisq dof] = LUXExponentialFit_framework(xdata, ydata, tol_settings, cut)
% This function fits an exponential (S1) to pulses in the first pass.
% [param chisq dof] = LUXExponentialFit_framework(xdata, ydata, cut)
%
% Inputs:
% 
%          xdata - x data (samples)
%          ydata - data trace
%            cut - logical mask for points to include in fit (all points
%                  are still included in chisq)
%   tol_settings - tolerance settings structure (optional)
%
% Outputs:
%   
%          param - amplitude, tau fall, location of max, tau rise
%          chisq - chi squared for fit
%            dof - degrees of freedom for fit

%
% 20110323 - JRV - Created
% 20110426 - JRV - Added do_plot option
% 20110601 - JRV - Now calls external function LUXExpFcn
% 20110907 - JRV - Now uses global vars to avoid multiple calls to optimset
% 20111029 - CHF - Output is now separating chisq and dof, to match new RQ1 variables
% 20111108 - PFS - Changed initial parameter for tau_rise
% 20111009 - JRV - Added cut input option
% 20120124 - PFS - subtract 3 samples from initial guess for param(3)
%					(offset). Appears to give more consistently reliable fits. 
%					Also, lower initial guess for rise time to 1.5 ns - this removes
%					a satellite local minimum in the tau populations
% 20120127 - PFS/RJG - IMPORTANT. added bounds on allowed parameter space in fit routine. 
%					Input data must now be POSITIVE, i.e. raw LUX data should be input as -data.
%					Fit stability seems improved. 
% 20130512 - JRV - Reduced tolerance values to 1e-2
%                  Set maximum iters to 100
% 20130529 - JRV - Now accepts tolerance settings

%% check inputs

if nargin < 2
    fprintf('LUXExponentialFit_framework requires at least two arguements\n');
    return;
end

if nargin < 4
    cut = true(1,length(xdata));
end

%% ensure that data is proper type

xdata = double(xdata);
ydata = double(ydata);

if sum(size(xdata) ~= size(ydata))
    ydata = ydata';
end

%% assign initial parameter estimates

% amplitude and mean

[ymin ymin_loc] = min(ydata);
[ymax ymax_loc] = max(ydata);


if abs(ymax) > abs(ymin)
    init_params(1) = ymax;
    init_params(3) = xdata(ymax_loc)-3;
else
    init_params(1) = ymin;
    init_params(3) = xdata(ymin_loc)-3;
end
% time constants

init_params(2) = 3; % ~29 ns
init_params(4) = 1.5; % this appears to be a better initial guess -pfs

%% use matlab non-linear least squares curve fit for exponential


persistent options_fit

 if isempty(options_fit)
    %configure the optimset for use with lsqcurvefit
    options_fit = optimset('lsqcurvefit');
   % options_fit.Algorithm = 'levenberg-marquardt';
	options_fit.Algorithm = 'trust-region-reflective';
    options_fit.Display = 'off';
    
    if exist('tol_settings','var')
        if isfield(tol_settings,'MaxFunEvals')
            options_fit.MaxFunEvals = tol_settings.MaxFunEvals;
        end
        if isfield(tol_settings,'MaxIter')
            options_fit.MaxIter = tol_settings.MaxIter;
        end
        if isfield(tol_settings,'TolFun')
            options_fit.TolFun = tol_settings.TolFun;
        end
        if isfield(tol_settings,'TolX')
            options_fit.TolX = tol_settings.TolX;
        end
    end
 end

% do the fitting 
% stored for later: [-inf 0 -inf 0],[inf 7 inf 3]
[param] = lsqcurvefit(@LUXExpFcn_framework,init_params,xdata(cut),ydata(cut),[0 0 -Inf 0],[Inf Inf Inf Inf],options_fit);

%% calculate chi squared / dof

% calculate fitted curve
yfit = LUXExpFcn_framework(param,xdata);
switch 0 % use 0
case 0
	% estimate error for each point on trace
	percent_of_max_for_baseline_error = 0.1;
	
	sd = sqrt(abs(ydata));
	sd((sd == 0) | (abs(ydata) <= percent_of_max_for_baseline_error*max(abs(ydata)))) ...
		 = std(ydata(abs(ydata) <= percent_of_max_for_baseline_error*max(abs(ydata))));
		  
	% calculate chi square and chi square per degree of freedom
	chisq = sum((ydata-yfit).^2./sd.^2);
case 1
	chisq_pieces =(yfit-ydata) + ydata.*log(ydata./yfit);
		chisq_pieces(isnan(chisq_pieces)) = 0;
		chisq_pieces(isinf(chisq_pieces)) = 0;		
 	chisq = 2*sum( chisq_pieces );
 
end

dof = length(xdata)-length(param);

%% plot output

if 0
    figure(783);
    clf; hold on;
    xplotdata = xdata(1):((xdata(2)-xdata(1))/10):xdata(end);
    plot(xplotdata,LUXExpFcn_framework(param,xplotdata),'r','linewidth',2);
    plot(xdata,ydata,'bo','markersize',7);
%    title(['tau fall: ' num2str(param(2)/100e6.*1e9) ' ns --- chi2dof: ' num2str(chisq_per_dof)])
%    xlabel('Time (samples)');
%    ylabel('Amplitude (phe)');
    hold off;
    set(gca,'ysc','log');
%    pause;
 xlabel('samples','fontsize',16);ylabel('phe / sample','fontsize',16);
end



end
