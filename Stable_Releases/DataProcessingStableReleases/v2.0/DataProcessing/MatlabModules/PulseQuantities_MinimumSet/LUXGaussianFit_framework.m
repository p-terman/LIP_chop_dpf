function [param chisq dof] = LUXGaussianFit_framework( xdata, ydata, tol_settings)
% This function fits a gaussian (S2) to pulses in the first pass.
% [param chisq dof] = LUXGaussianFit_framework(xdata, ydata)
%
% Inputs:
%          xdata - x data (samples)
%          ydata - data trace
%   tol_settings - tolerance settings structure (optional)
%
% Outputs:
%          param - amplitude, mean, sigma
%          chisq - chi squared for fit
%            dof - degrees of freedom for fit
%
%
% 20110315 - JRV - Created
% 20110706 - JRV - Now calls LUXGaussFcn, instead of defining Gaussian
%                  internally
% 20110907 - JRV - Now uses global vars to avoid multiple calls to optimset
% 20111029 - CHF - Output is now separating chisq and dof, to match new RQ1 variables
% 20111118 - JRV - LUXGaussFcn now calculates Jacobian. This is ~30%
%                  faster after taking advantage.
% 20130512 - JRV - Reduced tolerance values to 1e-2
%                  Set maximum iters to 100
% 20130529 - JRV - Now accepts tolerance settings

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
    init_params(2) = xdata(ymax_loc);
else
    init_params(1) = ymin;
    init_params(2) = xdata(ymin_loc);
end

% standard deviation

init_params(3) = 1e1;

%% use matlab non-linear least squares curve fit for gaussian

persistent options_fit_gauss

if isempty(options_fit_gauss)
    %configure the optimset for use with lsqcurvefit
    options_fit_gauss = optimset('lsqcurvefit');
    options_fit_gauss.Algorithm = 'levenberg-marquardt';
    options_fit_gauss.Display = 'off';
    options_fit_gauss.Jacobian = 'off';
    
    if exist('tol_settings','var')
        if isfield(tol_settings,'MaxFunEvals')
            options_fit_gauss.MaxFunEvals = tol_settings.MaxFunEvals;
        end
        if isfield(tol_settings,'MaxIter')
            options_fit_gauss.MaxIter = tol_settings.MaxIter;
        end
        if isfield(tol_settings,'TolFun')
            options_fit_gauss.TolFun = tol_settings.TolFun;
        end
        if isfield(tol_settings,'TolX')
            options_fit_gauss.TolX = tol_settings.TolX;
        end
    end
end


param = lsqcurvefit(@LUXGaussFcn_framework,init_params,xdata,ydata,[],[],options_fit_gauss);
param(3) = abs(param(3));

%% calculate chi squared / dof

% calculate fitted curve
yfit = LUXGaussFcn_framework(param,xdata);

% estimate error for each point on trace
percent_of_max_for_baseline_error = 0.10;

sd = sqrt(abs(yfit));
sd((sd == 0) | (abs(ydata) <= percent_of_max_for_baseline_error*max(abs(ydata)))) ...
     = std(ydata(abs(ydata) <= percent_of_max_for_baseline_error*max(abs(ydata))));
      
% calculate chi square and degree of freedom
chisq = sum((ydata-yfit).^2./sd.^2);
dof = length(xdata)-length(param);

%% plot output for testing

if 0
    figure(3);
    clf; hold on;
    xplotdata = xdata(1):0.01:xdata(end);
    plot(xdata,ydata,'bo','markersize',5);
    plot(xplotdata,LUXGaussFcn_framework(param,xplotdata),'r','linewidth',2);
    title([' chi2: ' num2str(chisq) ' dof: ' num2str(dof)])
    hold off;
end


