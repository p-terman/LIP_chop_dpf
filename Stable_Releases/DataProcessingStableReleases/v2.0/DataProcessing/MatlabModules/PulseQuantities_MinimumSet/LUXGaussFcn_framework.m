function y = LUXGaussFcn_framework(param,x)
%
%  Inputs:
%     param - [A, mean, sigma]
%         x - x data
%
% Outputs:
%         y - y data
%         J - Jacobian
%
% For LUXGaussianFit_framework
%
% 20110531 - JRV - Created
% 20111118 - JRV - Added Jacobian output

% ensure proper type

param = double(param);
x = double(x);

% calculate
A = param(1);
mean = param(2);
sigma = param(3);

y = A.*exp(-(x-mean).^2./(2.*sigma.^2));
        
if 0
    % calculate Jacobian 
    J = zeros(length(x),length(param));
    J(:,1) = exp(-(x-mean).^2./(2.*sigma.^2));
    J(:,2) = A.*(x-mean)./sigma.^2 .* J(:,1)';
    J(:,3) = J(:,2)'.*(x-mean)./sigma;
end

end