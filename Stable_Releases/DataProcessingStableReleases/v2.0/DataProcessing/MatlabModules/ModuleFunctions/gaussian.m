function [y, dydX]=gaussian(X,ydata,xdata,ydatasig)
%function [y, dydX]=gaussian(X,xdata)
%function [y]=gaussian(X,xdata)
%function [f]=gaussian(X,ydata,xdata,ydatasig)
%function [f]=gaussian(X,ydata,xdata)
%
%If there are 3 or 4 input parameters, function returns difference in standard 
%deviations between a gaussian function with parameters X (including linear 
%background, possibly) and observed data with values ydata (with optional 
%1-sigma uncertainties ydatasig) at points xdata.  If ydatasig is not input,
%uncertainties are assumed Poisson based on the expected, not actual, number.
%
%If there are two input parameters, function returns expected y-values (and 
%their dependences dy/dX, optionally) on the input parameters of the Gaussian.
%
% X(1) = amplitude
% X(2) = mean
% X(3) = standard deviation
% X(4) = offset of background (optional)
% X(5) = slope of background (optional; if included must include X(4) also)
%
% rws 27 Oct 97

if nargin == 2
	xdata = ydata(:);
elseif nargin >= 3
	xdata = xdata(:);ydata=ydata(:);
elseif nargin < 2
	disp([ 'Not enough arguments in gaussian.m (needs 3, has ' ...
			int2str(nargin) ')' ]);
	pause
end

if length(X)>=3
    if X(3)==0
        disp('WARNING: X(3)=0 in gaussian.m; replacing with eps');
        X(3)=eps;
    end
end	
if length(X) < 3
	disp('Not enough parameters for Gaussian in gaussian.m');
	pause
elseif length(X) == 3
	y=X(1)*exp(-(xdata-X(2)).^2/(2*X(3).^2));
elseif length(X) == 4
	y=X(1)*exp(-(xdata-X(2)).^2/(2*X(3).^2))+X(4);
else
	y=X(1)*exp(-(xdata-X(2)).^2/(2*X(3).^2))+X(4)+X(5)*xdata;
	if length(X) > 5
	disp('WARNING: using only first 5 parameters for Gaussian in gaussian.m');
	end	
end


if nargout == 2					% calculate only if dydX wanted
	dydX(:,1) = exp(-(xdata-X(2)).^2/(2*X(3).^2));
	dydX(:,2) = 2*(xdata-X(2))/(2*X(3).^2)* ...
				X(1).*exp(-(xdata-X(2)).^2/(2*X(3).^2));
	dydX(:,3) = 2*(xdata-X(2)).^2/(2*abs(X(3)).^3)* ...
				X(1).*exp(-(xdata-X(2)).^2/(2*X(3).^2));
	if length(X) > 3
		dydX(:,4) = ones(size(xdata));
	end
	if length(X) > 4	
		dydX(:,5) = xdata;
	end
end

if nargin > 2			% want difference between gaussian and data
	if nargin==3, ydatasig=[];end;
	if isempty(ydatasig)
		ydatasig = sqrt(y);
		ydatasig(y<=0)=eps*ones(sum(y<=0),1);
	end;
	y= (y - ydata)./ydatasig(:) ;	% use (:) to ensure column vector
%	disp(['chisq=' num2str(sum(y.*y)) ' for X=' num2str(X(1)) ' ' ...
%	  num2str(X(2)) ' ' num2str(X(3)) ' ' num2str(X(4)) ]) ;
end
