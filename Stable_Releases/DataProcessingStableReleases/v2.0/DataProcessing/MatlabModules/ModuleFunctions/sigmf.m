function h = sigmf(x, p)

% a nice smooth step function, like erf
% except here, the offset p(2) and width p(1)
% are easy inputs.
%     -SAH   april 15, 2013

h = (1 + exp(-p(1)*(x-p(2)))).^(-1);
