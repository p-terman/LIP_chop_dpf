function [boolout] = inrange( in , varargin )
%function [boolout] = inrange( in , varargin )
% rjg
% rjg 030221 allow more complex argument, ie if arg element is 2 value vector then
% this is a range in itself. Multiple arguments can be used for multiple OR'd ranges
% Still supports older format: First scalar arg is assumed a min, 2nd scalar is assumed max
%

flagmin = 1;
v = varargin;
boolout = false(size(in));

while length(v)>0
    switch length(v{1})
        case 1
            switch flagmin
                case 1
                    boolout = in>=v{1};
                    flagmin = 2;
                case 2
                    boolout = boolout & in<v{1};
                    flagmin = 3;
                otherwise
                    boolout = []; dis('inrange You can only have two scalar arguments');
            end
        case 2
            boolout = boolout | (in>=v{1}(1) & in<v{1}(2));
        otherwise
            boolout = []; dis('inrange format error');
    end
    v(1) = []; % Clear last argument
end

