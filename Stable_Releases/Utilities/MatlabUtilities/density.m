function h = density(xdata,ydata,xbine,ybine,varargin)
% h = density(xdata,ydata,xbine,ybine,options) Create a density plot of (x,y) data set
% Uses densitygrid( ) to render output since pcolor( ) and suf( ) get axis
% labels wrong
%
% *bine     A list of BIN EDGES (note lower edge is inclusive), or
%           scalar>0 sets width of bins, based on data min/max, or
%           scalar<0 sets # of bins, based on data min/max
%
% h         Is a handle to the graphics generated
%
% options   'binscale','log10' pair will give you a log10 scale
%           'binscale','percent' pair will give you a percent scale
%           'binscale','fraction' pair will give you a fraction scale
%           These operations can be chained
%
%           'zscale', ## pair will scale data by this factor , before
%           applying binscale operation
%
%
% Useful calls afterwards ....
%
% Use colorbar to get a bar scale
% 
% Use shading('flat') to remove black grid lines
%     shading('faceted') to keep black grid lines
%     shading('interp') to get interpolated colours
%     
% Use caxis([min max]) to select range for min & max colors
%       (Call colorbar again to update scale)
% 
% 030106 rjg First draft v1
% 050401 rjg Added scaling by factor

v = varargin;

bine{1} = xbine; bine{2} = ybine;
xy{1} = xdata(:); xy{2} = ydata(:);

% Bin edges handling for scalar input
for ii=1:2
    if(length(bine{ii}))==1
        minb = min(xy{ii});
        maxb = min(xy{ii});
        if maxb==minb
            maxb = minb+1; % If data is dumb then try to make reasonable choice
        end
        
        if bine{ii}<0 % Use this # of bins, and make sure we include max data
            bine{ii} = linspace(minb,maxb,abs(bine{ii}));
            bine{ii} = [bine{ii} bine{ii}(end)+diff(bine{ii}(1:2))];
        else % Use this as spacing, and make sure we include max data
            bine{ii} = minb:bine{ii}:maxb;
            bine{ii} = [bine{ii} bine{ii}(end)+diff(bine{ii}(1:2))];
        end
    end
end


[xbinc ybinc cts] = Make2DHist(xy{1},xy{2},bine{1},bine{2});

% Check to see if a 'zscale' arugument has been passed
ii=1;
while ii<=length(v)
    if strncmp( v{ii} , 'zscale' , length(v{ii}) )
        cts = v{ii+1} .* cts ;
        v(ii:ii+1)=[];
    else
        ii=ii+1;
    end
end

% THESE COMMAND DONE AFTER COMMANDS ABOVE
% Check to see if a 'binscale' aregument has been passed
ii=1;
while ii<=length(v)
    if strncmp( v{ii} , 'binscale' , length(v{ii}) )
        cts = feval(v{ii+1} , cts );
        v(ii:ii+1)=[];
    else
        ii=ii+1;
    end
end


h2 = densitygrid( bine{1} , bine{2} , cts' ); % requires cts to be transposed
% title(dis('Total counts in grid %7.2g', sum(sum(cts)) ) );

if nargout>0
    h = h2;
end


function y = percent(x)
    y = 100 .* x / sum(sum(x));


function y = fraction(x)
    y = x / sum(sum(x));

