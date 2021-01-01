function h = densitygrid(xe,ye,c)
% h = densitygrid(x,y,c) uses pcolor() to produce density like plot
% The x and y scales have been corrected since pcolor seems to get them
% wrong !
%
% xe, ye     bin edges - pcolor drops the squares at top and right edges when displaying for some
%            reason, so I'm going to add them back
%
% c          matrix of density values, top left is bottom left in density map
%
% Note size(c) = [length(x) length(y)] which means must take transpose from
%                Make2Dhist( ) output
% 030106 rjg First draft

% Add extra edges to c to get over pcolor( ) drop of top & right edges
c = [c;NaN*ones(1,size(c,2))];
c = [c NaN*ones(size(c,1),1)];

h2 = pcolor(xe,ye,c);
if( length(xe)>20 | length(ye)>20 )
    shading('flat');  % Remove the grid if it is going to be  a major distraction
                      % Can get it back with shading('faceted');
end

alpha(0.99); % This seems to fix a bug on MATLAB R13 Mac OS X
            % alpha(1) which is default causes square to be 1 box higher
            % than it should
            
%surf( xe , ye , c );
%view([0 0 1]);

if nargout>0
    h = h2;
end