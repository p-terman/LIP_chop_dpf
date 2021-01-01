% [cutmouse cutdata] = mousecut(xdata,ydata,cut,fignum)
%
% This function will allow you to make a logical cut to data based on a polygon drawn by the mouse on the figure window.
%
% INPUTS:
%               xdata, ydata
%    [optional] cut - apply cut to data. This preserves to global cut size, so you can combine this cut with others
%    [optional] fignum - plot the data for mousecutting in figure(fignum)
%
% OUTPUTS:
%               cutmouse - the logical cut inside the drawn polygon (shown in green in plot)
%
% 20090730 - CHF
%            Many thanks to Doug Hull from Mathworks for his insight
%
function [cutmouse cutdata] = mousecut(xdata,ydata,cut,fignum)

warning off
if nargin == 3
    fignum = 5555;
elseif nargin == 2
    fignum = 5555;
    cut = xdata./xdata > -1;
end

figure(fignum); clf
plt(xdata(cut),ydata(cut),'k.')

% answer1 = input('Log-Log? Y = 1: ');
% 
% if answer1 == 1
%     set(gca,'ysc','log')
%     set(gca,'xsc','log')
% end
% 
answer2 = input('Change axis? Y = 1: ');

if answer2 == 1
    newax = input('Enter new axis in [x1 x2 y1 y2] format: ');
    axis(newax);
end

fprintf('Click with your mouse on the vertices that bound the cut\n')
[x,y] = ginput;
k=convhull(x,y);
temp_cut=inpolygon(xdata,ydata,x(k),y(k));
cutmouse = cut & temp_cut;
hold on;
plot(x(k),y(k),'b-.');
plt(xdata(cutmouse),ydata(cutmouse),'g.');

if nargout == 2
    cutdata.x = x;
    cutdata.y = y;
    cutdata.k = k;
end

warning on