% function hrs = filename2hrs(filename)
%
% Transforms dataset number to hours since Jan, 1 2009 00:00
%
% CHF
function hrs = filename2hrs(filename)
%%

ind = strfind(filename,'_');

yyyy = str2double(filename(ind+1:ind+4));
mm = str2double(filename(ind+5:ind+6));
dd = str2double(filename(ind+7:ind+8));
hh = str2double(filename(ind+10:ind+11));
nn = str2double(filename(ind+12:ind+13));

hrs = ((yyyy-2009)*365*24) + (mm*30*24) + (dd*24) + hh + nn/60;