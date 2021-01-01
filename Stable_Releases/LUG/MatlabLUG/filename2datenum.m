% function out = filename2datenum(filename)
%
% Transforms dataset number to hours since Jan, 1 2009 00:00
%
% CHF
function out = filename2datenum(filename)
%%

ind = strfind(filename,'_');

yyyy = str2double(filename(ind+1:ind+4));
mm = str2double(filename(ind+5:ind+6));
dd = str2double(filename(ind+7:ind+8));
hh = str2double(filename(ind+10:ind+11));
nn = str2double(filename(ind+12:ind+13));
ss = 0;

out = datenum(yyyy,mm,dd,hh,nn,ss);