function PMTList = getPMTList
% PMTList = getPMTList
%
% For LUX0.1 use. It outputs a cell array with the PMT serial number. The
% index corresponds to the Ch. assignment.
%
% Example: PMTList = getPMTList;
%          PMTList{1} is Ch.1, or BA0208
%
% 20090529 CHF 

PMTList = {'BA0208';'BA0217';'BA0213';'BA0215'};