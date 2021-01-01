function [AP dates] = getAPfromPMTdb(PMTserial,limit,bias)
%
% Inputs PMTserial, outputs AP structure with fields:
%   AP.H
%   AP.He
%   AP.NO
%   AP.Ar
%   AP.Xe
%
% 20100208 CHF - Created

if nargin < 3
    biascheck = 0;
elseif nargin < 2
    limit = 1;
    biascheck = 0;
else
    biascheck = 1;
end

%% Load defaults
xmlsettings = XMLReader('defaultLUGSettings.xml');

%% Construct query

if biascheck
    query_string1 = sprintf('select AP_H, AP_He, AP_NO, AP_Ar, AP_Xe from lug.AP where pmt_serial = ''%s'' and pmt_bias_V = %d order by action_date desc limit %d',PMTserial,bias,limit);
    query_string2 = sprintf('select action_date from lug.AP where pmt_serial = ''%s'' and pmt_bias_V = %d order by action_date desc limit %d',PMTserial,bias,limit);
    query_string3 = sprintf('select pmt_bias_V from lug.AP where pmt_serial = ''%s'' and pmt_bias_V = %d order by action_date desc limit %d',PMTserial,bias,limit); %kludge, but works
else
    query_string1 = sprintf('select AP_H, AP_He, AP_NO, AP_Ar, AP_Xe from lug.AP where pmt_serial = ''%s'' order by action_date desc limit %d',PMTserial,limit);
    query_string2 = sprintf('select action_date from lug.AP where pmt_serial = ''%s'' order by action_date desc limit %d',PMTserial,limit);
    query_string3 = sprintf('select pmt_bias_V from lug.AP where pmt_serial = ''%s'' order by action_date desc limit %d',PMTserial,limit);
end

query_string1 = String_Replacement_Table(query_string1);
query_string2 = String_Replacement_Table(query_string2);
query_string3 = String_Replacement_Table(query_string3);

%% Query LUG

[result_cells1 query_time_LUG url_alive] = LUGQuery(query_string1,xmlsettings);
[result_cells2 query_time_LUG url_alive] = LUGQuery(query_string2,xmlsettings);
[result_cells3 query_time_LUG url_alive] = LUGQuery(query_string3,xmlsettings);

%% Reformat

AP_raw = str2num(result_cells1{1});
d = result_cells2{1};

b = str2num(result_cells3{1});

inds = strfind(d,' ');

for ww = 1:length(AP_raw)/5
    AP(ww).H  = AP_raw(1 + (ww-1)*5);
    AP(ww).He = AP_raw(2 + (ww-1)*5);
    AP(ww).NO = AP_raw(3 + (ww-1)*5);
    AP(ww).Ar = AP_raw(4 + (ww-1)*5);
    AP(ww).Xe = AP_raw(5 + (ww-1)*5);
    AP(ww).bias = b(ww);
end

for ii = 1:numel(inds)
    dates{ii}.label = d(inds(ii)-10:inds(ii)+8);
    dates{ii}.datenum = datenum(dates{ii}.label);
end


