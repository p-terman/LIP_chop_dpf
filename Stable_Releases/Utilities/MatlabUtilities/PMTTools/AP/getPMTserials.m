function serials = getPMTserials
%
% Returns cell array for R8778 PMT serials from LUG
%
% 20100208 CHF - Created

%% Load defaults
xmlsettings = XMLReader('defaultLUGSettings.xml');

%% Construct query
query_string = 'select pmt_serial from lug.PMT where pmt_serial REGEXP ''^B'' order by pmt_serial';
query_string = String_Replacement_Table(query_string);

%% Query LUG

[result_cells query_time_LUG url_alive] = LUGQuery(query_string,xmlsettings);

serials_tmp = result_cells{1};

%Reformat
inds = strfind(serials_tmp,'B');

for ii = 1:numel(inds)
    serials{ii} = serials_tmp(inds(ii):inds(ii)+5);
end