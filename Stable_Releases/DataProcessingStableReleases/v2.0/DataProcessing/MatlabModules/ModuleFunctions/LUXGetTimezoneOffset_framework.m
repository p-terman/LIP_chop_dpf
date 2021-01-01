function [ timezone_offset_hours ] = LUXGetTimezoneOffset( input , format)
%function [ timezone_offset_hours ] = LUXGetTimezoneOffset( input )
%
% Inputs:
%
%             input (optional) - Either filename prefix, unix, or datenum
%                                If no input then use current time
%   format (required if input) - 'filename', 'unix', or 'datenum'
%
% Outputs:
%
%   timezone_offset_hours - Lead, SD timezone offset in hours (-6 or -7)
%
%
% This function returns the timezone offset for Lead, SD given a date in
% unix time or a filename prefix. This function makes limited use of
% datenum to avoid CCV complications.
%
% ATTENTION: When using the filename prefix input there is an ambiguity
% between 2-3am on the first Sunday in Nov. This function gives all dates
% within this range a default offset of -7 hours. This could be avoided by
% querying the daq gui database for the acquistion start time in unix
% format; however, not all acquisitions are guaranteed to be in that
% database.
%
% 2012-03-26 - JRV - Created
% 2012-04-05 - JRV - Datenum day-of-week doesn't work on CCV - fixed.
% 2012-04-16 - JRV - Added support for datenum - format string required

%% Main function

% if no input then use current time
if nargin < 1
    [~,utstr] = unix('date +%s');
    input = str2num(utstr);
    format = 'unix';
end

if nargin == 1
    fprintf('ERROR: If an input exists you must include the format\n');
    return;
end

% check if we're currently in dst
try
    if strcmp(format,'filename')
        is_curr_dst = isdst_fnp(input);
    elseif strcmp(format,'unix')
        % use -7 for default non_dst_offset (MT)
        non_dst_offset = -7;
        is_curr_dst = isdst_unix(input, non_dst_offset);
    elseif strcmp(format,'datenum')
        is_curr_dst = isdst_datenum(input);
    end
catch
    fprintf('ERROR: Input incorrect.\n\n');
    timezone_offset_hours = [];
    return
end

% return timezone offset for Lead, SD
if is_curr_dst
    timezone_offset_hours = -6;
else
    timezone_offset_hours = -7;
end



%% Internal functions

% TODO: A lot of this code is reused; should be one function...

    function [ out ] = isdst_unix( unix_time, non_dst_offset )

        % unix time format
        year = str2num(datestr(datenum([1970 1 1 non_dst_offset 0 unix_time]),'yyyy'));
        month = str2num(datestr(datenum([1970 1 1 non_dst_offset 0 unix_time]),'mm'));
        day = str2num(datestr(datenum([1970 1 1 non_dst_offset 0 unix_time]),'dd'));
        hour = str2num(datestr(datenum([1970 1 1 non_dst_offset 0 unix_time]),'HH'));
        %dow = datestr(datenum([1970 1 1 non_dst_offset 0 unix_time]),'ddd')
        dow_num = get_dow(year,month,day);
        
        % get the number of the day of the week
        %dow_table = {'Sun','Mon', 'Tue', 'Wed','Thu','Fri','Sat'};
        %dow_num = find(strcmp(dow_table,dow))-1;

        % definitely not in dst
        if ( month < 3 || month > 11 )
            out = false;
        end

        % definitely in dst
        if (month > 3 && month < 11)
            out = true;
        end

        % case for March
        previous_sunday = day - dow_num;
        if (month == 3 && previous_sunday >=8)
            if ( dow_num == 0 && previous_sunday <= 14 && hour < 2)
                out = false;
            else
                out = true;
            end
        elseif (month == 3)
            out = false;
        end

        % case for November
        if (month == 11)
            if (previous_sunday < 0)
                out = true;
            elseif ( dow_num == 0 && day <= 7 && hour < (2-1) ) % (2-1) because of init. dst assumption 
                out = true;
            else
                out = false;
            end
        end
    end

    function [ out ] = isdst_fnp( filename_prefix )

        % filename prefix format
        year = str2num(filename_prefix(7:10));
        month = str2num(filename_prefix(11:12));
        day = str2num(filename_prefix(13:14));
        hour = str2num(filename_prefix(16:17));
        min = str2num(filename_prefix(18:19));
        %dow = datestr(datenum([year month day hour min 0]),'ddd')
        dow_num = get_dow(year,month,day);
        
        % get the number of the day of the week
        %dow_table = {'Sun','Mon', 'Tue', 'Wed','Thu','Fri','Sat'};
        %dow_num = find(strcmp(dow_table,dow))-1;

        % definitely not in dst
        if ( month < 3 || month > 11 )
            out = false;
        end

        % definitely in dst
        if (month > 3 && month < 11)
            out = true;
        end

        % case for March
        previous_sunday = day - dow_num;
        if (month == 3 && previous_sunday >=8)
            if ( dow_num == 0 && previous_sunday <= 14 && hour < 2)
                out = false;
            else
                out = true;
            end
        elseif (month == 3)
            out = false;
        end

        % case for November
        if (month == 11)
            if (previous_sunday < 0)
                out = true;
            elseif ( dow_num == 0 && day <= 7 && hour < 2 ) 
                out = true;
            else
                out = false;
            end
        end
    end

    function [ out ] = isdst_datenum( matlab_datenum )

        % matlab datenum time format
        year = str2num(datestr(matlab_datenum,'yyyy'));
        month = str2num(datestr(matlab_datenum,'mm'));
        day = str2num(datestr(matlab_datenum,'dd'));
        hour = str2num(datestr(matlab_datenum,'HH'));
        dow_num = get_dow(year,month,day);
        
        % get the number of the day of the week
        %dow_table = {'Sun','Mon', 'Tue', 'Wed','Thu','Fri','Sat'};
        %dow_num = find(strcmp(dow_table,dow))-1;

        % definitely not in dst
        if ( month < 3 || month > 11 )
            out = false;
        end

        % definitely in dst
        if (month > 3 && month < 11)
            out = true;
        end

        % case for March
        previous_sunday = day - dow_num;
        if (month == 3 && previous_sunday >=8)
            if ( dow_num == 0 && previous_sunday <= 14 && hour < 2)
                out = false;
            else
                out = true;
            end
        elseif (month == 3)
            out = false;
        end

        % case for November
        if (month == 11)
            if (previous_sunday < 0)
                out = true;
            elseif ( dow_num == 0 && day <= 7 && hour < 2 ) 
                out = true;
            else
                out = false;
            end
        end
    end

    % datestr dow doesn't work on CCV...
    % https://en.wikipedia.org/wiki/Determination_of_the_day_of_the_week
    function dow = get_dow(y,m,d)
        t = [0, 3, 2, 5, 0, 3, 5, 1, 4, 6, 2, 4];
        y = y - (m < 3);
        dow = mod((y + floor(y/4) - floor(y/100) + floor(y/400) + t(m) + d),7);
    end
                

end

