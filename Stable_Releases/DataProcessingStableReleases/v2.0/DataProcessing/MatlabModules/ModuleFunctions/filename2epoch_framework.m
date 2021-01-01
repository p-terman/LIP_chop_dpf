function epoch_time = filename2epoch(filename_prefix)
%
% epoch_time = filename2epoch(filename_prefix)
%
% This function will convert a filename_prefix into epoch (unix) time
%
% INPUTS:
%   filename_prefix
%
% OUTPUTS:
%   epoch_time      - epoch (unix) time corresponding to filename_prefix date
%                     this has to be used for Slow Control / LUG queries
%
% 20111202 CHF - Created based on Adam Bradley's GetDetectorState.m code
% 20111206 CHF - Minor fix
% 20120313 CHF - Changed hardcoded timezone for -6... kludge
% 20120327 JRV - Now uses LUXGetTimezoneOffset

dt = zeros(1,6);
dt(1) = str2num(filename_prefix(7:10));
dt(2) = str2num(filename_prefix(11:12));
dt(3) = str2num(filename_prefix(13:14));
dt(4) = str2num(filename_prefix(16:17));
dt(5) = str2num(filename_prefix(18:19));
dn = datenum(dt(1),dt(2),dt(3),dt(4),dt(5),dt(6));

secperday = 86400;      % Constant used in date math below
% zerotime is tricky since matlab is using local timezone, not GMT
% This is hard coded for a matlab session using Mountain Daylight Time
% 6 hours behind 

% switch dt(1)
%     case 2010
%         dst = {'Mar 14, 2010 02:00' 'Nov 07, 2010 02:00'};
%     case 2011
%         dst = {'Mar 13, 2011 02:00' 'Nov 06, 2011 02:00'};
%     case 2012
%         dst = {'Mar 11, 2012 02:00' 'Nov 04, 2012 02:00'};
% end;

% this part seems to break due to 'matlab -r' + datenum
%if (dn < datenum(dst{1}))
%    tzhouroffset = -7; % Must be <0, or code needs modifying
%elseif (dn > datenum(dst{1})) && (dn < datenum(dst{2}))
%    tzhouroffset = -6; % Must be <0, or code needs modifying
%elseif (dn > datenum(dst{2}))
%    tzhouroffset = -7; % Must be <0, or code needs modifying
%end

% get the timezone offset
tzhouroffset = LUXGetTimezoneOffset_framework(filename_prefix,'filename');

%zerotime = 719529 + tzhouroffset/24; % MUST hardcode this because the datenum command below crashes in certain circumstances
zerotime = datenum([1969 12 31 24+tzhouroffset 0 0]); 

epoch_time = ceil((dn - zerotime) * secperday); % to be used in SlowControlQuery
