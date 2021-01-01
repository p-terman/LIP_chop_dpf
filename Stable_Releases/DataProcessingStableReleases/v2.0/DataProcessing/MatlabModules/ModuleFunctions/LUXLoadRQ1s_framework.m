function [dp_struct livetime] = LUXLoadRQ1s_framework(file_name, data_path, evt_list, options)

% function [dp_struct livetime] = LUXLoadRQ1s_framework(file_name, data_path, evt_list, options)
%
% Reads the contents of RQ1 files into structure dp.
%
% Takes input filename. filename is either an RQ1 file, a merged RQ1 file,
% or a directory containing either RQ1 or merged RQ1 files
%
% Inputs:
%       file_name  - RQ1 file (.rq1 or .rq1.mrg) to load.
%       data_path - path to filename
%       evt_list  - event list to load. Can be a number, or an array to
%                   load (default: all)
%       options   - structure containing options
%         .RQ1s_list         - cell array with list of RQs to load (default: load all)
%         .no_structure_flag - if set to 1, this function outputs dp_struct
%                              with full RQ1 tree-level names
%                              (i.e. admin.filename_prefix is
%                              admin__filename_prefix )
%                              This option is NEEDED for LUXLoad........
%                              (default: 0)
%        .load_mat           - check for .rq.mat file in ./matfiles/* and if it exists,
%                              load that instead of the .rq file
%
%
% Outputs:
%       dp_struct - 1st pass RQs, structure (i.e. dp.fieldname)
%
%   2008-05-12 CED
%   2009-04-23 LdV Renamed it from ReadNtupleFile to LUX01LoadRQ1s; minor
%                  modifications, limits the number of events read
%   2009-04-24 ChF Using filename_prefix and RQ1s_path separately. Check if
%                  system is gsk-10 to use default path. Outputs dp with
%                  correct dimensions. Looks for .mrg files in folder.
%   2009-04-27 JJC Changed nb_evts to evt_list so that it can be an array
%                  not starting at 1.
%   2009-04-27 JJC Implemented RQ1s_list input. Also with evt_list starting
%                  above event 1, remove empty events.
%   2010-02-04 JJC Fixed endianness so it always reads accurately
%                  regardless of what machine you are on, or the file was
%                  written on.
%   2010-05-27 JJC v5, fopen_endian_safe now reads the settings as well and
%                  so you start at the beginning of the block header string
%                  size, as normal.
%   2011-04-02 DCM Revamped / simplified user interface
%   2011-04-11 JJC Devamped / unsimplified user interface
%   2011-04-11 JJC outputs livetime (block 3)
%   2011-05-17 CHF The function now reshapes the rqs into structure
%                  (takes long variable name with '__' tree level
%                  separators and outputs corresponding fields.
%                  Function also adds .admin field, with settings
%   2011-05-20 CHF Removed RQ1s_list from input, now is a field in
%                  'options', a new input.
%                  This function now has a mode for loading RQ1s with the
%                  full tree-level name ('__' instead of '.') so that
%                  variables can be concated easily.
%   2011-05-23 CHF Added livetime structure in admin.
%   2011-08-15 CHF Since .rq1 files now have the raw RQ1 name, this loader
%                  is now using LUXGenerateRQ1Structure to convert to
%                  tree-level structure
%   2011-08-17 JJC No longer MEXed, just vectorized everything without
%                  looping over events. Old mexed version is now LUXLoadRQ1s_mex.m
%   2011-09-08 JJC no longer casts int8 when freading, now is *int8 -
%                  faster and better with memory management.
%   2011-10-31 LdV Loads RQ1 files in matlab format (.rq1.mat)
%   2011-11-14 JRV dp_struct is always assigned from dp_tmp
%   2011-11-17 PS  rather than crash - automatically unzip files if they are 
%                  in .gzip format, then recompress them after loading
%   2011-11-21 CHF 'flag' is a reserved Matlab function name. Changed to
%                  'myflag', otherwise code in gsk-30 crashes.
%                  Added myflag.gzip=0 clause - otherwise
%                  variable is undefined if statement isn't true
%   2012-06-20 JJC logical data type (1 byte) now cast correctly (needed
%                  for ROOT analysis files).
%   2012-12-11 CHF Livetime is now being read from BLOCK 3, not XML header
%   2013-06-19 pfs ensure transparency when loading .mat files instead of .rq
%                  added new options.load_mat (default is to load rq)
%   2014-05-22 AL  Changed type 'float' in the variable_type_strings_mfriendly
%                  array to 'single', as float is a function and not a type in
%                  Matlab and so would cause crashes when trying to cast
%   2014-06-26 AC  Block 3 (livetime) reading to handle multiple  data lines
%                  (since RQ files can contain multiple DAQ sequences)
%   20190226 PAT Fixed once and for all the issue with one event rqs being
%               loaded in transpose. 
%
dp = [];
livetime = [];
data = cell(1,1);
names = {{}};
ioerror = 0;
filename_prefix = file_name(1:19);

if exist('options','var')
    if ~isempty(options)
        field_names = fieldnames(options);
        for ii=1:length(field_names)
            eval([field_names{ii} ' = options.' field_names{ii} ';']);
        end
    end
end

if ~exist('file_name','var')
    fprintf('LUXLoadRQ1s: ERROR: requires a file name to load\n');
    return
end

myflag.gzip=0;

file_name_full = [data_path filesep file_name];

if ~exist(file_name_full,'file') % then check if it is gz
	fprintf('The requested .rq1 file could not be found on the specified path.');
	fprintf('checking for .gz ...');
	if exist([file_name_full '.gz'],'file')==2
		fprintf('\ngzip -d (de-compressing .rq1 file)...\n');		
		[s,w]=unix(['gzip -d ' file_name_full '.gz']);
		myflag.gzip=1;
    end
end

if ~exist('evt_list','var')
    evt_list = [];
end

if ~exist('RQ1s_list','var')
    RQ1s_list = [];
end
if ~exist('no_structure_flag','var')
    no_structure_flag = 0;
end
if ~exist('load_mat','var')
    load_mat = 0;
end


variable_type_strings={ ...
    'char', ...
    'int8', ...
    'int16', ...
    'int32', ...
    'int64', ...
    'uint8', ...
    'uint16', ...
    'uint32', ...
    'uint64', ...
    'single', ...
    'double', ...
    'float', ...
    'float128', ...
    'logical', ... %
    'unsigned int', ...
    'unsigned short', ...
    'unsigned long long', ...
    'bool', ...
    'int', ...
    'short', ...
    'long long' ...
    };

variable_type_bytes=[1 1 2 4 8 1 2 4 8 4 8 4 16 1 4 2 8 1 4 2 8];

variable_type_strings_mfriendly={ ...
    'char', ...
    'int8', ...
    'int16', ...
    'int32', ...
    'int64', ...
    'uint8', ...
    'uint16', ...
    'uint32', ...
    'uint64', ...
    'single', ...
    'double', ...
    'single' ...
    'float128', ...
    'logical', ... %
    'uint32', ...
    'uint16', ...
    'uint64', ...
    'logical', ...
    'int32', ...
    'int16', ...
    'int64' ...
    };

fprintf('Loading %s\n', file_name_full);


%if strcmp(file_name_full(end-3:end),'.mat') || exist([file_name_full '.mat'],'file')
if exist([data_path 'matfiles/' file_name '.mat'],'file') & load_mat==1
    fprintf('Loading the matlab version\n');
    %if exist([file_name_full '.mat'],'file')
    %    file_name_full = [file_name_full '.mat'];
    %end

    file_name_full = [data_path 'matfiles/' file_name '.mat'];
	dp_struct = load(file_name_full); % dump it into a structure
    
    %if isfield(dp,'livetime_end_samples')
    %    livetime = dp.livetime_end_samples - dp.livetime_latch_samples;
    %elseif isfield(dp.admin,'livetime')
    %    livetime = dp.admin.livetime.end_samples - dp.admin.livetime.latch_samples;
    %else
    %    livetime = 0;
    %end
    %dp_struct = dp;
    
else
    [headerinfo, fid] = LUXSuperLoader_framework(file_name, data_path);
    filename = [data_path filesep file_name];
    
    % livetime should not be in the XML header! CHF
    if isfield(headerinfo,'livetime')
        headerinfo = rmfield(headerinfo,'livetime');
    end
    
    %fid = fopen(filename,'rb');
    %fread(fid,1,'uint32')
    %xmllength = fread(fid,1,'uint32')
    %xmlsettings = char(fread(fid,xmllength,'char')')
    header1length = fread(fid,1,'uint16');
    header1 = char(fread(fid,header1length,'char')');
    nblines1 = fread(fid,1,'int32');
    delims=[0 regexp(header1,'[\.:;]')];
    linewords = (length(delims)-1)/3;
    names{1}=cell(1,linewords);
    vartype=cell(1,linewords);
    varsize=cell(1,linewords);
    varbytes=zeros(1,linewords);
    
    
    totalbytes=0;
    for nw=1:linewords
        names{1}{nw}=header1((delims(3*nw-2)+1):(delims(3*nw-1)-1));
        vartype{nw}=header1((delims(3*nw-1)+1):(delims(3*nw)-1));
        varsize{nw}=sscanf(header1((delims(3*nw)+1):(delims(3*nw+1)-1)),'%d,')';
        found_vartype = find(strcmp(variable_type_strings,vartype{nw}),1,'first');
        varbytes(nw) = variable_type_bytes(found_vartype);
        vartype{nw} = variable_type_strings_mfriendly{found_vartype};
        totalbytes = totalbytes + varbytes(nw)*prod(varsize{nw});
        data{1}.(names{1}{nw}) = cast(zeros([varsize{nw} nblines1]), vartype{nw});
    end
    
    block1data = fread(fid,totalbytes,'*int8');
    
    start=1;
    for nw=1:linewords
        if strcmp(vartype{nw},'logical') || strcmp(vartype{nw},'bool') || strcmp(vartype{nw},'boolean')
            data{1}.(names{1}{nw}) = logical(block1data(start:(start-1+prod(varsize{nw})*varbytes(nw)))');
        elseif strcmp(vartype{nw},'char')
            data{1}.(names{1}{nw}) = cast(block1data(start:(start-1+prod(varsize{nw})*varbytes(nw)))', vartype{nw});
        else
            data{1}.(names{1}{nw}) = typecast(block1data(start:(start-1+prod(varsize{nw})*varbytes(nw))), vartype{nw});
        end
        start = start+prod(varsize{nw})*varbytes(nw);
    end
    
    
    header2length = fread(fid,1,'uint16');
    header2 = char(fread(fid,header2length,'char')');
    nblines2 = fread(fid,1,'int32');
    delims=[0 regexp(header2,'[\.:;]')];
    linewords = (length(delims)-1)/3;
    names{2}=cell(1,linewords);
    vartype=cell(1,linewords);
    varsize=cell(1,linewords);
    varbytes=zeros(1,linewords);
    
    
    totalbytes=0;
    for nw=1:linewords
        names{2}{nw}=header2((delims(3*nw-2)+1):(delims(3*nw-1)-1));
        vartype{nw}=header2((delims(3*nw-1)+1):(delims(3*nw)-1));
        varsize{nw}=sscanf(header2((delims(3*nw)+1):(delims(3*nw+1)-1)),'%d,')';
        found_vartype = find(strcmp(variable_type_strings,vartype{nw}),1,'first');
        varbytes(nw) = variable_type_bytes(found_vartype);
        vartype{nw} = variable_type_strings_mfriendly{found_vartype};
        totalbytes = totalbytes + varbytes(nw)*prod(varsize{nw});
        data{2}.(names{2}{nw}) = cast(zeros([varsize{nw} nblines2]), vartype{nw});
    end
    
    block2data = fread(fid,totalbytes*nblines2,'*int8');
    block2reshaped = reshape(block2data,totalbytes,nblines2);
    
    start=1;
    for nw=1:linewords
        if strcmp(vartype{nw},'char') || strcmp(vartype{nw},'logical')
            data{2}.(names{2}{nw}) = cast(block2reshaped(start:(start-1+prod(varsize{nw})*varbytes(nw)),:)', vartype{nw});
        else
            data{2}.(names{2}{nw}) = reshape(typecast(reshape(block2reshaped(start:(start-1+prod(varsize{nw})*varbytes(nw)),:),1,prod(varsize{nw})*varbytes(nw)*nblines2), vartype{nw}),prod(varsize{nw}),nblines2);
        end
        
        if length(varsize{nw})==2
            data{2}.(names{2}{nw}) = reshape(data{2}.(names{2}{nw}),varsize{nw}(1),varsize{nw}(2),nblines2);
        end
        
        start = start+prod(varsize{nw})*varbytes(nw);
    end
    
    
    header3length = fread(fid,1,'uint16');
    header3 = char(fread(fid,header3length,'char')');
    nblines3 = fread(fid,1,'int32');
    delims=[0 regexp(header3,'[\.:;]')];
    linewords = (length(delims)-1)/3;
    names{3}=cell(1,linewords);
    vartype=cell(1,linewords);
    varsize=cell(1,linewords);
    varbytes=zeros(1,linewords);
    
    
    totalbytes=0;
    for nw=1:linewords
        names{3}{nw}=header3((delims(3*nw-2)+1):(delims(3*nw-1)-1));
        vartype{nw}=header3((delims(3*nw-1)+1):(delims(3*nw)-1));
        varsize{nw}=sscanf(header3((delims(3*nw)+1):(delims(3*nw+1)-1)),'%d,')';
        found_vartype = find(strcmp(variable_type_strings,vartype{nw}),1,'first');
        varbytes(nw) = variable_type_bytes(found_vartype);
        vartype{nw} = variable_type_strings_mfriendly{found_vartype};
        totalbytes = totalbytes + varbytes(nw)*prod(varsize{nw});
        data{3}.(names{3}{nw}) = cast(zeros([varsize{nw} nblines3]), vartype{nw});
    end
    
    block3data = fread(fid,totalbytes*nblines3,'*int8');
    block3reshaped = reshape(block3data,totalbytes,nblines3);
    
    start=1;
    for nw=1:linewords
        if strcmp(vartype{nw},'char')
            data{3}.(names{3}{nw}) = cast(block3reshaped(start:(start-1+prod(varsize{nw})*varbytes(nw)),:)', vartype{nw});
        else
            data{3}.(names{3}{nw}) = reshape(typecast(reshape(block3reshaped(start:(start-1+prod(varsize{nw})*varbytes(nw)),:),1,prod(varsize{nw})*varbytes(nw)*nblines3), vartype{nw}),prod(varsize{nw}),nblines3);
        end
        start = start+prod(varsize{nw})*varbytes(nw);
    end
    
    % Output dp with the correct dimensionality
    fn = fieldnames(data{2}); %Get field names from data{2}
    for ii = 1:numel(fn)
        data{2}.(char(fn{ii})) = squeeze(data{2}.(char(fn{ii})));
        vals = data{2}.(char(fn{ii}));
        sizes = size(vals);
        if length(size(data{2}.(char(fn{ii})))) == 2 %transpose if 2-dimensional
            dp_tmp.(char(fn{ii})) = vals;
        elseif length(size(data{2}.(char(fn{ii})))) == 3 %permute dimensions [1 2 3]->[2 3 1] if 3-dimensional
            dp_tmp.(char(fn{ii})) = squeeze(permute(vals,[1 2 3]));
        end
        if sizes(2) == 1
            dp_tmp.(char(fn{ii})) = vals';
        end
        %dp.(char(fn{ii})) = dp.(char(fn{ii}))(:,evt_list);
    end
    
    
    fclose(fid);
        
    dp_struct = dp_tmp; % just output the raw names
    dp_struct.admin = headerinfo;
    
    % Save livetime FROM BLOCK 3 (NOT HEADERINFO) to admin field
    for ii = 1:linewords
        dp_struct.admin.livetime.(names{3}{ii}) = data{3}.(names{3}{ii});
    end
    
    % this is a *temporary* kludge to fix RQ files written by the root RQ loader
    if ~isfield(dp_struct.admin,'filename_prefix')
        dp_struct.admin.filename_prefix = filename_prefix;
    end
    
    if isfield(dp_struct, 'event_number')
        nb_evts_in_file = length(dp_struct.event_number);
        fields = fieldnames(dp_struct);
        for qq = 1: length(fields)
            current_field = fields{qq};
            if (nb_evts_in_file == 1) && (size(dp_struct.(current_field),1) == 1)...
             && (size(dp_struct.(current_field), 2) ~= 1) && (ndims(dp_struct.(current_field))==2)  
                dp_struct.(current_field) = dp_struct.(current_field)';
                fprintf('\n~~~~~~ transposing dimensionality of rq: %s\n',current_field);
            end
        end
    end

end

if myflag.gzip
	fprintf('gzip (re-compressing .rq1 file) ...\n');		
	[s,w]=unix(['gzip ' data_path '/' file_name]);
end

