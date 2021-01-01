function return_status = LUXBinaryRQWriter_framework(settings, dp, filename, output_path, event_list, livetime)
% return_status = LUXBinaryRQWriter_framework(settings, dp, filename, output_path, event_list, livetime)
%
% settings is the structure of settings to write to header
% livetime is the return from LUXGetLivetime for that file
% dp is structure of rq1s to be written to file.
% filename is the name of the file to write.
% output_path is where to write the file.
% To be called from LUXFirstPass
% Write binary RQ file to disk
%
% 2009-04-21 JJC
% 2010-02-01 JJC, DCM
% 2011-05-02 CHF Made event_list an input
%                (previously dp.evt_list)
% 2011-05-17 CHF Now turns dp.() structure into tree-structured variable
%                names, where the tree-level separator is '__'
%                i.e.      pulse__timing__t_10l_samples for
%                          pulse.timing.t_10l_samples
%                LUXLoadRQ1s loads the names as such and turns it back into
%                a structure.
%                Also added daq, event builder and first pass settings into
%                xml header section.
% 2011-08-15 CHF Now writes raw variable names (without the '__' tree-level
%                structure (e.g. baseline_mV instead of
%                pulse__stats__baseline_mV). LUXGenerateRQ1Structure.m now
%                takes care of putting it back into the tree-level format
%
% 2011-09-08 JJC vectorized writing to make is super fast - writes one block
%                cast as bytes ~x300 faster
% 2011-11-15 PS  added error trapping for non-existent output_path
% 2011-11-16 JJC, JRV
%                Fixed typecast bug for char rqs
% 2011-11-17 CHF Minor fix for livetime bugs
% 2012-10-19 CHF Fixed livetime names, changed to livetime_latch_samples, livetime_end_samples
% 2012-11-02 CHF Added livetime as an input. Code checks if it is input, if
%                not, look for it inside dp.admin, for smooth transition.
%                We have to stop storing livetime inside dp.admin soon.
% 2012-11-16 CHF NEW COPY, renamed LUXBinaryRQWriter_framework
%                Made adaptations for new framework
% 2012-12-20 CHF Fixed issue that skipped 1 x nb_evts variables
% 2013-02-16 CHF Converted subfunctions to _framework variety
% 2013-04-13 pfs modified logic in preparing to write, so that we can handle files with only 
%                n=1 events, and files with n=5 events (in addition to of course all other 
%                n>0 cases. have not yet tested how this works on the n=0 case.)
% 2013-04-13 pfs patched so files with zero events no longer crash the writer. 
%                return_status 0 in this case, so we know not to run subsequent modules
% 2013-05-03 JRV Removed keyboard in try/catch line ~318
% 2014-05-22 AL  Changed type 'float' in the variable_type_strings_mfriendly
%                array to 'single', as float is a function and not a type in
%                Matlab and so would cause crashes when trying to cast
% 2014-05-28 AC  Changed livetime block to contain 2 1-dimensional RQs with
%                nEntries=nSequences (as opposed to 2 nSequences-dimensional 
%                RQs with nEntries=nSequences which don't get filled properly -- bug 61)



%%

if ~exist('livetime','var')
    % if livetime is NOT an input, then look inside dp.admin
    % in the future, we need to make all calls to this function have
    % livetime as an input
    
    % livetime struct changed - 2011-11-14 - JRV
    livetime.livetime_latch_samples = dp.admin.livetime.livetime_latch_samples;
    livetime.livetime_end_samples = dp.admin.livetime.livetime_end_samples;
end

%warning off;
%fprintf('Writing Binary RQ File\n');
if isempty(event_list); event_list=0; end; % simple fix for files with zero events, see also EOF

return_status = 0;
first_evt_in_file = event_list(1);
nb_evts_in_file = length(event_list);

if exist(output_path,'dir')~=7 % then make it
    [s,w]=unix(['mkdir -p ' output_path]);
end
RQ1_filename = [output_path, filesep, filename];
bin_fid = fopen(RQ1_filename, 'wb');

if isempty(bin_fid) || bin_fid<0
    fprintf('\nCould not open binary file to write');
    return;
end

dp.admin = settings;

%% Write XML Settings Headers
endianness = hex2dec('01020304');
fwrite(bin_fid, endianness, 'uint32');
xmlstring = MakeXMLString_framework(settings);
fwrite(bin_fid,length(xmlstring),'uint32');
fwrite(bin_fid,xmlstring,'int8');

%%% Write Block 1: File Header %%%

%% Write Block 1 Header %%

b1_header_string = 'dataset_name;char;19;first_evt_in_file;uint32;1;nb_evts_in_file;uint32;1;';
b1_header_string_size = length(b1_header_string);
fwrite(bin_fid, b1_header_string_size,'uint16');
fwrite(bin_fid, b1_header_string, 'char');
b1_nb_data_lines = 1;
fwrite(bin_fid, b1_nb_data_lines, 'int32');

%% Write Block 1 Data Lines %%
fwrite(bin_fid, settings.filename_prefix, 'char');
fwrite(bin_fid, first_evt_in_file, 'uint32');
fwrite(bin_fid, nb_evts_in_file, 'uint32');



%%% Write Block 2: RQs %%%

%% Write Block 2 Header %%

dpn = dp;

var_name = fieldnames(dp);
write_raw_var_name = var_name;

% This is not needed anymore in new framework

% temp_var_names = fieldnames(dp);
%
% for nn = 1:length(temp_var_names)
%     dotinds = strfind(temp_var_names{nn},'.');
%     varname1 = temp_var_names{nn}(dotinds+1:end);
%     % this is the name that will be written into the RQ1 file - no
%     % structure to it.
%     write_raw_var_name{nn} = temp_var_names{nn}(dotinds(end)+1:end);
%     var_name{nn} = strrep(varname1,'.','__');
%     var_name{nn} = strrep(var_name{nn},'()','');
%     eval(sprintf('dpn.%s = [%s];',var_name{nn},temp_var_names{nn}))
% end

nb_RQs = length(var_name);

per_evt = zeros(nb_RQs,1);
var_types = cell(1,nb_RQs);
var_size = cell(1,nb_RQs);
var_size_str = cell(1,nb_RQs);

% this is the old logic that fails if we only have a single event in the file
if 0%for rq=1:nb_RQs
    
    var_types{rq} = class(dpn.(var_name{rq}));
    
    % Added this to make sure RQs with dim 1 x n_evts are written
    sz = size(dpn.(var_name{rq}));
    ndmz = length(sz);
    
    % Fixed an issue here CHF
    if ndmz == 2 && sz(2) == 1;
        dpn.(var_name{rq}) = dpn.(var_name{rq})';
		fprintf('\n****** flipping rq: %s',var_name{rq});
    end
    
    if size(dpn.(var_name{rq}),ndmz) == nb_evts_in_file ...
            && ~strcmp(var_types{rq},'char')
        dimz_end = (ndims(dpn.(var_name{rq}))-1);
        per_evt(rq) = 1;
    else
        dimz_end = ndims(dpn.(var_name{rq}));
        per_evt(rq) = 0;
    end
    
    
    
    for dimz=1:dimz_end
        var_size{rq}(dimz) = size(dpn.(var_name{rq}),dimz);
    end
    
    var_size_str{rq} = '';
    
    for dimz=1:dimz_end
        var_size_str{rq} = strcat(var_size_str{rq}, num2str(var_size{rq}(dimz)),',');
    end
    
    var_size_str{rq} = var_size_str{rq}(1:(end-1));
    
    if per_evt(rq) == 1
        
        if strcmp(var_types(rq),'logical')
            var_types{rq} = 'uint8';
        end
        
        var_header{rq} = sprintf('%s;%s;%s;', write_raw_var_name{rq}, var_types{rq}, var_size_str{rq});
        %         var_header{rq} = sprintf('%s;%s;%s;', var_name{rq}, var_types{rq}, var_size_str{rq});
        
        
    end
    
end

% this is the new logic -pfs
for rq=1:nb_RQs
    
    var_types{rq} = class(dpn.(var_name{rq}));
    
    % transpose 2d array if only 1 event, since in this case craplab takes the unauthorized liberty of flipping it for you (we are just undoing their bs)
    % want event_number 1xn 
    %if (size(dpn.(var_name{rq}),1) == nb_evts_in_file) & (length(size(dpn.(var_name{rq})))<3) 
	% basically, only flip dimensionality of event_number, on the first pass through (!?)
	%if (size(dpn.(var_name{rq}),1) == nb_evts_in_file) & (length(size(dpn.(var_name{rq})))<3)  && (strcmp(var_name{rq}(1:4),'even'))
%{
    PAT removed this problem is fixed at the loader level, writer won't
    actually do anything.
	if (nb_evts_in_file == 1) && (size(dpn.(var_name{rq}),1) == 1)... %PAT modified 190226 so that does what it says
         && (size(dpn.(var_name{rq}), 2) ~= 1) && (ndims(dpn.(var_name{rq}))==2)  % but this really doesn't matter - the problem is in the loader
        dpn.(var_name{rq}) = dpn.(var_name{rq})';
		fprintf('\n~~~~~~ transposing dimensionality of rq: %s\n',var_name{rq});
    end
%}
	if ~isstruct(dpn.(var_name{rq}))
        if nb_evts_in_file > 1
			dimz_end = (ndims(dpn.(var_name{rq}))-1);
		elseif nb_evts_in_file == 1
			dimz_end = (ndims(dpn.(var_name{rq})));		
		else
		end
        per_evt(rq) = 1;
    else
        dimz_end = ndims(dpn.(var_name{rq}));
        per_evt(rq) = 0;
    end
    
    for dimz=1:dimz_end
        var_size{rq}(dimz) = size(dpn.(var_name{rq}),dimz);
    end
    
    var_size_str{rq} = '';
    
    for dimz=1:dimz_end
        var_size_str{rq} = strcat(var_size_str{rq}, num2str(var_size{rq}(dimz)),',');
    end
    
    var_size_str{rq} = var_size_str{rq}(1:(end-1));
    
    if per_evt(rq) == 1
        if strcmp(var_types(rq),'logical')
            var_types{rq} = 'uint8';
        end
        
        var_header{rq} = sprintf('%s;%s;%s;', write_raw_var_name{rq}, var_types{rq}, var_size_str{rq});
    end
    
end
%per_evt
%keyboard

b2_header_string = [var_header{:}];
b2_header_string_size = length(b2_header_string);
fwrite(bin_fid, b2_header_string_size, 'uint16');
fwrite(bin_fid, b2_header_string, 'char');
fwrite(bin_fid, nb_evts_in_file, 'int32');

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
    'int8', ...
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

totalbytes=0;
for nw=1:nb_RQs
    if per_evt(nw)==1
        found_vartype = find(strcmp(variable_type_strings,var_types{nw}),1,'first');
        var_bytes(nw) = variable_type_bytes(found_vartype);
        var_types{nw} = variable_type_strings_mfriendly{found_vartype};
        totalbytes = totalbytes + var_bytes(nw)*prod(var_size{nw});
    end
end
rq_data_block = zeros(totalbytes, nb_evts_in_file, 'int8');


byte_start=1;
byte_end=0;
for rq = 1:nb_RQs
    if per_evt(rq)==1
        byte_end = byte_end + var_bytes(rq)*prod(var_size{rq});
        try
        rq_data_block(byte_start:byte_end,:) = reshape(typecast(reshape(dpn.(var_name{rq}),1,prod(var_size{rq})*nb_evts_in_file),'int8'),(var_bytes(rq)*prod(var_size{rq})),nb_evts_in_file);
        catch ;end % removed keyboard here!
        byte_start = byte_end+1;
        
    end
end
write_me = reshape(rq_data_block,1,numel(rq_data_block));
fwrite(bin_fid,write_me,'int8');






%keyboard
%% Write Block 2 %%
% for evt=1:nb_evts_in_file
%  %   keyboard
%     % Loop over every RQ %
%     for rq=1:nb_RQs
%
%
%
%         RQ_data = dpn.(var_name{rq});
%         if per_evt(rq) == 1
%
%             if length(var_size{rq})==1
%                 for dimz=1:var_size{rq}(1)
%                     fwrite(bin_fid, RQ_data(dimz,evt), var_types{rq});
%                 end
%             else % if multi-dimensional
%                 for dimz1=1:var_size{rq}(end) % last dimension
%
%                     for dimz2=1:var_size{rq}(end-1) % next to last dimension
%
%                         if length(var_size{rq})==2
%                             fwrite(bin_fid, RQ_data(dimz2,dimz1,evt), var_types{rq});
%                         else
%                             for dimz3=1:var_size{rq}(end-2) % next to next to last dimension - only good for three dimensions + events right now
%                                 fwrite(bin_fid, RQ_data(dimz3,dimz2,dimz1,evt), var_types{rq});
%                             end
%                         end
%
%                     end
%
%                 end
%
%             end
%
%         end
%
%
%
%     end
%
% end

%% Write Block 3: Livetime

if ~isempty(livetime)
    %nb_seqs = 1;%settings.settings.daq_settings.sis3301.global.nb_seqs_per_file;
    % it turns out this variable was also called settings... oh well
    nb_seqs = numel(livetime(1).livetime_latch_samples);
    b3_header_string = 'livetime_latch_samples;uint64;1;livetime_end_samples;uint64;1;';
    %b3_header_string = sprintf('livetime_latch_samples;uint64;%d;livetime_end_samples;uint64;%d;', nb_seqs, nb_seqs);      
    b3_header_string_size = length(b3_header_string);
    fwrite(bin_fid, b3_header_string_size, 'uint16');
    fwrite(bin_fid, b3_header_string, 'char');
    fwrite(bin_fid, nb_seqs, 'int32');
    
    for iseq=1:nb_seqs
        fwrite(bin_fid, livetime(1).livetime_latch_samples(iseq),'uint64');
        fwrite(bin_fid, livetime(1).livetime_end_samples(iseq),'uint64');
    end
end

%% Set file permissions to 755 for lux user
[status, result] = unix(['chmod 755 ' RQ1_filename]);
if status ~= 0
    disp(['Warning: unable to set file permissions for ' RQ1_filename]);
    fprintf('Message: %s', result);
end

clear var_size;
clear per_evt;
clear var_header;
clear var_size_str;
clear dimz*;
clear var_name;
clear var_types;
clear RQ_data;
fclose(bin_fid);
if all(event_list==0)
	fprintf('~~~~~~ LUXBinaryRQWriter_framework: exist status 0 !! (no events in file)\n');
	return_status = 0; % fail
else
	return_status = 1;
end

fprintf('\nwrote file:\n %s\n',[output_path filesep filename]); drawnow;


