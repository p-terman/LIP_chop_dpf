% LUXFirstPass_BinaryRQWriter
% 
% DO NOT CALL THIS MANUALLY, IT DOESN'T MAKE SENSE
% 
% To be called from LUXFirstPass
% Write binary RQ file to disk
% 
% 2010-02-01 JJC, DCM
% 


%
% %%%%% Binary RQ Writer 2009-04-21 JJC %%%%%%
%
%


fprintf('Writing Binary RQ File\n');

if nb_evts_in_file==0
    first_evt_in_file = 0;
else
    first_evt_in_file = evt_list(1);
end

RQ1_filename = sprintf([output_path, filesep, output_filename]);
bin_fid = fopen(RQ1_filename, 'wb');

if isempty(bin_fid) || bin_fid<0
    fprintf('Could not open binary file to write');
    return
end

%% Write XML Settings Headers
endianness = hex2dec('01020304');
fwrite(bin_fid, endianness, 'uint32');
xmlstring = MakeXMLString(s);
fwrite(bin_fid,length(xmlstring),'uint32');
fwrite(bin_fid,xmlstring,'int8');


%%% Write Block 1: File Header %%%

%% Write Block 1 Header %%

%endianness = hex2dec('01020304');
%fwrite(bin_fid, endianness, 'uint32');

b1_header_string = 'dataset_name;char;19;first_evt_in_file;uint32;1;nb_evts_in_file;uint32;1;';
b1_header_string_size = length(b1_header_string);
fwrite(bin_fid, b1_header_string_size,'uint16');
fwrite(bin_fid, b1_header_string, 'char');
b1_nb_data_lines = 1;
fwrite(bin_fid, b1_nb_data_lines, 'int32');

%% Write Block 1 Data Lines %%
fwrite(bin_fid, filename_prefix, 'char');
fwrite(bin_fid, first_evt_in_file, 'uint32');
fwrite(bin_fid, nb_evts_in_file, 'uint32');

%%% Write Block 2: RQs %%%

%% Write Block 2 Header %%

var_name = fieldnames(dp{1}); % get field names
nb_RQs = length(var_name);
per_evt = zeros(nb_RQs,1);
for rq=1:nb_RQs
    var_types{rq} = class(getfield(dp{1},cell2mat(var_name(rq)))) ;
    if size(getfield(dp{1},cell2mat(var_name(rq))),ndims(getfield(dp{1},cell2mat(var_name(rq))))) == nb_evts_in_file
        dimz_end = (ndims(getfield(dp{1},cell2mat(var_name(rq))))-1);
        per_evt(rq) = 1;
        %var_name(rq)
    else
        dimz_end = ndims(getfield(dp{1},cell2mat(var_name(rq))));
        per_evt(rq) = 0;
    end
    for dimz=1:dimz_end
        var_size{rq}(dimz) = size(getfield(dp{1},cell2mat(var_name(rq))),dimz);
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
        var_header{rq} = sprintf('%s;%s;%s;', var_name{rq}, var_types{rq}, var_size_str{rq});
    end
end


b2_header_string = [var_header{:}];
b2_header_string_size = length(b2_header_string);
fwrite(bin_fid, b2_header_string_size, 'uint16');
fwrite(bin_fid, b2_header_string, 'char');
fwrite(bin_fid, nb_evts_in_file, 'int32');

%% Write Block 2 %%
for evt=1:nb_evts_in_file
    % Loop over every RQ %
    for rq=1:nb_RQs
        
        
        
        RQ_data = getfield(dp{1},cell2mat(var_name(rq)));
        if per_evt(rq) == 1
            
            if length(var_size{rq})==1
                for dimz=1:var_size{rq}(1)
                    fwrite(bin_fid, RQ_data(dimz,evt), var_types{rq});
                end
            else % if multi-dimensional
                for dimz1=1:var_size{rq}(end) % last dimension
                    
                    for dimz2=1:var_size{rq}(end-1) % next to last dimension
                        
                        if length(var_size{rq})==2
                            fwrite(bin_fid, RQ_data(dimz2,dimz1,evt), var_types{rq});
                        else
                            for dimz3=1:var_size{rq}(end-2) % next to next to last dimension - only good for three dimensions + events right now
                                fwrite(bin_fid, RQ_data(dimz3,dimz2,dimz1,evt), var_types{rq});
                            end
                        end
                        
                    end
                    
                end
                
            end
        end
        
        
        
    end
    
end

%% Write Block 3: Livetime
if ~isempty(livetime)
    nb_seqs = s.daq_settings.sis3301.global.nb_seqs_per_file;
    b3_header_string = 'timestamp_latch;uint64;1;timestamp_end;uint64;1;';
    b3_header_string_size = length(b3_header_string);
    fwrite(bin_fid, b3_header_string_size, 'uint16');
    fwrite(bin_fid, b3_header_string, 'char');
    fwrite(bin_fid, nb_seqs, 'int32');
    
    for iseq=1:nb_seqs
        fwrite(bin_fid, livetime(1).latch(iseq),'uint64');
        fwrite(bin_fid, livetime(1).end(iseq),'uint64');
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



