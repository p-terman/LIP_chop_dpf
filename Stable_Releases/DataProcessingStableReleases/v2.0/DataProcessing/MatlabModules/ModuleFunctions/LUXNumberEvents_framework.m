function [evt_list] = LUXNumberEvents_framework(filename,data_path)
% [evt_list] = LUXNumberEvents_framework(filename,data_path)     
%
%
% Input: filename - name of .evt file that you want to look at
%        data_path - path to the aforementioned file
% Output: evt_list - Array with event numbers
%
% 2010-05-27 jjc - created
% 2011-01-11 jjc - finally updated help to reflect actual inputs and
% outputs.
% 2011-08-15 jjc - now works with new file format for .evt files
% 2013-02-16 CHF Converted subfunctions to _framework variety
% 2013-04-02 pfs - added back the fclose(fid). this may solve the "too many files open"
%                  error that I have been seeing
%
%%

options.load_settings=0;
gunzip_loc = './';

 temp_loc=0;
    if strcmp(filename((end-1):end),'gz')
        fprintf('gunzipping file %d...\n',str2num(filename(22:30)));
        %[stat mess] = fileattrib([data_path '/' files_to_load(ff).name]);
        %if mess.GroupWrite==1
        %    gunzip([data_path '/' files_to_load(ff).name]);
        %    files_to_load(ff).name = files_to_load(ff).name(1:(end-3));
        %    temp_loc=0;
        %else
        %fprintf('I do not have write privileges for gunzipping in the data_path,\nso I will put the gunzipped file in your current directory.\nI will delete it when I am done.\n');
        if ~exist([gunzip_loc '/' filename],'file')
            copyfile([data_path '/' filename], gunzip_loc);
        end
        temp_loc=1;
        gunzip([gunzip_loc filename]);
        filename = filename(1:(end-3));
        %end
    end
    
    if temp_loc==1
        [settings, fid] = LUXSuperLoader_framework(filename,gunzip_loc,options);
    else
        [settings, fid] = LUXSuperLoader_framework(filename,data_path,options);
    end


%[settings fid] = LUXSuperLoader(filename, data_path, options);

file_header = fread(fid,4,'uint32');
nb_evts_in_file = file_header(4);
   
evt_index_block = fread(fid,2*nb_evts_in_file,'uint32');
events_in_file_numbers = evt_index_block(logical(mod(1:(2*nb_evts_in_file),2)));
%events_in_file_positions = evt_index_block(~logical(mod(1:(2*nb_evts_in_file),2))); 

evt_list = events_in_file_numbers;

    %%  "Counts" the number of events for the a single file
    
    %filefid = OpenBinaryFile([datapath '/' filename]);
  %  options.load_settings=0;
  %  [settings fid] = LUXSuperLoader(filename, data_path, options);
  %  if fid>0
  %      fseek(fid,4,'cof');
  %      header_gid = fread(fid, 4, 'uint32'); % Note that the first ulong (endianess) was already read by OpenBinaryFile
  %      evt_list = header_gid(3) + [1:(header_gid(4))]; % event numbering starts at 0
  %      status = fclose(fid);
   % end
   fclose(fid);
end

