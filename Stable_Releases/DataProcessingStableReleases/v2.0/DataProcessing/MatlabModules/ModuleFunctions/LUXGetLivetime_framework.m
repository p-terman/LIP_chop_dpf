function livetime = LUXGetLivetime_framework(filename, data_path)
% [livetime] LUXGetLivetime(filename, data_path);
% returns livetime.latch, .end
%	latch - is the timestamp of when the acquisition began
%	end - is the tiemstamp of the last sample that ended the acquisition (multi-event mode)
% 2009-06-24 JJC created
% 2011-06-27 LdV Creates a blank livetime before it begins to read the livetimes
% 2011-08-15 JJC modified to work with new evt file format (index of events
% in header)
% 2011-09-06 JJC can read .evt or .dat file
% 2012-11-02 CHF Corrected livetime names to include _samples, as agreed.
% 2012-12-18 PfS added some verbosity in case file_header = []
% 2012-12-20 CHF Copied as new file name (_framework). Changed outputs to
%                livetime_latch_samples, livetime_end_samples
% 2013-02-16 CHF Converted subfunctions to _framework variety

livetime = [];
options.load_settings = 0;
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
        [s, fid] = LUXSuperLoader_framework(filename,gunzip_loc,options);
    else
        [s, fid] = LUXSuperLoader_framework(filename,data_path,options);
    end

%[s fid] = LUXSuperLoader(filename,data_path,options);
%keyboard
if strcmp(filename((end-2):(end)),'evt')
file_header = fread(fid,4,'uint32');
if isempty(file_header); disp('\n*** file_header = []. This is going to crash LUXGetLivetime.m. You can probably avert the crash by unzipping the .evt files...'); end
nb_evts_in_file = file_header(4);
   
fseek(fid,nb_evts_in_file*8,'cof'); % skip over index

elseif strcmp(filename((end-2):(end)),'dat')
    file_header = fread(fid,8,'uint32');
elseif strcmp(filename((end-2):(end)),'rq1')
    fprintf('Better off just loading the data file, see LUXLoadRQ1s\n');
    return;
else
   fprintf('I do not recognize the filetype you are trying to open.\nChoose a .dat or .evt file\n');
   return;
end
nb_seqs_in_file = fread(fid,1,'uint16');
 for i_seq=1:nb_seqs_in_file
     livetime.livetime_latch_samples(i_seq) = fread(fid,1,'uint64');
     livetime.livetime_end_samples(i_seq) = fread(fid,1,'uint64');
 end
 fclose(fid);

