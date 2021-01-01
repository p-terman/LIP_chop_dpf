function [evt_list] = LUXNumberEvents(filename,data_path)
% nb_evts = LUXNumberEvents(xml_settings)       
%
%
% Input: structure with the dataset settings, created with LUXLoadSettings
% Output: Number of Events in Dataset
%
% 2010-05-27 jjc - created


    %%  "Counts" the number of events for the a single file
    
    %filefid = OpenBinaryFile([datapath '/' filename]);
    options.load_settings=0;
    [settings fid] = LUXSuperLoader(filename, data_path, options);
    if fid>0
        fseek(fid,4,'cof');
        header_gid = fread(fid, 4, 'uint32'); % Note that the first ulong (endianess) was already read by OpenBinaryFile
        evt_list = header_gid(3) + [1:(header_gid(4))]; % event numbering starts at 0
        status = fclose(fid);
    end
end

