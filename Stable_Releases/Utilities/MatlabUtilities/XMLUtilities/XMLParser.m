% function xml = XMLParser(filename)
%
% Reads the xml file and outputs it as a structure
%
% Input:  file   -   dataset structure with name and file fields filled
% Output: xml   -   return dataset structure with xml fields
%
% Eg -- <fieldname> <dataname1> datastring1 </dataname1> <dataname2> datastring2 </dataname2> </fieldname>
% creates fields dataset.xml.fieldname.dataname1 = datastring1
%                dataset.xml.fieldname.dataname2 = datastring2
%
% NOTE:  content of all xml subfields are converted to numerics is possible
%        THIS MEANS:  12/13/06 12:02:31 is interpreted as the array
%        [.154 12 14 16 18 20 22 24 26 28 30]
% NOTE2:  if a fieldname appears multiple times at the same level, an 
%         array of structs, or cell array of strings is formed.
% 10/03/2005, ED
% Release Version 1.0

function [xml restofstring] = XMLParser(xmlstring, start_tag)

if nargin < 2
    start_tag = '';
end;

% initialize variables
currentfieldname = '';
current_content = xmlstring;
dataholder = true;

%% recursive loop
while true
    
    xml_bits = regexp([' ' current_content],'^((.(?!</?\w+>))*.?)<(/?\w+)>(.*)$','tokens');
    
    if isempty(xml_bits)
        if ~exist('xml','var')
            xml = struct();
        end;
        if ~exist('restofstring','var')
            restofstring = '';
        end;
        return;
    end;
    
    if strcmp(xml_bits{1}{2},['/' start_tag])
        if dataholder
            datastring = xml_bits{1}{1}(2:end);
            if (~isempty(regexp(datastring,'[^0-9eE\[\]:;,\-\.\+ infaINFA]','once')) ...
                    ...&& ...
                    ...~strcmpi(datastring,'inf') && ...
                    ...~strcmpi(datastring,'nan') ...
                    ) ...|| ...
                    %isempty(str2num(datastring))
                xml = datastring;
            else
                xml = str2num(datastring);
                if isempty(xml)
                    xml = datastring;
                end;
            end;
        end;
        restofstring = xml_bits{1}{3};
        return;
    else
        dataholder = false;
        [newdata current_content] = XMLParser(xml_bits{1}{3}, xml_bits{1}{2});
        
        if ~exist('xml','var') || ~isfield(xml, xml_bits{1}{2})
            xml.(xml_bits{1}{2}) = newdata;
        else
            
            if isstruct(xml.(xml_bits{1}{2}))
                olddata = xml.(xml_bits{1}{2});
                if isstruct(newdata)
                    newfields = fieldnames(newdata);
                    for f=1:length(newfields)
                        xml.(xml_bits{1}{2})(length(olddata)+1).(newfields{f}) = newdata.(newfields{f});
                    end;
                else
                    xml.(xml_bits{1}{2}) = {olddata(1)};
                    for n=2:length(olddata)
                        xml.(xml_bits{1}{2}){n} = olddata(n);
                    end;
                end;
                    
            elseif iscell(xml.(xml_bits{1}{2}))
                xml.(xml_bits{1}{2}){end+1} = newdata;
            else
                xml.(xml_bits{1}{2}) = {xml.(xml_bits{1}{2})};
                xml.(xml_bits{1}{2}){end+1} = newdata;
            end;
        
        end;
    end;

end;
