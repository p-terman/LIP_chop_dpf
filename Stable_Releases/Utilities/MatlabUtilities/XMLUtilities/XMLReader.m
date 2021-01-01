% function xml = XMLReader(filename)
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

function xml = XMLReader(filename)


%% open xml file and load its content into the variable filecontent
xmlfid = fopen(filename,'r');
filecontent = char(fread(xmlfid,inf,'uchar')');
fclose(xmlfid);

%% Parse xml file

xml = XMLParser(filecontent);
