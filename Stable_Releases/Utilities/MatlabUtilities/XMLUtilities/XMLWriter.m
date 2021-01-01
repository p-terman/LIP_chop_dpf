% function success = XMLWriter(filename,xmlstruct)
%
% XMLWriter takes a structure and constructs an XML file out of it.
% It does not print a first line identifying the file as an xml file,
% 'cause I don't know what the info in that line means -- if you like, you
% can add that into the code below -- there's a spot for it.
%
% 10/05/2005 ED
% Release Version 1.0
% 2010-04012 JJC
% Took out firstline, return length of xmlstring 

function [success errormessage xmlstringlength] = XMLWriter(filename,xmlstruct)

%% set defaults
success = false;
errormessage = '';

%% generate xmlstring
[xmlstring success] = MakeXMLString(xmlstruct);

if ~success
    errormessage = 'Error making xml string in XMLWriter';
    return;
end;

%% write xml file

xmlfid = fopen(filename,'w');
if xmlfid<0
    errormessage = 'Failed to open xmlfile for writing in XMLWriter';
    fclose(xmlfid);
    return;
end;
%firstline = '<?xml version="1.0"?>';
%fprintf(xmlfid,'%s\n',firstline);
fprintf(xmlfid,'%s\n',xmlstring);
fclose(xmlfid);
success = true;
xmlstringlength = length(xmlstring);