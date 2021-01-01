% function xmlstring = MakeXMLString_framework(xmlstruct,indentlevel)
%
% MakeXMLString is a recursive function that makes an XMLString from a
% structure.  Called by XMLWriter
%
% 10/05/2005 ED
% Release Version 1.0

function [xmlstring success] = MakeXMLString_framework(xmlstruct,indentlevel)

%% default indentlevel is zero
if nargin<2
    indentlevel = 0;
end;

%% error handlers
if ~isstruct(xmlstruct)
    disp('passed a non-struct to MakeXMLString');
    success = false;
    xmlstring = '';
    return;
end;

if prod(size(xmlstruct)) ~= 1
    disp('passed a non-scalar struct to MakeXMLString');
    success = false;
    xmlstring = '';
    return;
end;

%% create line indent string
indent_string = '';
for n=1:indentlevel
    indent_string = [indent_string '    '];
end;

%% file the fields in xmlstruct
fields = fieldnames(xmlstruct);

%% initialize xmlstring
xmlstring = '';

%% loop through fields recursively
for f=1:length(fields)
    if isstruct(xmlstruct.(fields{f}))
        if prod(size(xmlstruct.(fields{f}))) == length(xmlstruct.(fields{f}))
            for s=1:length(xmlstruct.(fields{f}))
                xmlstring = [xmlstring indent_string sprintf('<%s>\n',fields{f})];
                [subxmlstring subsuccess] = MakeXMLString_framework(xmlstruct.(fields{f})(s),indentlevel+1);
                if subsuccess
                    xmlstring = [xmlstring subxmlstring indent_string sprintf('</%s>\n',fields{f})];
                else
                    xmlstring = '';
                    success = false;
                    return;
                end;
            end
        else
            disp('struct passed to MakeXMLString contains multi-dimensional struct array');
            xmlstring = '';
            success = false;
            return;
        end;
    elseif iscell(xmlstruct.(fields{f}))
        if prod(size(xmlstruct.(fields{f}))) == length(xmlstruct.(fields{f}))
            for s=1:length(xmlstruct.(fields{f}))
                if isnumeric(xmlstruct.(fields{f}){s}) || islogical(xmlstruct.(fields{f}){s})
                    datastring = num2str(xmlstruct.(fields{f}){s});
                    size_datastring = size(datastring);
                    xmldatastring = '';
                    for r=1:size_datastring(1)
                        xmldatastring = [xmldatastring '; ' datastring(r,:)];
                    end;
                    xmldatastring = xmldatastring(3:end);
                else
                    xmldatastring = char(xmlstruct.(fields{f}){s});
                end;
                xmlstring = [xmlstring indent_string sprintf('<%s>%s</%s>\n',fields{f},xmldatastring,fields{f})];
            end;
        else
            disp('struct passed to MakeXMLString contains multi-dimension cell array');
            xmlstring = '';
            success = false;
            return;
        end;
    else
        if isnumeric(xmlstruct.(fields{f})) || islogical(xmlstruct.(fields{f}))
            datastring = num2str(xmlstruct.(fields{f}));
            size_datastring = size(datastring);
            xmldatastring = '';
            for r=1:size_datastring(1)
                xmldatastring = [xmldatastring '; ' datastring(r,:)];
            end;
            xmldatastring = xmldatastring(3:end);
        else
            xmldatastring = char(xmlstruct.(fields{f}));
        end;
        xmlstring = [xmlstring indent_string sprintf('<%s>%s</%s>\n',fields{f},xmldatastring,fields{f})];
    end;
end;
success = true;