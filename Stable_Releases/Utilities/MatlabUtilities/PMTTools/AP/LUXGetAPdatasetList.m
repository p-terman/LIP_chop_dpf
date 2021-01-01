function [ APdatasetList commentlist] = LUXGetAPdatasetList( directory )
%
% function [ APdatasetList commentlist] = LUXGetAPdatasetList( directory )
%
% Returns a list of the names of files in a directory
%
% Inputs:
%
%           directory - path of diretory
%
% Outputs: 
%
%       APdatasetList - list of datasets
%
% JRV - 20091028


d = dir(directory);

[sorted_names, sorted_index] = sortrows({d.name}');
APdatasetList = sorted_names;

if length(sorted_index) > 0
    commentlist(length(sorted_index)) = ' ';
else
    commentlist = [];
end