function [PMT_Radius PMTG] = LUXLRF_GroupDefinition


% [PMT_Radius PMTG] = LUXLRF_GroupDefinition
% 
% Groups de top array according to the radius
% 
% Inputs:
%      None
%
% Outputs:
%  
%  PMT_Radius - the radius of each PMT rounded to cm
%  PMTG - a vector with the radius of each group
%  

if nargout <1
    fprintf('You should try the following\n');
    fprintf('[PMT_Radius PMTG] = LUXLRF_GroupDefinition\n');
end



[pmt_pos sextant_arrangement] = LUXPMTArray;
for i=1:122,
    PMT_r(i,:)=pmt_pos(find(sextant_arrangement==i),:);
end

PMT_Radius = round(sqrt(PMT_r(1:122,2).^2 + PMT_r(1:122,1).^2)-PMT_r(1:122,3)+24.4000);
PMTG = sort(unique(PMT_Radius([1:60 121])'), 'Ascend'); 

