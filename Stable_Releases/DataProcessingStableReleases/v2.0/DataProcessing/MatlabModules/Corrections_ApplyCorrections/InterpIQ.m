function [iq_value] = InterpIQ(dataset_name, iq_times, iq_values)
%
%
%  _____           _                          _____   ___      
% |_   _|         / |_                       |_   _|.'   `.    
%   | |   _ .--. `| |-'.---.  _ .--.  _ .--.   | | /  .-.  \   
%   | |  [ `.-. | | | / /__\\[ `/'`\][ '/'`\ \ | | | |   | |   
%  _| |_  | | | | | |,| \__., | |     | \__/ |_| |_\  `-'  \_  
% |_____|[___||__]\__/ '.__.'[___]    | ;.__/|_____|`.___.\__| 
%                                   [__|                       
%
%
%
%  function [iq_value] = InterpIQ(dataset_name, iq_times, iq_values)
%
%  This function interpolates between the available IQs given the times of
%  the available IQs and the current dataset name.
%
%
%  Inputs
%     dataset_name  - the name of the dataset (i.e. from d.admin.filename_prefix).
%     
%     iq_times - an array containing the unix time of the available IQs.
%
%     iq_values - an array containing the IQ values to be interpolated between.
%
%  Outputs
%
%     iq_value - the interpolated value of the IQ.
%
%  Created BE 06/03/2013


  [a b] = size(iq_times);
  for i=1:b
    id(i) = i;
  end

  [first_iq_time first_iq_index] = min(iq_times);
  [last_iq_time last_iq_index] = max(iq_times);
  
  this_dataset_time = filename2epoch_framework(dataset_name);
  if (this_dataset_time<first_iq_time)
    iq_value = iq_values(first_iq_index); 
  elseif (this_dataset_time>last_iq_time)
    iq_value = iq_values(last_iq_index);
  else
    iq_value = interp1(iq_times,iq_values,this_dataset_time,'linear'); 
  end
 
  if b==0
    fprintf('Found no acceptable IQs.....problem');
  end
  