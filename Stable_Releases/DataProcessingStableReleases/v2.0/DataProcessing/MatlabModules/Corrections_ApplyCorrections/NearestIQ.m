function [index iq_time] = NearestIQ(dataset_name, iq_times)
%
%
%     _   __                          __  ________ 
%    / | / /__  ____ _________  _____/ /_/  _/ __ \
%   /  |/ / _ \/ __ `/ ___/ _ \/ ___/ __// // / / /
%  / /|  /  __/ /_/ / /  /  __(__  ) /__/ // /_/ / 
% /_/ |_/\___/\__,_/_/   \___/____/\__/___/\___\_\ 
%
%
%
%  function [index iq_time] = NearestIQ(dataset_name, iq_times)
%
%  This function finds the closest available IQ given the times of
%  the available IQs and the current dataset name.
%
%
%  Inputs
%     dataset_name  - the name of the dataset (i.e. from d.admin.filename_prefix).
%     
%     iq_times - an array containing the unix time of the available IQs.
%
%
%  Outputs
%     
%     index - the index (from the iq_times array) of the closest time.
%
%     iq_time - the IQ time closest to that of the current dataset.
%
%     iq_number - the LUG IQ entry of the IQ chosen.
%
%  Created BE 06/03/2013


  [a b] = size(iq_times);
  for i=1:b
    id(i) = i;
  end

  [first_iq_time first_iq_index] = min(iq_times);
  [last_iq_time last_iq_index] = max(iq_times);
  
  this_dataset_time = filename2epoch_framework(dataset_name);

  if b>1 
    index = interp1(iq_times,id,this_dataset_time,'nearest');

    if (this_dataset_time<first_iq_time)
      index = first_iq_index; 
 %     fprintf('This dataset is from before the first calculated correction. \n');
 %     fprintf('Using the first IQ value. \n');
    end    
    if (this_dataset_time>last_iq_time)
      index = last_iq_index;  
  %    fprintf('This dataset is from after the last calculated correction. \n');
  %    fprintf('Using the last IQ value. \n');
    end
    
  elseif b==1
    index = 1;
 %   fprintf('Only one acceptable IQ found - using that!');

  elseif b==0
    fprintf('Found no acceptable IQs.....problem');
  end
    
    

  
  iq_time = iq_times(index);
  
  
  