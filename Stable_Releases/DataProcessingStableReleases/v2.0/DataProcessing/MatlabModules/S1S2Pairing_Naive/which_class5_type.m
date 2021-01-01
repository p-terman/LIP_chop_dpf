function dp = which_class5_type (dp)

% this takes a dp with rq information and determines whether or not a given
% else/class 5 pulse is s1 like or s2 like
% an s1 like pulse is where the class 5 is at the beginning of the pulse
% (before an s2) and there is no previously classified s1.
% an s2 like pulse is where there has been a previously classified s2 and
% there are class 5 pulses that follow it. 
% No pulse is given these designations unless it is at least 10 phe.
%
%   update 20190220 PAT - change output to uint8 

%initialize the new matricies
s1_like_class5 = false(size(dp.pulse_classification ,1), size(dp.pulse_classification, 2), 'logical');
s2_like_class5 = false(size(dp.pulse_classification ,1), size(dp.pulse_classification, 2), 'logical');

for ii = 1: size(dp.pulse_classification , 2)
    
    % find the indices for S1 and S2 pulses
    s1_inds = find(dp.pulse_classification(:,ii) == 1);
    s2_inds = find(dp.pulse_classification(:,ii) == 2); % | dp.pulse_classification(:,ii) == 4); % treating all SE as s2
    class5_inds = find(dp.pulse_classification(:,ii) == 5);
    
    % how many pulses of each type are there?
    num_s1 = length(s1_inds);
    num_s2 = length(s2_inds);               
    num_class5 = length(class5_inds);
    
    
    % are there one or more pulses of each type?
    has_s1 = num_s1 > 0; 
    has_s2 = num_s2 > 0;
    has_class5 = num_class5 > 0;

    
    if has_s1 %s1 is found  - the any else pulses cannot be s1 like, but might be s2_like
      
        if has_s2 %s2 is found
            s2_first_start = dp.pulse_start_samples(s2_inds(1), ii); % the first s2 start time
            
           for j = 1:num_class5 %looping over all class 5
              k=class5_inds(j); %k is index in class 5 space
              if dp.pulse_start_samples(k,ii)> s2_first_start %then this is after the first s2 and we will treat as s2 like class5 
                  s2_like_class5 (k, ii) = 1;

              end
           end %for loop 
        end %if has  s2
    end %if has_s1
    
    
    if ~has_s1 %no s1 found, so we will look for an 's1-like' class 5 in the first class 5 pulse- or see if it is 's2-like'
        if has_class5 && has_s2
            s2_first_start = dp.pulse_start_samples(s2_inds(1), ii); % the first s2 start time
            class5s_before_s2 = find(dp.pulse_start_samples(class5_inds, ii) < s2_first_start);%corrected with find 180622
            %find all class 5 pulses before first s2 and then take the
            %largest as the s1
            class5s_before_s2_ind = class5_inds(class5s_before_s2);
            [~ , largest_class5] = max(dp.pulse_area_phe(class5s_before_s2_ind, ii));
            largest_class5_ind= class5_inds(largest_class5);
            s1_like_class5(largest_class5_ind, ii) = 1;
           
           %the above dealt with first class 5 index. Next look at other indicies 
           %check remaining pulses
           for j = largest_class5_ind:num_class5
              k=class5_inds(j); %k is index in class 5 space
              %s1_like_class5(k,ii)= 0; % we know it is not s1-like if it is past first class 5 
              if dp.pulse_start_samples(k,ii)> s2_first_start %then this is after the first s2 and we will treat as s2 like class5 
                  s2_like_class5 (k, ii) = 1;
              end
           end
            
        end % has s2 and class 5
        
        %now what to do when there is a class 5 but no s2
        if ~has_s2 && has_class5 % also remember no s1
            [~ , largest_class5] = max(dp.pulse_area_phe(class5_inds, ii));
            largest_class5_ind= class5_inds(largest_class5);
            s1_like_class5(largest_class5_ind, ii) = 1;
            
 
            %}
        end    
    end % if doesn't have s1
  
end %main for loop

%removed due to redundancy 
% a= dp.pulse_classification~=5; %where there is no class 5 pulse at all
% s1_like_class5(a) = 0;
% s2_like_class5(a) = 0;
 
%there is a problem where a small class 5 pulse will get -100 position flag
%and then carry this over to the corrections, ruining the event. 
a = find(dp.pulse_area_phe < 10);
s1_like_class5(a) = 0;
s2_like_class5(a) = 0;

%save them
dp.s1_like_class5 = uint8(s1_like_class5);
dp.s2_like_class5 = uint8(s2_like_class5);

end


