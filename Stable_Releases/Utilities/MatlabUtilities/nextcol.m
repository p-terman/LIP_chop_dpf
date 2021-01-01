function [str] = nextcol(varargin);
% Gives next colour in cyclic list
% Provide a cell list if want to set colors in list
% rjg Mar99,Oct99
% nextcol(-2) show current list
% nextcol(-1) reinit default
% nextcol(n>=0) set to this part of list

global nextcol_list nextcol_index nextcol_readme
% ... global brings variables into existance

% If globals don't yet exist then set them
   if( isempty(nextcol_list) | isempty(nextcol_index) )
      % DEFAULT VALUES 
disp('nextcol: Default colour list set');
      %nextcol_list = {'r+','yx','go','bs'}; % ,'m','c','k'};
      %nextcol_list = {'r','m','g','b','y','c','k'};
      nextcol_list = {'r','m','g','b','c','k'};
      %nextcol_list
      nextcol_index = 0;
      nextcol_readme = ...
'Function nextcol provides an endless cyling list of colours, see help nextcol';
   end


 if( length(varargin)==0 )
   % Provide next element in list
   str = nextcol_list{nextcol_index+1};   % +1 because list starts at one
   nextcol_index = mod(nextcol_index+1,length(nextcol_list));

 else
   if( isa(varargin{1},'double') )
     if(varargin{1}==-1)
       clear global nextcol_index; 
       nextcol; % Causes default to be reset
       str = nextcol_list;
       return; 
     elseif(varargin{1}==-2) 
       str = nextcol_list;
       return; 
     else
       % Use number to set position in list
       nextcol_index = mod(varargin{1},length(nextcol_list));
     end
   else
     nextcol_list = varargin{1}; % First argument
     nextcol_index = 0;
   end
   % return an argument, in case it is needed, however, notice that we have
   % not incremented index !! so we will get same color next time we call
   % without an argument
   str = nextcol_list{nextcol_index+1};   % +1 because list starts at one
 end

















