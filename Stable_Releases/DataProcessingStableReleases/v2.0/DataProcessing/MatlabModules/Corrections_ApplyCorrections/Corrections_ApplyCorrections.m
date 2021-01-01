function status = ApplyCorrections(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
% status = ApplyCorrections_Basic(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%
%
%
%
%
%
% Required RQs:
lifetime_values = [];
dataset_times = [];
filename_prefixs = {};
cp_number = {};

%
%
%
%
%
% Versioning:
%   20130306 BE - Created
%   20180131 PAT - allowed for class 4 and 5 pulses to be counted as s2 and
%   s1 respectively. 
%
%
%
%
% RQ versions:
%
%
%
%  Required RQs:
%
%  - event_number, 
%
%  
%


%
%% Load .rq file
if strcmp(filename_rq(1),'.') % check for '.' infornt of the file name
    filename_prefix=filename_rq(2:20);
else
    filename_prefix=filename_rq(1:19);
end

dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);

settings.evt_settings = dp.admin.evt_settings;
settings.daq_settings = dp.admin.daq_settings;
settings.filename_prefix = dp.admin.filename_prefix;

%% Bookkeeping

myname = 'ApplyCorrections';
fprintf('\n\n *** Starting module %s\n',myname);

fprintf('\n Data processing xml file = %s\n',data_processing_xml_path);

dp_settings_xml = XMLReader_framework(data_processing_xml_path);
lug_iqs_xml = XMLReader_framework(iq_xml_path);

module_names = {dp_settings_xml.data_processing_settings.module.module_name};
index_temp = strfind(module_names,myname);
index_module = find(not(cellfun('isempty', index_temp)));

mymodule_settings = dp_settings_xml.data_processing_settings.module(index_module).parameters;

max_num_pulses = dp_settings_xml.data_processing_settings.global.max_num_pulses;


%% Initialize variables

pmt_chs = 1:122; % must take from daq settings

N = length(dp.event_number);


%dp = which_class5_type(dp);  %180205 PAT added 

pulse_event_size = [max_num_pulses N];

%livetime = LUXGetLivetime_framework(filename_evt,data_path_evt);
livetime = dp.admin.livetime;

%  Correcting possible issue with top_bottom_ratio

[a1 b1] = size(dp.pulse_area_phe);  
[a2 b2] = size(dp.top_bottom_ratio);  adif = a1-a2;

blanks = [];

for i=1:b1
   blanks = [blanks 0]; 
end

if adif>0
   for i=1:adif
     dp.top_bottom_ratio = vertcat(dp.top_bottom_ratio,blanks);    
   end
end

% Initializing corrected values of area

dp.z_corrected_pulse_area_all_phe = dp.pulse_area_phe;
dp.xyz_corrected_pulse_area_all_phe = dp.pulse_area_phe;

dp.z_corrected_pulse_area_bot_phe = dp.pulse_area_phe./(1+dp.top_bottom_ratio);
dp.xyz_corrected_pulse_area_bot_phe = dp.pulse_area_phe./(1+dp.top_bottom_ratio);

dp.correction_electron_lifetime = zeros(max_num_pulses,N, 'single');
dp.correction_s1_z_dependence = zeros(max_num_pulses,N, 'single');
dp.correction_s1_xyz_dependence = zeros(max_num_pulses,N, 'single');
dp.correction_s2_xy_dependence = zeros(max_num_pulses,N, 'single');

s1_ref_z_ns = mymodule_settings.detector_centre_z_ns;
allowed_gs = mymodule_settings.allowed_gs;


%-----------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%% Reading iq values from the IQs xml ----------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------- 
% electron_lifetime - Interpolate - DONE
% z_dep_par_all - Nearest - DONE
% z_dep_par_bot - Nearest - DONE
% x_bins,y_bins,z_bins,s1_map_all - Nearest
% x_bins,y_bins,z_bins,s1_map_bot - Nearest
% x_bins,y_bins,s2_map_all - Nearest
% x_bins,y_bins,s2_map_bot - Nearest

[a num_iqs] = size(lug_iqs_xml.iq);

%-----------------------------------
% Finding the electron lifetime
%-----------------------------------

lifetime_values = [];
dataset_times = [];
filename_prefixs = {};
cp_number = {};

 for i=1:num_iqs
  if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

    if (isfield(lug_iqs_xml.iq(i).correction.fit,'electron_attenuation_us').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
 %      if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
 %          fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
           lifetime_values = [lifetime_values lug_iqs_xml.iq(i).correction.fit.electron_attenuation_us];
           dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix) ];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
 %      else
 %          fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
 %      end
              
    end
  end  
 end
 
current_data_set_time=filename2epoch_framework(filename_prefix );
 
 if inrange(current_data_set_time,[min(dataset_times), max(dataset_times)] )
     [electron_lifetime] = InterpIQ(filename_prefix,dataset_times,lifetime_values);
 else    
     [index electron_lifetime_time] = NearestIQ(filename_prefix,dataset_times);
     electron_lifetime = lifetime_values(index);
 end

%{  
/// Commented Out, maybe use Z bins for S2 correction in the future -AD

%------------------------------------------
% Finding the S2 Z correction map values in each Z bin
%------------------------------------------

s2_all_z0_vals = [];
s2_bot_z0_vals = [];
dataset_times = [];
filename_prefixs = {};
cp_number = {};

 for i=1:num_iqs
  if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

    if (isfield(lug_iqs_xml.iq(i).correction.fit,'s2_bottom_z0_phe').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
%       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
%           fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
           s2_all_z0_vals = [s2_both_z0 lug_iqs_xml.iq(i).correction.fit.s2_both_z0_phe];
           s2_bot_z0_vals = [s2_both_z0 lug_iqs_xml.iq(i).correction.fit.s2_bottom_z0_phe];
           dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
 %      else
 %          fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
 %      end
              
    end
  end  
 end
 
current_data_set_time=filename2epoch_framework(filename_prefix );
 
 if inrange(current_data_set_time,[min(dataset_times), max(dataset_times)] )
     [s2_all_z0] = InterpIQ(filename_prefix,dataset_times,s2_all_z0_vals);
     [s2_bot_z0] = InterpIQ(filename_prefix,dataset_times,s2_bot_z0_vals);
 else    
     [index z0_times] = NearestIQ(filename_prefix,dataset_times); 
     s2_all_z0 = s2_all_z0_vals(index);
     s2_bot_z0 = s2_bot_z0_vals(index);
 end
   

s2_z_index = [];
dataset_times = [];
filename_prefixs = {};
cp_number = {};

 for i=1:num_iqs
  if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

    if (isfield(lug_iqs_xml.iq(i).correction.fit,'bin_means_s2_bottom').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
%       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
%           fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
           s2_z_index = [s2_z_index i];
           dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
 %      else
 %          fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
 %      end
              
    end
  end  
 end
 
   [index iq_time_s2zDep] = NearestIQ(filename_prefix,dataset_times);
      
   s2_z_bins = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.bin_means_s2_bottom;
   s2_z_all = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.means_s2_bottom;
   s2_z_bottom = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.means_s2_both;
   
 %}
 
 
%-----------------------------------
% Finding the S1 Z-dep values
%-----------------------------------

z_dep_both_values = zeros(0,3);
z_dep_bottom_values = zeros(0,3);
dataset_times = [];
filename_prefixs = {};
cp_number = {};

 for i=1:num_iqs
  if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

    if (isfield(lug_iqs_xml.iq(i).correction.fit,'s1_both_zdep_quad_fit').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
%       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
%           fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
           z_dep_both_values = vertcat(z_dep_both_values,lug_iqs_xml.iq(i).correction.fit.s1_both_zdep_quad_fit);
           z_dep_bottom_values = vertcat(z_dep_bottom_values,lug_iqs_xml.iq(i).correction.fit.s1_bottom_zdep_quad_fit);
           
           dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix )];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
%       else
%           fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
%       end
              
    end
  end  
 end
 
   [index iq_time_zDep] = NearestIQ(filename_prefix,dataset_times);
   
   z_dep_par_all = z_dep_both_values(index,:);
   z_dep_par_bot = z_dep_bottom_values(index,:);
      

%------------------------------------------
% Finding the S2 xy correction map values
%------------------------------------------

s2_xy_index = [];
dataset_times = [];
filename_prefixs = {};
cp_number = {};

 for i=1:num_iqs
  if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

    if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s2_both').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1        
%       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
%           fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
           s2_xy_index = [s2_xy_index i];
           dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
 %      else
 %          fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
 %      end
              
    end
  end  
 end
 
   [index iq_time_s2xyDep] = NearestIQ(filename_prefix,dataset_times);
      
   s2_x_bins = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.x_bin_center;
   s2_y_bins = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.y_bin_center;
   s2_map_all = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.norm_s2_both;
   s2_map_bottom = lug_iqs_xml.iq(s2_xy_index(index)).correction.fit.norm_s2_bottom;
   
%------------------------------------------
% Finding the S1 xy correction map values
%------------------------------------------
   
s1_xy_index = [];
dataset_times = [];
filename_prefixs = {};
cp_number = {};

 for i=1:num_iqs
  if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

    if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s1_both').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1       
%       if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
 %          fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);
           s1_xy_index = [s1_xy_index i];                      
           dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix ) ];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
 %      else
 %          fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
 %      end
              
    end
  end  
 end
 
   [index iq_time] = NearestIQ(filename_prefix,dataset_times);
      
   s1_x_bins = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.x_bin_center;
   s1_y_bins = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.y_bin_center;
   s1_map_all = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.norm_s1_both;
   s1_map_bottom = lug_iqs_xml.iq(s1_xy_index(index)).correction.fit.norm_s1_bottom;

%------------------------------------------
% Finding the S1 xyz correction map values
%------------------------------------------
   
s1_xyz_index = [];
dataset_times = [];
filename_prefixs = {};
cp_number = {};

 for i=1:num_iqs
  if isfield(lug_iqs_xml.iq(i).correction,'fit')==1   

    if (isfield(lug_iqs_xml.iq(i).correction.fit,'norm_s1_both_xyz').*strncmp(lug_iqs_xml.iq(i).global.algorithm_name,'LUXkrypCal',10))==1       
 %      if ismember(lug_iqs_xml.iq(i).global.gs_number,allowed_gs)==1;
 %          fprintf('Its allowed %d  - %s \n',lug_iqs_xml.iq(i).global.gs_number,lug_iqs_xml.iq(i).global.filename_prefix);

           s1_xyz_index = [s1_xyz_index i];                      
           dataset_times = [dataset_times filename2epoch_framework(lug_iqs_xml.iq(i).global.filename_prefix )];
           filename_prefixs = [filename_prefixs lug_iqs_xml.iq(i).global.filename_prefix];
           cp_number = [cp_number lug_iqs_xml.iq(i).global.cp_number];
%       else
%           fprintf('Its not allowed %d \n',lug_iqs_xml.iq(i).global.gs_number);
%       end
              
    end
  end  
 end
 
   [index iq_time_s1xyzDep] = NearestIQ(filename_prefix,dataset_times);
      
   s1_xyz_x_bins = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.x_bin_center;
   s1_xyz_y_bins = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.y_bin_center;
   s1_xyz_z_bins = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.z_bin_center;
   
   xx = size(s1_xyz_x_bins); yy = size(s1_xyz_y_bins); zz = size(s1_xyz_z_bins);
   
   s1_xyz_map_all = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_both_xyz; s1_xyz_map_all = reshape(s1_xyz_map_all,xx(2),yy(2),zz(2));
   s1_xyz_map_bottom = lug_iqs_xml.iq(s1_xyz_index(index)).correction.fit.norm_s1_bottom_xyz;   s1_xyz_map_bottom = reshape(s1_xyz_map_bottom,xx(2),yy(2),zz(2));
    


%-----------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%% Calculating corrections ----------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------- 

%% Finding the drift time, x and y associated with S1

% Finding the drift time and xy position from the largest S2 pulse in the
% pairing


[a b] = size(dp.z_drift_samples); % b is number of events in file ??

drift = dp.z_drift_samples;
drift(find(isnan(drift))) = 0.0;
s1_drift_ns = +(dp.pulse_classification==1 | dp.s1_like_class5); % changed from only class 1 PAT 180131 - 180205 changed to have s1 and d2 like class 5
s1_x_cm = +(dp.pulse_classification==1 | dp.s1_like_class5); % same
s1_y_cm = +(dp.pulse_classification==1 | dp.s1_like_class5); %same
s2_phe =  dp.pulse_area_phe.*(dp.pulse_classification==2 | dp.pulse_classification==4 | dp.s2_like_class5); % changed from only class 2
s2_phe(find(isnan(s2_phe))) = 0.0;

%s1s = (sum(dp.pulse_classification==1)==1) & (sum(dp.pulse_classification==2)>0) & (sum(dp.s1s2_pairing==1)>1);

 s1s = (sum(dp.s1s2_pairing==1)>1); % tells you if its an event with at least one s1 and one s2.
  
 for i=1:length(s1s)  % for each event in the file
     if s1s(i)>0 % if the event had an S1
       if s1s(i)==1
           %... pair first S1 with first S2 -commented out by AD, July 21
           %2014
         %[r c v] = find(drift(:,i),1,'first');
         %[c1 r1 v1] = find(dp.pulse_classification(:,i)==1,1,'first');
         %[c2 r2 v2] = find(dp.pulse_classification(:,i)==2,1,'first');
         %s1_drift_ns(c1,i) = 10.*v;
         %s1_x_cm(c1,i) = dp.x_cm(c2,i);
         %s1_y_cm(c1,i) = dp.y_cm(c2,i);   
        
         [v r] = max(s2_phe(:,i));% value and index (r for row) of Max S2
         [c1 r1 v1] = find(dp.pulse_classification(:,i)==1 | dp.s1_like_class5(:,i),1,'first');
         %this line (regarding s1) removed by PAT 180131 -r1 and v1
         %appeared to be unused/dummy vars so I did not reproduce their
         %info - line added back and modified 180205
         %{
           removed 180205
         s1_inds = dp.pulse_classification==1  ;% allows for s1 or class 5 pulses
         [~, c1] = max(dp.pulse_area_phe(s1_inds(:,i), i)) ; %the ~ is just a dummy var, we want the index - added by PAT to take the largest s1 as the 'first'
         %}
         %[c2 r2 v2] = find(dp.pulse_classification(:,i)==2,1,'first');
         s1_drift_ns(c1,i) = 10.*drift(r,i); %row c1,r and column i
         s1_x_cm(c1,i) = dp.x_cm(r,i); %it seems that this is taking the location of the max s2 pulse and taking that to be the location of the s1
         %this can make sense since caluclated s1 pos is useless, so we
         %need a basis in the s2 to have a pos for s1
         s1_y_cm(c1,i) = dp.y_cm(r,i);  
         
       else %no s1 or class 5
          s1_drift_ns(:,i) = 0;
          s1_x_cm(:,i) = 0;
          s1_y_cm(:,i) = 0;
       end  
     end
 end

 
 
fprintf('Done finding S1 depths :) \n');

s1_drift_ns(find((s1_drift_ns==0))) = nan;
s1_x_cm(find((s1_x_cm==0))) = nan;
s1_y_cm(find((s1_y_cm==0))) = nan;



%--------------------------------------------------------------------------
% Calculate Z corrections
%--------------------------------------------------------------------------

% Calculate S1 Z-correction

s1_z_correction = polyval(z_dep_par_all,s1_ref_z_ns./1000)./polyval(z_dep_par_all,s1_drift_ns./1000);
s1_z_correction(find(isnan(s1_z_correction))) = 1.0;

s1_z_correction_bot = polyval(z_dep_par_bot,s1_ref_z_ns./1000)./polyval(z_dep_par_bot,s1_drift_ns./1000);
s1_z_correction_bot(find(isnan(s1_z_correction_bot))) = 1.0;

% Calculate electron lifetime correction (S2 Z-correction)
% Reading the values of electron lifetime from the LUG and interpolating
% between the two either side

electron_lifetime_correction = exp(double(dp.z_drift_samples)./(100.*electron_lifetime));
electron_lifetime_correction(find(isnan(electron_lifetime_correction))) = 1.0;
 
%--------------------------------------------------------------------------
% Calculate XYZ corrections
%--------------------------------------------------------------------------

% Calculate S1 xyz corrections from map stored in IQ

%s1_x_cm(isnan(s1_x_cm)) = 0; s1_y_cm(isnan(s1_y_cm)) = 0; s1_drift_ns(isnan(s1_drift_ns)) = 0;

s1xyz_correction = interp3(s1_xyz_x_bins,s1_xyz_y_bins,s1_xyz_z_bins,s1_xyz_map_all,s1_x_cm,s1_y_cm,s1_drift_ns./1000,'spline');
s1xyz_correction(find(isnan(s1xyz_correction))) = 1.0;

s1xyz_correction_bot = interp3(s1_xyz_x_bins,s1_xyz_y_bins,s1_xyz_z_bins,s1_xyz_map_bottom,s1_x_cm,s1_y_cm,s1_drift_ns./1000,'spline');
s1xyz_correction_bot(find(isnan(s1xyz_correction_bot))) = 1.0;


% Calculate S1 XY-only corrections

% s2x = dp.x_cm.*(+(dp.pulse_classification==2)); s2x(find(s2x==0)) = nan;
% s2y = dp.y_cm.*(+(dp.pulse_classification==2)); s2y(find(s2y==0)) = nan;
% 
% s2xy_correction = interp2(s1_x_bins,s1_y_bins,s2_map_all,dp.x_cm,dp.y_cm,'spline');
% s2xy_correction = s2xy_correction.*(+dp.pulse_classification==2);
% s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;
% 
% s2xy_correction_bot = interp2(s1_x_bins,s1_y_bins,s2_map_bottom,dp.x_cm,dp.y_cm,'spline');
% s2xy_correction_bot = s2xy_correction_bot.*(+dp.pulse_classification==2);
% s2xy_correction_bot(find(s2xy_correction_bot==0))=1.0;  s2xy_correction_bot(find(isnan(s2xy_correction_bot)))=1.0;




% Calculate S2 XY corrections

s2x = dp.x_cm.*(+(dp.pulse_classification==2 | dp.pulse_classification==4 | dp.s2_like_class5)); s2x(find(s2x==0)) = nan; %PAT changed to allow for class 4 to be counted as s2, and also s2 like class 5
s2y = dp.y_cm.*(+(dp.pulse_classification==2 | dp.pulse_classification==4 | dp.s2_like_class5)); s2y(find(s2y==0)) = nan; %180131 same - these seem to not be used??

s2xy_correction = interp2(s2_x_bins,s2_y_bins,s2_map_all,dp.x_cm,dp.y_cm,'spline');
s2xy_correction = s2xy_correction.*(+(dp.pulse_classification==2 | dp.pulse_classification==4 | dp.s2_like_class5)); % PAT changed
s2xy_correction(find(s2xy_correction==0))=1.0;  s2xy_correction(find(isnan(s2xy_correction)))=1.0;

s2xy_correction_bot = interp2(s2_x_bins,s2_y_bins,s2_map_bottom,dp.x_cm,dp.y_cm,'spline');
s2xy_correction_bot = s2xy_correction_bot.*(+(dp.pulse_classification==2 | dp.pulse_classification==4 | dp.s2_like_class5));% PAT modified
s2xy_correction_bot(find(s2xy_correction_bot==0))=1.0;  s2xy_correction_bot(find(isnan(s2xy_correction_bot)))=1.0;



% add RQs of correction factors

dp.correction_electron_lifetime = single(electron_lifetime_correction);
dp.correction_s1_z_dependence = single(s1_z_correction);
dp.correction_s1_z_dependence_bot = single(s1_z_correction_bot);
dp.correction_s1_xyz_dependence = single(s1xyz_correction);
dp.correction_s1_xyz_dependence_bot = single(s1xyz_correction_bot);
dp.correction_s2_xy_dependence = single(s2xy_correction);
dp.correction_s2_xy_dependence_bot = single(s2xy_correction_bot);

dp.admin.corrections.electron_lifetime = electron_lifetime;
dp.admin.corrections.s1_z_dependence.iq_time = iq_time_zDep;
dp.admin.corrections.s1_z_dependence.all = z_dep_par_all;
dp.admin.corrections.s1_z_dependence.bottom = z_dep_par_bot;
dp.admin.corrections.s1_xyz_dependence.iq_time = iq_time_s1xyzDep;
dp.admin.corrections.s2_xy_dependence.iq_time = iq_time_s2xyDep;


%-----------------------------------------------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%% Applying corrections ----------------------------------------------------------------------------
%-----------------------------------------------------------------------------------------------------------------
%----------------------------------------------------------------------------------------------------------------- 

%--------------------------------------------------------------------------
% Apply Z corrections
%--------------------------------------------------------------------------

dp.z_corrected_pulse_area_all_phe = dp.z_corrected_pulse_area_all_phe.*electron_lifetime_correction.*s1_z_correction; 
dp.z_corrected_pulse_area_bot_phe = dp.z_corrected_pulse_area_bot_phe.*electron_lifetime_correction.*s1_z_correction_bot; 

%--------------------------------------------------------------------------
% Apply XYZ corrections
%--------------------------------------------------------------------------

dp.xyz_corrected_pulse_area_all_phe = dp.xyz_corrected_pulse_area_all_phe.*electron_lifetime_correction.*s2xy_correction.*s1xyz_correction;
dp.xyz_corrected_pulse_area_bot_phe = dp.xyz_corrected_pulse_area_bot_phe.*electron_lifetime_correction.*s2xy_correction_bot.*s1xyz_correction_bot;




%dp

%% Write Output File

status = LUXBinaryRQWriter_framework(settings, dp, filename_rq, data_path_rq, dp.event_number, livetime); % Should add livetime input at the end, remove event_number as settings field


