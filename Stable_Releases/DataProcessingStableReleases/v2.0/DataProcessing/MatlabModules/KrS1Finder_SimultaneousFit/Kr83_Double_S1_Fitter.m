function fit_struct = Kr83_Double_S1_Fitter(time_samples,data_phe_per_sample, minAmp, risetime_samples, falltime_samples, evt, fit_options)
%
% fit_struct = Kr83_Double_S1_Fitter(time_samples,data_phe_per_sample)
%
% This function operates on the S1 pulses for Kr-83. It does the following:
%
% 1) It first performs a very sensitive pulse-finding algorithm to idenfity
%    the locations of S1a and S1b. 
% 2) It then performs a simultaneous double S1 fit.
% 
% It gives the area, height, decay/rise times, as well as chisq and dof.
% 
% Inputs:
%   time_samples - the sumpod time in samples for the Kr-83 S1 pulses
%   data_phe_per_sample - the sumpod data in phe/sample for the Kr-83 S1 pulses
%
% Outputs:
%   fit_struct
%      .Kr83_s1a_pulse_area_phe
%      .Kr83_s1b_pulse_area_phe
%      .Kr83_s1a_amplitude_phe_per_sample
%      .Kr83_s1b_amplitude_phe_per_sample
%      .Kr83_s1a_time_offset_samples
%      .Kr83_s1b_time_offset_samples
%      .Kr83_s1ab_chisq                  (global fit chisq)
%      .Kr83_s1ab_dof                    (global fit degrees of freedom)
% 
% Versioning:
%   v1.0 20140827 CHF&AC created
%   v1.1 20140827 CHF&AC fixed issue where S1b could happen before S1a,
%                        tightened fit bounds
%   v2   20141020 SAH replaced fitted-shape version with fixed-shape version (risetimes and falltimes set in xml)



% just flags for the inspection plots
plot_flag = 0;
save_flag = 0;




%% Flatten the data and filter it
% Flatten data (ie, when abs(data)< 1phe/samp, data->0)
data_phe_per_sample_flat = data_phe_per_sample;
data_phe_per_sample_flat(abs(data_phe_per_sample_flat) < 1) = 0;

% Low-pass filter data trace with 2nd order Butterworth
[b,a] = butter(2,0.40,'low');
data_phe_per_sample_filt = filtfilt(b,a,data_phe_per_sample_flat);

% Get diff of filtered trace and flatten it again (when abs(diff(data))<0.5, diff(data)->0   )
diff_data_phe_per_sample_filt = [0 diff(data_phe_per_sample_filt')];
diff_data_phe_per_sample_filt(abs(diff_data_phe_per_sample_filt) < 0.5) = 0;

% Find transitions from positive diff to negative diff (ie, pulse tops)
positive_inds = diff_data_phe_per_sample_filt > 0;
negative_inds = diff_data_phe_per_sample_filt <= 0;
pulse_mask = (positive_inds*1 + negative_inds*-1 + 1)/2;

%% Find the initial locations of S1a and S1b

% Get pulse-finding box edges, left and right
left_edges  = find(diff(pulse_mask)== +1);
right_edges = find(diff(pulse_mask)== -1);
nothingfound = 0;
if isempty(left_edges)
    [val left_edges] = max(data_phe_per_sample);
    nothingfound = 1;
end

if length(left_edges) > 1   % If we found multiple pulse tops

    % Initialize vector of max values for each region
    max_vals = zeros(1,length(left_edges));
    max_inds = zeros(1,length(left_edges));

    % Loop for each left edge, find max value within pulse-finding box
    for ii = 1:length(left_edges)
        if length(right_edges) < ii
            right_edges(ii) = left_edges(ii)+1;
        end
        range = left_edges(ii):right_edges(ii);
        [max_vals(ii), ind] = max(data_phe_per_sample(range));
        max_inds(ii) = range(ind);    
    end
    [max_vals_ordered,I] = sort(max_vals,'descend');   % Sort pulses by max val
    
    if max_vals_ordered(2)>6          % if the second-largest pulse is not just a little blip
      inds_kr83 = max_inds(I(1:2));   % select the largest two pulses
    else                                              % if the second-largest pulse is tiny
      inds_kr83 = [max_inds(I(1))  max_inds(I(1))+3]; % select only the largest pulse (twice)
    end
    
    
elseif length(left_edges) ==1 && nothingfound == 0  % if one and only one pulse was found
    
    range = left_edges:right_edges;
    [max_val, ind] = max(data_phe_per_sample(range));
    max_ind = range(ind);
    inds_kr83 = [max_ind max_ind+3];
    
else % if zero pulses were found
    
    inds_kr83 = [left_edges left_edges+3];
end


%% Now do fits
%   (risetime and falltime values FIXED in the xml settings)

% Initial parameters

% S1a
init_params(1) = data_phe_per_sample(inds_kr83(1)); % S1a max
init_params(2) = time_samples(inds_kr83(1)-1); % S1a time of max

% S1b
init_params(3) = data_phe_per_sample(inds_kr83(2)); % S1b max
init_params(4) = time_samples(inds_kr83(2)-1); % S1b time of max

% we have to combine the time and the pulse shape quantities into a single input quantity here
input.tt = time_samples;
input.risetime_samples = risetime_samples;
input.falltime_samples = falltime_samples;

% min and max fit ranges
param_mins  = [minAmp*1.5           ...
               min(time_samples)+1   ...
               minAmp               ...
               min(time_samples)+1 ];
param_maxes = [Inf                   ...
               max(time_samples)-5   ...
               Inf                   ...
               max(time_samples)-5 ];


% Do the fit
[temp_param, resnorm, residual, exitflag] = lsqcurvefit(@LUXDoubleExpFcn_fixedshape_framework,...
                                                           init_params,...
                                                           input,...
                                                           data_phe_per_sample,...
                                                           param_mins, param_maxes,...
                                                           fit_options);

if temp_param(2) > temp_param(4)   % if the S1b time is before the S1a time, swap the incorrect S1a S1b assignments
    param(1:2) = temp_param(3:4);
    param(3:4) = temp_param(1:2);
else
    param = temp_param;
end
    
% turn the fit parameters into pulses
S1a_param = [param(1), falltime_samples, param(2), risetime_samples];
S1b_param = [param(3), falltime_samples, param(4), risetime_samples];
s1a_fit = LUXExpFcn_framework(S1a_param , time_samples);
s1b_fit = LUXExpFcn_framework(S1b_param , time_samples);
s1t_fit = s1a_fit + s1b_fit;

% Get maxes and areas
[max_s1a,maxind_s1a] = max(s1a_fit);
[max_s1b,maxind_s1b] = max(s1b_fit);
area_s1a = trapz(s1a_fit);
area_s1b = trapz(s1b_fit);

% Calculate chisq and dof
ct = inrange( time_samples, time_samples(maxind_s1a)-10,time_samples(maxind_s1b)+30);
s1t_cut = s1t_fit(ct);
percent_of_max_for_baseline_error = 0.1;
cutlowval = abs(s1t_cut) <= percent_of_max_for_baseline_error*max(abs(s1t_cut));
sd = sqrt(abs(s1t_cut));
sd((sd == 0) | cutlowval) = std(s1t_cut(cutlowval));
chisq = sum( (data_phe_per_sample(ct) - s1t_fit(ct)).^2  ./ sd.^2 );
dof = length(ct) - length(param); 

% Reshape fit parameters into data structure, including chisq and dof
fit_struct.Kr83fit_s1a_area_phe   = area_s1a;
fit_struct.Kr83fit_s1b_area_phe   = area_s1b;
fit_struct.s1a_t_samples          = time_samples(maxind_s1a);
fit_struct.s1b_t_samples          = time_samples(maxind_s1b);
fit_struct.Kr83fit_dt_samples     = time_samples(maxind_s1b) - time_samples(maxind_s1a);
fit_struct.Kr83fit_chisq      = chisq;
fit_struct.Kr83fit_dof        = dof;
fit_struct.Kr83fit_exitflag   = exitflag;

% ----------------- All fitting ---------------- %





%% Plot or save figure if needed

if plot_flag

    figure(1);
    clf
    
    % pulse finding
    subplot(2,1,1)
    plot(time_samples, 10*pulse_mask, 'm.', 'markersize', 50, 'color', [1 0.8 0.5]);
    hold on;
    plot([time_samples(inds_kr83(1)) time_samples(inds_kr83(1))],[-0.1*max(data_phe_per_sample),1.4*max(data_phe_per_sample)],'r-','linewidth', 4);
    plot([time_samples(inds_kr83(2)) time_samples(inds_kr83(2))],[-0.05*max(data_phe_per_sample),1.35*max(data_phe_per_sample)],'b-','linewidth', 4);
    plot(time_samples,data_phe_per_sample,'ko','markers',6,'markerface',[0.7 0.7 0.7], 'color', [0.7 0.7 0.7]);
    plot(time_samples,data_phe_per_sample_filt,'k-')
    plot(time_samples,diff_data_phe_per_sample_filt,'g-');
    hold off;
    xlim([min(time_samples),min(time_samples)+100])
    ylim([-0.1*max(data_phe_per_sample),1.4*max(data_phe_per_sample)])
    xlabel('samples')
    ylabel('phe/sample')
    title(['event ' num2str(evt) ]);
    legend('pulse mask (x10)', 'initial guess for S1a peak','initial guess for S1b peak', 'all data', 'all data, low-pass filtered', 'derivative, with some funny business');
    
    % pulse fitting
    subplot(2,1,2)
    plot(time_samples,data_phe_per_sample,'ko','markers',6,'markerface',[0.7 0.7 0.7], 'color', [0.7 0.7 0.7]);
    hold on
    plot(time_samples,s1t_fit,'k-', 'color', [0.2 0 0.2]);
    plot(time_samples,s1a_fit,'r-', 'color', [1 0.2 0.2]);
    plot(time_samples,s1b_fit,'b-', 'color', [0.4 0.4 1]);
    grid on
    xlim([min(time_samples),min(time_samples)+100])
    ylim([-0.1*max(data_phe_per_sample),1.4*max(data_phe_per_sample)])
    xlabel('samples')
    ylabel('phe/sample')
    title(['S1b amplitude = ' num2str(S1b_param(1)) '[phe/samp]    chisq/dof=' num2str(fit_struct.Kr83fit_chisq/fit_struct.Kr83fit_dof) ]);

           
end

if save_flag
    
set(gcf,'PaperUnits'        , 'centimeters',...
        'PaperSize'         , [30, 30],...
        'PaperPosition'     , [0, 0, 30, 30],...
        'PaperPositionMode' , 'manual');
print('-dpdf','-r100',['/home/scott/analysis/LUX/Kr_and_H3_withfieldvariation/run_module/fitplot_' num2str(10000+evt) '.pdf']);

end