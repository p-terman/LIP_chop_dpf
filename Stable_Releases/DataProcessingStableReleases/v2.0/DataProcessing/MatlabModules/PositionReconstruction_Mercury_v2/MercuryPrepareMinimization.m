function [dp rec_set] = MercuryPrepareMinimization(dp, rec_set, lrf_iq)

% The top array will be equalized and an estimative of the position
% will be obtained
%
% [dp rec_set] = MercuryPrepareMinimization(dp, rec_set, lrf_iq)
%
% Inputs:
%         dp     - the data array
% Outputs:
%         dp     - the data array with the field dp.reconstructed
% Example usage:
% 
%
% 20130213 CFPS - Created

%%
%% Input Check
%%


%%
%% Some initial definitions
%%

NEVTS = size(dp.peak_area_phe, 3);
NPULS = size(dp.peak_area_phe, 1);
NPMTS = min(size(dp.peak_area_phe, 2), 122);

[pmt_pos sextant_arrangement] = LUXPMTArray; 
for i=1:122,
    PMT_r(i,:)=pmt_pos(find(sextant_arrangement==i),:);
    the_sextant(i,:) = Find_Sextant(i, pmt_pos, sextant_arrangement);
end

%%
%% Array Equalization
%%
dp.peak_area_eq = zeros(NPULS,NPMTS, NEVTS);

% Correct with the gains
if rec_set.verbose > 1
    fprintf(['MercuryPrepareMinimization: The array will be equalized\n']);
end

gains = lrf_iq.QE.QE_Values(1:NPMTS);    
dp.peak_area_eq = max(dp.peak_area_phe(:,1:NPMTS,:)./(repmat(repmat(gains, [NPULS 1]), [1 1 NEVTS])), 0);

% The PMTs that we are not going to use we are going set them to zero.
dp.peak_area_eq(:,rec_set.PMTS_To_Use==0,:) = 0;

% It is just a precausion if some of these values are nan or inf
dp.peak_area_eq(isinf(dp.peak_area_eq) | isnan(dp.peak_area_eq)) = 0;

Cut = isinf(dp.peak_area_eq) | isnan(dp.peak_area_eq);

if Cut(:)>0 & rec_set.verbose > 1
    fprintf(['MercuryPrepareMinimization: In the array equalization ' ...
             'we have some infinite/nan values for the peak areas\' ...
             'n']);
end
    

% Check the pulse height. If the peak height is smaller than the
% threshold of the struck the peak area should be zero.

if isfield(dp, 'peak_height_mV')
    if ~isfield(lrf_iq, 'peak_height_mV_threshold')
        lrf_iq.peak_height_mV_threshold = 1.5;
    end
    dp.peak_area_eq(dp.peak_height_mV < 1.5) = 0;
elseif isfield(dp, 'peak_height_phe_per_sample')
    dp.peak_area_eq(dp.peak_height_phe_per_sample < 0.15) = 0;  
end        

% The PMTs that we are not going to use we are going set them to zero.
dp.peak_area_eq(:,rec_set.PMTS_To_Use==0,:) = 0;

% It is just a precausion if some of these values are nan or inf
dp.peak_area_eq(isinf(dp.peak_area_eq) | isnan(dp.peak_area_eq)) = 0;

Cut = isinf(dp.peak_area_eq) | isnan(dp.peak_area_eq);

if Cut(:)>0 & rec_set.verbose > 1
    fprintf(['MercuryPrepareMinimization: In the array equalization ' ...
             'we have some infinite/nan values for the peak areas\n']);
end

% Check the pulse height. If the peak height is smaller than the
% threshold of the struck the peak area should be zero.


%% Definition of the total area

top_area = squeeze(sum(dp.peak_area_eq(:,[1:60 121],:), 2));
bot_area = squeeze(sum(dp.peak_area_eq(:,[61:120 122],:), 2));
dp.pulse_area_eq = top_area+bot_area;

dp.peak_area_eq = double(dp.peak_area_eq);
dp.pulse_area_eq = double(dp.pulse_area_eq);

if ~isfield(rec_set, 'getapproximatedposition')
    rec_set.getapproximatedposition = 'GuessTri';
    if rec_set.verbose > 1
        fprintf(['MercuryPrepareMinimization: We will get an ' ...
                 'approximated position using the GessTri method']);
    end
end


if ~isfield(rec_set, 'overwritepositions')
    overwritepositions = 1;
else
    overwritepositions = rec_set.overwritepositions;
end

%%
%% Initialization of the variables
%%

if ~isfield(dp, 'x_cm') | ~isfield(dp, 'y_cm') 
    overwritepositions = 1;
    dp.x_cm(NPULS,NEVTS)  = 0; dp.y_cm(NPULS,NEVTS)  = 0;
end    

if ~isfield(dp, 'chi2') | ~isfield(dp, 'rec_dof') 
    dp.rec_dof(NPULS,NEVTS)  = 0; dp.chi2(NPULS,NEVTS)  = 0;
end

if ~isfield(dp, 's2_rec')
    dp.s2_rec(NPULS,NEVTS)  = 0;
end

if rec_set.compute_sd
    dp.sd_radius_inf = zeros([NPULS, NEVTS], 'single');
    dp.sd_radius_sup = zeros([NPULS, NEVTS], 'single');
    dp.sd_phiXR = zeros([NPULS, NEVTS], 'single');
end

if ~isfield(dp, 'peak_area_rec')
    dp.peak_area_rec = ones(NPULS, 122, NEVTS).*(-100);
end

%%
%%
%%%
if overwritepositions
    
    NumberOfTopPMTsWithSignal = squeeze(sum(dp.peak_area_phe(:, [1:60 121],:)>1, 2));
    for e = 1:NEVTS,
        for p = 1:NPULS,
            peak_area_event = squeeze(dp.peak_area_eq(p, :,e));
            peak_area_event(32) = mean(peak_area_event([32 33 35 36]));
            peak_area_event(5) = mean(peak_area_event([1 2 6 8 57 54]));
            peak_area_event(41) = mean(peak_area_event([34 45 42]));
            
            if dp.reconstructed(p, e) == 1 & NumberOfTopPMTsWithSignal(p, e) > 7                
                %% This is to solve the problem of the PMTs
                %% that are kaputt.
                try 1;
                    [dp.x_cm(p, e) dp.y_cm(p, e) Status] = LUXLRF_Guess_Tri(peak_area_event, PMT_r, the_sextant, 3, 1);
                catch err
                    if rec_set.verbose > 1
                        fprintf(['MercuryPrepareMinimization: The ' ...
                                 'standard method gave an error and I will try a different thing (provisory - right now I am putting these events in the centre of the chamber)\n']);
                    end
                    dp.x_cm(p,e)  = 2*rand-1;
                    dp.y_cm(p,e)  = 2*rand-1;
                end
                %%                if dp.x_cm(p,e)^2+dp.y_cm(p,e)^2 > 600
                %%    peak_area_event(peak_area_event<1) = 0;
                %%    dp.x_cm(p,e)  = mean((PMT_r([1:60 121],1)'.*peak_area_event([1:60 121]))./(sum(peak_area_event([1:60 121]))));
                %%    dp.y_cm(p,e)  = mean((PMT_r([1:60 121],2)'.*peak_area_event([1:60 121]))./(sum(peak_area_event([1:60 121]))));
                %%                        
                %% end
            elseif dp.reconstructed(p, e) == 1 & NumberOfTopPMTsWithSignal(p, e) < 7

                peak_area_event(peak_area_event<1) = 0;
                NPMTSE = sum(peak_area_event([1:60 121])>0);
                if NPMTSE>0
                    dp.x_cm(p,e)  = sum(PMT_r([1:60 121],1)'.*peak_area_event([1:60 121]))./(sum(peak_area_event([1:60 121]))*NPMTSE);
                    dp.y_cm(p,e)  = sum(PMT_r([1:60 121],2)'.*peak_area_event([1:60 121]))./(sum(peak_area_event([1:60 121]))*NPMTSE);
                else
                    dp.x_cm(p,e)  = 2*rand-1;
                    dp.y_cm(p,e)  = 2*rand-1;
                end
            end
        end
    end
    % remove infinities
    dp.x_cm(isinf(dp.x_cm)) = 2*rand-1;
    dp.y_cm(isinf(dp.y_cm)) = 2*rand-1;
end
    
    dp.x_cm(~dp.reconstructed) = -100;
    dp.y_cm(~dp.reconstructed) = -100;
    dp.chi2(~dp.reconstructed) = -100;
    dp.s2_rec(~dp.reconstructed) = -100;
    dp.rec_dof(~dp.reconstructed) = -100;


