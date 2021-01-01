function [cor_xy_iq dp] = LUX_GetTheCorrectionsFromKrypton(dp, options)

% dtsetReconstruct = PositionReconstruction_MercuryGetXY(dp, rec_set, fit_set)
%
% Inputs:
% dp - The data structure with the pure krypton signal in the second dimension. 
% options  
%
%Outputs:
% factor_of_correction(n_dtime_cuts, n_phi_cuts, n_radius_cuts) - the factors of correction 
%    
%EXAMPLE OF USAGE:    
%
%factor_of_correction = LUX_GetTheCorrectionsFromKrypton(dp)    2
%
%Versioning:
% v0.00 - % 20130725 Claudio created this
% v1.10 - % 20130729 Translation in XY in such a way that we have
% the same number of events in each quadrant
% v1.20 - % 20130830 The code was modified to be more user friendly
% and to include also TAXY reconstruction. It selects automatically
% krypton events and sends the information to the LUG and saves
% automatically the table.
% v1.30 - % 20140108 Overlapped layers included to descrfibe better the corrections with X and Y                    

if nargin <2 
    options.plot = 0;
end

if ~isstruct(options)
    disp(sprintf(['LUX_GetTheCorrectionsFromKrypton: options is not ' ...
                  'a structure']));
    return
end

if ~isfield(options, 'plot')
        options.plot = 0;
end
    
if ~isfield(options, 'algorithm')
    options.algorithm = 'Mercury';
end
    
if ~isfield(options, 'computed_by')
    options.computed_by = 'Unknown';
end

if ~isfield(options, 'ComputeKryptonPeak')
    options.ComputeKryptonPeak = 0;
end

if ~isfield(options, 'Contraction')
    options.Contraction = 0.5;
end
    
if ~isfield(options, 'max_drift_time')
    options.max_drift_time = 320;
end

if ~isfield(options, 'translationXY')
    options.translationXY = 0;
end


%%
%% Sehr Unvollständig

Step_Translation = 0.01;

if ischar(dataset)
    disp(sprintf('LUXCheckReconstruction: Loading the data from %s', ...
                 dataset));
    vars_to_load = {'event_number', 'pulse_area_phe', 'peak_area_phe', 'prompt_fraction', ...
                'x_cm', 'y_cm', 'chi2', 'pulse_classification',  ...
                'hft_t10l_samples','hft_t10r_samples', 'hft_t0_samples', 'hft_t1_samples', 'hft_t50l_samples', ...
                'pulse_start_samples', 'x_cm_tmplt', 'y_cm_tmplt', 'xy_sigma_cm'}
    dp = LUXLoadMultipleRQMs_framework(dataset, vars_to_load);
elseif isstruct(dataset)
    dp = dataset;
else
    disp(sprintf('LUX_GetTheCorrectionsFromKrypton: The variable that you gave is not a path or a data structure.'));
end

if options.ComputeKryptonPeak

    %% Select single scatterer events
    
    Cut = sum(dp.pulse_classification==1,1) == 1 & sum(dp.pulse_classification==2,1) == 1;
    dp = LUX_filter(dp, Cut);
    dp = LUX_redpulses(dp, dp.pulse_classification==1, dp.pulse_classification==2);
    topchs = [1:60,121];
    botchs = [61:120,122];
    topS1 = squeeze(sum(dp.peak_area_phe(1, topchs ,:),2));
    botS1 = squeeze(sum(dp.peak_area_phe(1, botchs ,:),2));
    erS1 = dp.pulse_area_phe(1,:);

    Theta = 0.27;
    WeightedArea = topS1*cos(Theta)+botS1*sin(Theta);
    ZPosition = (botS1*cos(Theta)-topS1*sin(Theta))./WeightedArea;
    
    cutKr = inrange(WeightedArea, 80, 170);

    for j = 1:2,
        dp.dtime_us(j,:) = (dp.hft_t10l_samples(j,:) - dp.hft_t10r_samples(1,:))./100;
    end    
    
    dp = LUX_filter(dp, cutKr');
    
end

%%
if ~isfield(dp, 'dtime_us')
    dp.dtime_us = dp.z_drift_samples./100;
end

STEP_DT = 10;
STEP_DTA = 20;

dp.x_result = dp.x_cm(2,:).*0;
dp.y_result = dp.y_cm(2,:).*0;

for DT = 1:(ceil(options.max_drift_time/STEP_DT))
    Cut_DT = inrange(dp.dtime_us(2,:), max(DT*STEP_DT-0.5*STEP_DT-0.5*STEP_DTA, 0), max(min(DT*STEP_DT-0.5*STEP_DT+0.5*STEP_DTA, 365), 30));
    max(DT*STEP_DT-0.5*STEP_DT-0.5*STEP_DTA, 0)
    max(min(DT*STEP_DT-0.5*STEP_DT+0.5*STEP_DTA, 365), 30)
    if strcmp(options.algorithm, 'TAXY');
        XRec = dp.x_cm_tmplt(2,Cut_DT);
        YRec = dp.y_cm_tmplt(2,Cut_DT);
    else
        options.algorithm = 'Mercury';
        XRec = dp.x_cm(2,Cut_DT);
        YRec = dp.y_cm(2,Cut_DT);
    end
    
    clear x_cm_cor y_cm_cor;
    
    XRec_Init = XRec;
    YRec_Init = YRec;
    
    if options.translationXY
    
        TransX = Step_Translation*(1-2*(sum(XRec>0)>sum(XRec<0)));
        Step = TransX;
        
        while (sum(XRec<0)>sum(XRec>0) & sign(Step)==1) |  (sum(XRec<0)<sum(XRec>0) & sign(Step)==-1)
            XRec = XRec_Init+TransX;
            TransX = TransX+Step;
        end 
        
        TransY = Step_Translation*(1-2*(sum(YRec>0)>sum(YRec<0)));
        Step = TransY;
        while (sum(YRec<0)>sum(YRec>0) & sign(Step)==1) | (sum(YRec<0)<sum(YRec>0) & sign(Step)==-1)
            YRec = YRec_Init+TransY;
            TransY = TransY+Step;
        end 
        
        Translation(1,DT) = TransX;
        Translation(2,DT) = TransY; 
    else
        Translation(1,DT) = 0;
        Translation(2,DT) = 0; 
    end       
    
    phi = atan2(YRec, XRec).*180/pi;
    %% Find the distance to the wall
    Phi_Rec_pr = atan(YRec./XRec);
    PhiCut_pr = inrange(abs(Phi_Rec_pr), pi*1./6, pi/2);
    PhiSextant = abs(abs(Phi_Rec_pr)-PhiCut_pr*pi/3);
    WALL = (24.5*cos(pi/12)./cos(abs(abs(abs(Phi_Rec_pr)-PhiCut_pr*pi/3)-pi/12)));
    RATIO = 24.5./WALL;
    
    SHRINK = 600./(800-options.Contraction.*STEP_DT*DT);
    
    x_cm_Ak = XRec.*RATIO.*sqrt(SHRINK);
    y_cm_Ak = YRec.*RATIO.*sqrt(SHRINK);
    
    Radius = x_cm_Ak.^2+y_cm_Ak.^2;
    Radius_ceil = ceil((x_cm_Ak.^2+y_cm_Ak.^2));
    
    if options.plot
        figure(99)
        clf
        hold on
        PlotTopArray
        plot(XRec, YRec, 'b.','MarkerS',2)
        text(-24, 24, sprintf('%d', DT))
        plot(x_cm_Ak(1), y_cm_Ak(1), 'r.','MarkerS',5);
    end
    
    for pe = -179:1:180
        pea = pe+180;
        %% The radius with a cut in phi
        CutPhi = inrange(phi, pe-3, pe+3);

        if pe > -178 & pe < 178
            CutPhi = inrange(phi, pe-3, pe+3);
            adding = 0;
            while sum(CutPhi)<750;
                adding = adding + 1;
                CutPhi = inrange(phi, pe-3-adding*0.25, pe+3+adding*0.25);
            end
        elseif pe < -177
            CutPhi = inrange(phi, -180, pe+3) | inrange(phi, pe-3+360, 180);
            
            adding = 0;
            while sum(CutPhi)<750;
                adding = adding + 1;
                CutPhi = inrange(phi, -180, pe+3+0.25*adding) | inrange(phi, pe-3-0.25*adding+360, 180);
            end
            
        elseif pe > 177
            CutPhi = inrange(phi, pe-3, 180)  | inrange(phi, -180, pe+3-360);
            
            adding = 0;
            while sum(CutPhi)<750;
                adding = adding + 1;
                CutPhi = inrange(phi, pe-3-0.25*adding, 180)  | inrange(phi, -180, pe+3+0.25*adding-360);
            end
        end
        Radius_CutInPhi = Radius(CutPhi);
        
        %% The number of values that will be used
        NVALS = floor((numel(Radius_CutInPhi)-150)./600).*600;
        %% The radius orded 
        Radius_ascend = sort(Radius_CutInPhi(end-NVALS+1:end), 'ascend');
        
        clear position_in_Radius_ascend
        r2_instep = 1;
        for i=1:NVALS
            if ceil(Radius_ascend(i)) > r2_instep & r2_instep<601;
                position_in_Radius_ascend(r2_instep) = i;
                r2_instep = r2_instep+1;
            end
            i = i+1;
        end
        
        position_in_Radius_ascend(end:600) = NVALS;
        step_points = round(NVALS./600);
        
        factor_of_correction(DT, pea, :) = sqrt((step_points:step_points:NVALS)./position_in_Radius_ascend)./sqrt(SHRINK);
        factor_of_correction(DT, pea, 1) = factor_of_correction(DT, pea, 2);
        Ratio_EVT = Radius_ceil.*0+1;
        Ratio_EVT(Radius_ceil > 0 & Radius_ceil <601) = factor_of_correction(DT, pea, Radius_ceil(Radius_ceil > 0 & Radius_ceil <601));
        
        x_cm_cor(CutPhi) = XRec(CutPhi)./Ratio_EVT(CutPhi); 
        y_cm_cor(CutPhi) = YRec(CutPhi)./Ratio_EVT(CutPhi);
        
        if options.plot
            deletelastchild
            plot(x_cm_cor(CutPhi), y_cm_cor(CutPhi), 'r.','MarkerS',5)
            text(24, 24, sprintf('%d', NVALS))
            %pause
            deletelastchild
        end
    end
    dp.x_result(Cut_DT) = x_cm_cor;
    dp.y_result(Cut_DT) = y_cm_cor;
end
   
cor_xy_iq.factor_of_correction = factor_of_correction;
cor_xy_iq.Translation = Translation;

cor_xy_iq.filename_prefix= 'lux10_20130510T1300';
cor_xy_iq.global.computed_by= options.computed_by;
data_agora = clock;
cor_xy_iq.global.computed_date= sprintf('%d%02d%02d', data_agora(1), data_agora(2), data_agora(3));
cor_xy_iq.global.algorithm_name= 'xycorrections';
cor_xy_iq.global.algorithm_version= 1.2;
cor_xy_iq.global.algorithm_reconstruction = options.algorithm;
cor_xy_iq.global.drift_time= '4 to 320 microseconds';
cor_xy_iq.global.source= 'Kr83';
cor_xy_iq.global.source_id= 'Kr83';
cor_xy_iq.global.iq_type= 'xy_cor';
cor_xy_iq.global.method_used= 'ExpandingAroundZero';
cor_xy_iq.global.max_drift_time = options.max_drift_time*100;
cor_xy_iq.global.min_drift_time = 4*100;

cor_xy_iq.global.Contraction = options.Contraction;

Filename = ['cor_xy_iq_' cor_xy_iq.global.computed_date '_' cor_xy_iq.global.algorithm_reconstruction];
save(Filename, 'cor_xy_iq');