function dtsetReconstruct = PositionReconstruction_MercuryGetXY(dp, rec_set, lrf_iq)

% dtsetReconstruct = PositionReconstruction_MercuryGetXY(dp, rec_set, lrf_iq)
%
% Inputs:
% dp - The data structure of our dataset. 
% rec_set - Reconstruction settings.
% lrf_iq - LRF interisting quantities.
%
%Outputs:fit

% dtsetReconstruct - The dataset with the positions   
%    
%EXAMPLE OF USAGE:    
%
%dtsetReconstruct = MercurySim(dtset, fit)    
%    
%
%Versioning:
% v1.00 - CFPS All this file was reestructured
% 20130403 pfs - Multiple patches to obtain smooth operation in the new DP:
%                lines 171-175 (added)
%                lines 364, 386, 397, 480 (defined a previously undefined variable)
% 20130411 CFPS some modifications 
% 20140625 CFPS Version 2.0.1 small patches 	
% 20140712 CFPS Version 2.1 The LRFs can now  be function of the radius 
% 20140714 CFPS Version 2.1.1 Small Patch to use different saturation points for each PMT
% 20140716 CFPS Version 2.1.2 Small Patch. Check if the variable spike_count exists
% 20140723 CFPS Version 2.1.3 Small Patch. Change of the variable type
% of the variables reconstructed and s2_rec.
% 20140724 CFPS Version 2.1.4 Bug related with the photon counting corrected.
% 20140725 CFPS Version 2.1.5 Bug related with the computation of
% uncertainties corrected
% 20140730 CFPS Version 2.1.6 Bug related with the computation of the peak areas corrected. It only
%                             afected the PMT 121.
% 20140808 CFPS Version 2.1.7 Bug related with the computation of the uncertainties in the photon counting mode.

if rec_set.verbose > 1
    sprintf('Creating ntuple tables\n') 
end
tic


%%
%% Part 0.1 Definition and Initialization of some variables
%%
Min_Handle = @fminsearch; %Min_Handle = @fminunc; %Min_Handle = @fminsearch; Min_Handle = @fmincon;
minim_opt = optimset('Display','off', 'LargeScale','off');

[NPULS NPMTS NEVTS] = size(dp.peak_area_phe); %% How the number of events is chosen in each module

if ~isfield(lrf_iq.ref, 'variabledependence')
    fit_set.ref.variabledependence = zeros(1, 61);
end
    
if ~isfield(dp, 'peak_area_rec') & rec_set.get_estimated_peak_areas == 1
    dp.peak_area_rec(NPULS,1:122, NEVTS) = 0;
end
gains = lrf_iq.QE.QE_Values(1:NPMTS);
gains(rec_set.PMTS_Working==0) = 1;

dp.x_cm_old = dp.x_cm; % The old values of XRec and YRec are copied to a new variable.
dp.y_cm_old = dp.y_cm;
dp.minimization_flag = zeros([NPULS, NEVTS], 'uint8');

%This is because the name of a variable changed
if isfield(dp, 'pmt_2pct_saturation_flag')
    pmt_sat_flag = 'pmt_2pct_saturation_flag';
elseif isfield(dp, 'pmt_saturation_flag')
    pmt_sat_flag = 'pmt_saturation_flag';
else
    sprintf(['NO SATURATION FLAG. THIS WILL NOT WORK WELL WITHOUT THAT']);
    pmt_sat_flag = 'pmt_saturation_flag_prov';
    dp.(pmt_sat_flag) = zeros(NPULS, NPMTS,NEVTS);
end

chi2_old = dp.chi2;

if ~isfield(dp, 'chi2') | ~isfield(dp, 'rec_dof')
    dp.chi2 = zeros([NPULS, NEVTS], 'single');
    dp.rec_dof = zeros([NPULS, NEVTS], 'single');
end    

dp.chi2(~dp.reconstructed) = -1000;

dp.hft_t10r_samples=int32(dp.hft_t10l_samples);

if isfield(dp, 'pulse_width')
    pulse_width = dp.pulse_width;
elseif isfield(dp, 'hft_t10r_samples')
    pulse_width = dp.hft_t10r_samples-dp.hft_t10l_samples;
elseif isfield(dp, 'pulse_length_samples')
    pulse_width = dp.pulse_length_samples;
else
    pulse_width = dp.pulse_end_samples-pulse_start_samples;
end

if ~isfield(dp, 'spike_count') & rec_set.PhotonCount_ML > 5
    disp(sprintf(['\n **** \n PositionReconstruction_Mercury: NO SPIKE COUNT RQ. Please change the settings variable rec_set.PhotonCount_ML to -1\n or run the spike count module previously\n The reconstruction will proceed using the old pulse areas method \n *****\n']));
    rec_set.MLM_maxphe = rec_set.PhotonCount_ML; rec_set.PhotonCount_ML = -1;
end

%%
%% Part 0.2 Definition of the PMTs geometry
%%

%Creation of the arrays of the PMTs and of the sextant. These
%arrays are numbered according to the PMTs. These Variables are
%defined as global below this point.

global topchs
topchs = [1:60,121];

% This creates a PMT Positions vector map and the sextant maps
[pmt_pos sextant_arrangement] = LUXPMTArray; 
global PMT_r;
global the_sextant;
for i=1:122,
    PMT_r(i,1:3)=pmt_pos(find(sextant_arrangement==i),:);
    the_sextant(i,:) = Find_Sextant(i, pmt_pos, sextant_arrangement);
end

global PMT_Radius;
global PMTG;
[PMT_Radius PMTG] = LUXLRF_GroupDefinition;

%%
%% Part 0.3 LRF initiation
%%

% some definitions of the LRF functions
lrf_iq.ref.Tau_Decay_Function_In = eval('inline(''Amp*exp(-x.^2*InvSig) + K '',''x'', ''Amp'',''InvSig'', ''K'');');
lrf_iq.ref.Amplitude_Function = eval(lrf_iq.ref.Amplitude_Function_String); %% This is the third function
if strcmp(lrf_iq.rad.FunctionToFit,'bivariate_cauchy')
    lrf_iq.rad.CallFunction_In = eval('inline(''((1+x.^2*par(2)).^(-3/2)).*par(1)+x.*par(3)+par(4)'', ''x'', ''par'')');
elseif strcmp(lrf_iq.rad.FunctionToFit,'extended_cauchy')
    lrf_iq.rad.CallFunction_In = eval('inline(''((1+x.^2*par(2)).^(-3/2)).*par(1)+exp(-x.*par(4)).*par(3)+x.*par(5)+par(6)'', ''x'', ''par'')');
elseif strcmp(lrf_iq.rad.FunctionToFit,'extended_cauchy2')
    lrf_iq.rad.CallFunction_In = eval('inline(''((1+x.^2*par(2)).^(-3/2)).*par(1)+exp(-x.^2*par(4)).*par(3)+x.*par(5)+par(6)'', ''x'', ''par'')');
elseif strcmp(lrf_iq.rad.FunctionToFit,'extended_cauchy3')
    lrf_iq.rad.CallFunction_In = eval('inline(''((1+x.^2*par(2)).^(-3/2)).*par(1)+exp(-x.*par(4)).*par(3)+exp(-x.*par(6)).*par(5)+x.*par(7)+par(8)'', ''x'', ''par'')');
end

% Definition of the maximum likelihood function
if rec_set.MLM_maxphe > 0 %% Code necessary for the
    %% definitions
    Gaussian_MLM = inline('exp(-((x-N).^2)./(2*(sigmas.*N + sigma_zero)))./sqrt(2*pi*(sigmas.*N + sigma_zero))', 'x', 'sigmas', 'N', 'sigma_zero');
end     

% failsafe method
if isfield(rec_set, 'fail')
    if rec_set.fail == 1
        dp.x_cm(:) = -100;
        dp.y_cm(:) = -100;
        dp.chi2(:) = -100;
        return
    end
end

if rec_set.verbose > 1
    disp(sprintf('******************************************************************************************************************'));
    disp(sprintf('Minimization Started'));
    disp(sprintf('******************************************************************************************************************'));
end
    
if strcmp(rec_set.algorithm, 'mercuryi_splines')
    % Get the splines look-up table
    
    Delta_Phi_sp = 1:1800;
    lrfmat.spval = zeros(61, 160, numel(Delta_Phi_sp));
    for Radius_sp = 1:160
        for PMTGroup = 2:9,
            cutPMTs = PMT_Radius == PMTG(PMTGroup);
            lrfmat.spval(cutPMTs, Radius_sp, 1:numel(Delta_Phi_sp)) = repmat(ppval(lrf_iq.sp(PMTGroup,Radius_sp), Delta_Phi_sp.*pi/(3600)-pi/(3600)), [sum(cutPMTs) 1]);
        end
    end    
    lrfmat.sp = lrf_iq.sp;
    
    %% The direct component
    RhoSquaredT10 = 1:36000;
    Rho = 60-sqrt(RhoSquaredT10./10);
    lrfmat.Rho = lrf_iq.rad.CallFunction_In(Rho, lrf_iq.rad.Parameters);
    
    [pmt_pos sextant_arrangement] = LUXPMTArray; 
    for i=1:122,
        lrfmat.PMT_r(i,1:3)=pmt_pos(find(sextant_arrangement==i),:);
    end
    
    lrfmat.PMT_phi = atan2(lrfmat.PMT_r(:,2), lrfmat.PMT_r(:,1));
    
    lrfmat.topchs = [1:60 121];
    lrfmat.topcut = PMT_r(:,3) > 0;
    lrfmat.sp = lrf_iq.sp;
    
    chi2_Handle = @chi2_BSplines;
    [lrfmat.PMT_Radius lrfmat.PMTG] = LUXLRF_GroupDefinition;
    
elseif strcmp(rec_set.algorithm, 'mercuryi_functional')
    
    if ~isfield(lrf_iq.ref, 'Amp2') | ~isfield(lrf_iq.ref, 'Gam2')
        lrf_iq.ref.Amp2 = zeros(122, 1);
        lrf_iq.ref.Gam2 = zeros(122, 1);
    end
    
    success = 0;
    if exist([char(rec_set.file_lrfmat) '_Started'], 'file') == 2
        max_t = 1;
        while exist([char(rec_set.file_lrfmat) '_Finished'], 'file') ~= 2 & max_t < 40;
            pause(1);
            disp(sprintf('waiting for the file with the look up table...\n', rec_set.file_lrfmat));
            max_t = max_t +1;
        end
        
        if max_t < 40
            if rec_set.verbose > 1;
                disp(sprintf('loading the file %s with the look-up-table information.\n', rec_set.file_lrfmat));
            end
            load(char(rec_set.file_lrfmat));
            if exist('lrfmat')
                if ~isequal(lrf_iq, lrfmat.lrf_iq) | ~isequal(rec_set, lrfmat.rec_set)
                    disp(sprintf('The stored LRF is different from the current one. The loaded file will not be used\n'));
                    success = 0;
                    delete(rec_set.file_lrfmat);
                    delete([char(rec_set.file_lrfmat) '_Started']);
                        delete([char(rec_set.file_lrfmat) '_Finished']);
                else
                    success = 1;
                end
            end
        else
            success = 0;
        end
    else
        success = 0;
    end
    
    if ~success
        
        started.a = 'initial_flag';
        save([char(rec_set.file_lrfmat) '_Started'], 'started');
        
        lrfmat.lrf_iq = lrf_iq;
        lrfmat.variabledependence = lrf_iq.ref.variabledependence;
        lrfmat.rec_set = rec_set;

        if rec_set.verbose > 1;
            disp(sprintf('A new look-up table is being calculated.\n'));
        end

        RadiusSqT10 = 1:65000;
        lrfmat.Tau = zeros(61, size(RadiusSqT10, 2));
        lrfmat.Amp = zeros(61, size(RadiusSqT10, 2));
        dis_event_wall = 24.5-sqrt(RadiusSqT10./100);
        
        lrfmat.Amp2 = lrf_iq.ref.Amp2;
        lrfmat.Gam2 = lrf_iq.ref.Gam2;
        
        % PMTs defined by groups or individually
        if size(lrf_iq.ref.Exp_Amplitude, 1)<10
            for PMTGroup = 1:9,
                cutPMTs = PMT_Radius(topchs) == PMTG(PMTGroup);
                TauD = repmat(lrf_iq.ref.Tau_Decay_Function_In(dis_event_wall, lrf_iq.ref.Tau_Decay(PMTGroup,1), lrf_iq.ref.Tau_Decay(PMTGroup,2), lrf_iq.ref.Tau_Decay(PMTGroup,3)), [122 1]);        
                AmpD = repmat(lrf_iq.ref.Amplitude_Function(dis_event_wall, lrf_iq.ref.Exp_Amplitude(PMTGroup,:)), [122 1]);
                lrfmat.Tau(cutPMTs,:) = TauD(cutPMTs,:);
                lrfmat.Amp(cutPMTs,:) = AmpD(cutPMTs,:);
            end
            PMTGroup
        else
            for PMTSingle = [1:60]
                lrfmat.Tau(PMTSingle,:) = lrf_iq.ref.Tau_Decay_Function_In(dis_event_wall, lrf_iq.ref.Tau_Decay(PMTSingle,1), lrf_iq.ref.Tau_Decay(PMTSingle,2), lrf_iq.ref.Tau_Decay(PMTSingle,3));        
                lrfmat.Amp(PMTSingle,:) = lrf_iq.ref.Amplitude_Function(dis_event_wall, lrf_iq.ref.Exp_Amplitude(PMTSingle,:));
            end
        end
        
        RhoSquaredT10 = 1:36000;
        Rho = 60-sqrt(RhoSquaredT10./10);
        for PMTSingle = [1:60 121]
            if abs(sum(lrf_iq.rad.Parameters_EachPMT(PMTSingle,:)))>0,
                lrfmat.Rho(min(PMTSingle,61),:) = lrf_iq.rad.CallFunction_In(Rho, lrf_iq.rad.Parameters_EachPMT(PMTSingle,:));
            else
                lrfmat.Rho(min(PMTSingle,61),:) = lrf_iq.rad.CallFunction_In(Rho, lrf_iq.rad.Parameters);
            end
        end
        
        [pmt_pos sextant_arrangement] = LUXPMTArray; 
        for i=1:122,
            lrfmat.PMT_r(i,1:3)=pmt_pos(find(sextant_arrangement==i),:);
        end
        
        lrfmat.topchs = [1:60 121];
        lrfmat.topcut = PMT_r(:,3) > 0;
        lrfmat.notwor(1:61) = rec_set.PMTS_Working(topchs)==0;
        lrfmat.radial_only = zeros(1, 61);

        R2 = 1:6200;
        lrfmat.chi2C = 0;
        
        phi = 0:pi*0.0001:pi*0.5;
        lrfmat.disphi = 24.5*cos(pi/12)./cos(abs(abs(abs(phi)-inrange(abs(phi), pi*1./6, pi/2)*pi/3)-pi/12));
        lrfmat.disphi(end) = 24.5;
        save(rec_set.file_lrfmat, 'lrfmat');
        started.a = 'final_flag';
        save([char(rec_set.file_lrfmat) '_Finished'], 'started');
    else
        if rec_set.verbose > 1;
            disp(sprintf(['***WARNING*** To improve speed the code ' ...
                          'is using a pre-calculated lookup table ' ...
                          'stored in the file \n %s. ***WARNING***\n'], rec_set.file_lrfmat));
        end
    end
    
    %% Miminization handles definition
    if (rec_set.only_radial_component) chi2_Handle = @chi2_OnlyRadial; else chi2_Handle = @chi2_Matricial; end
    MLM_Handle = @MLM_Matricial;
    PC_Handle = @MLM_PhotonCount;
    
    if (rec_set.MLM_maxphe > 0) sphes = lrf_iq.sphe(topchs).^2; end
    
else
    lrfmat = lrf_iq;
    chi2_Handle = @chi2_FuncPar; % This is the default function
end
    
if rec_set.verbose > 1
    disp(sprintf('Finishing Initialization in %.2f s\n', toc));
end

t_0 = cputime; % CPU time, used to control the time of the minimization
tic

%%
%% Start of the reconstruction code
%%

if rec_set.verbose ==1
    fprintf(sprintf('Reconstructing ----  \t%d%%', 00.00));
end

EVTS = 0; % Counting the number of variables
TotEVTS = sum(dp.reconstructed(:)>0);
for e = 1:NEVTS,     
    for p = 1:NPULS
        if dp.reconstructed(p, e)  >= 1 
            %% Step 0
            %% Printing initial Information
            %%
            
            if rem(EVTS,50)==0 & rec_set.verbose>0
                disp(sprintf('---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'));
                disp(sprintf(['%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\' ...
                              't%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\' ...
                              't%s\t%s\t%s\t%s\t%s\t%s'], 'Event', 'Pulse','MT', 'CPU','Sum' ,'Amp_O','Amp_M','R_Ori','R_Min','Phi_O', 'Phi_M','x_Or', 'y_Or','x_cm', 'y_cm','S_R I','S_R S', 'S_Phi_I','S_Phi_S', 'Flag', 'Num', 'N_PMTs', 'Chi2_O', 'Chi2_N'));
                disp(sprintf('%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s', '#',      '#',   '',     '(s)', 'phe', 'phe',   'phe', '(cm)',' (cm)','(cm)',' (cm)','(deg)','(deg)', '(cm)','(cm)', '(cm)','(cm)','(deg)','(deg)', '--', 'Iter','#', 'red.', 'red.'));
                disp(sprintf('---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'));
                valchi = dp.chi2(1:NPULS,1:e); valchi_old = chi2_old(1:NPULS,1:e); valrecos = dp.reconstructed(1:NPULS,1:e);
                disp(sprintf('Event Number: %d - Old Mean of Chi2: %.3f New Mean of the chi2: %.3f\n', e, mean(valchi_old(valrecos>0)), mean(valchi(valrecos>0))));
                disp(sprintf('---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------'));
            end
            
                        
            EVTS = EVTS + 1;
            clear p_min;

            if rec_set.verbose ==1
                if rem(EVTS,10)==1
                    %fprintf('\b\b\b\b\b');
                    fprintf('%2.2f%%\n', 100*EVTS/(TotEVTS));
                end
            end
            
            
            %% Step 
            %% Computation of Some Basic Variables and the initial value for the fit
            %%
            
            % Some variables that are needed

            PEAK_AREA = dp.peak_area_eq(p,topchs,e );
            PULSE_AREA = dp.pulse_area_eq(p,e);
            
            if (abs(dp.x_cm_old(p,e))<0.01) dp.x_cm_old(p,e) = dp.x_cm_old(p,e)+2*rand-1; end
            if (abs(dp.y_cm_old(p,e))<0.01) dp.y_cm_old(p,e) = dp.y_cm_old(p,e)+2*rand-1; end

            
            if dp.rec_energy_flag(p,e)
                In_Guess = [min(max(dp.x_cm_old(p,e), -24.5), 24.5) min(max(dp.y_cm_old(p,e), -24.5), 24.5) PULSE_AREA];
            elseif rec_set.search_for_double_scatterer
                PEAKA = dp.peak_area_eq(p,:,e );
                if PMT_r(find(PEAKA==max(PEAKA)), 3) > 0
                    In_Guess = [PMT_r(find(PEAKA==max(PEAKA)), 1) PMT_r(find(PEAKA==max(PEAKA)), 2)];
                    PULSE_AREA = 10*max(PEAKA);
                else
                    disp(sprintf('The maximum of the peak area is smaller than 3\n'));
                    In_Guess = [dp.x_cm_old(p,e) dp.y_cm_old(p,e)];
                end
                rec_set.maximum_distance_PMTevent = 7;
                rec_set.PMT_MinNum = 9;
            else
                In_Guess = [dp.x_cm_old(p,e) dp.y_cm_old(p,e)];
            end
            
            rho_start = sqrt((In_Guess(1)-PMT_r(topchs,1)).^2+(In_Guess(2)-PMT_r(topchs,2)).^2);
                        
            % To solve the problem that sometimes the initial value is just wrong
            
            [val pos] = max(PEAK_AREA);
            if rho_start(pos) > 10
                if (pos == 61) pos = 121; end
                In_Guess(1:2) = PMT_r(pos,1:2)+[2*rand-1, 2*rand-1];
                rho_start = sqrt((In_Guess(1)-PMT_r(topchs,1)).^2+(In_Guess(2)-PMT_r(topchs,2)).^2);
            end                
          
            %% Step 
            %% Selection of the PMTs involved in the fit - The number of PMTs does not change in each minimization
            %%

            cutEv = rho_start < rec_set.maximum_distance_PMTevent & rec_set.PMTS_To_Use(topchs)' == 1;
            DISTANCE = rec_set.maximum_distance_PMTevent;
            while sum(cutEv)<rec_set.PMT_MinNum & DISTANCE <31
                DISTANCE = DISTANCE + 2.5;
                cutEv = rho_start < DISTANCE  & rec_set.PMTS_To_Use(topchs)' == 1;
            end
            if sum(cutEv) < 3
                cutEv = rec_set.PMTS_To_Use(topchs)' == 1;
            end
                
            %% To increase the speed we create these two variables
            
            data_not_eq = double(dp.peak_area_phe(p,topchs,e));
            
            if dp.pulse_area_phe(p,e) > 5000
                val_minimum_data = 2.5;
            else
                val_minimum_data = 1;
            end
            Chi2_divisor = (max(abs(data_not_eq(cutEv)'), val_minimum_data).*(1+lrf_iq.sphe(cutEv).^2)');

            %% Step 
            %% Actual Minimization (Three different methods encoded)
            %%
            
            %% Method I Chi 2 

            if PULSE_AREA > 10 & PULSE_AREA > rec_set.MLM_maxphe & sum(cutEv) > 2 & PULSE_AREA > rec_set.PhotonCount_ML
                if numel(In_Guess) > 2, dp.reconstructed(p, e) = 5; Method_Of_Reconstruction = 'recE'; 
                else dp.reconstructed(p, e) = 1; Method_Of_Reconstruction = 'chi2'; end;
                [p_min, fval, exitflag, output] = Min_Handle(@(p_min) chi2_Handle(p_min, PEAK_AREA, Chi2_divisor, PULSE_AREA, cutEv, lrfmat),double(In_Guess), minim_opt);
                
            %% Method II Maximum Likelihood
            elseif PULSE_AREA <= rec_set.MLM_maxphe & sum(cutEv) > 2 & PULSE_AREA > rec_set.PhotonCount_ML
                dp.reconstructed(p, e) = 2; %% Meaning the pulse was  reconstructed  using the Maximum Likelihood method
                Method_Of_Reconstruction = 'MLM';
                
                PEAK_AREA_MLM = PEAK_AREA.*lrf_iq.QE.QE_Values(topchs);
                % Number of points calculated = num_values*2 +1
                num_values = 12;
                ML_Methode.PEAK_AREA_CutEv = round(max(repmat(PEAK_AREA_MLM(cutEv), [2*num_values+1 1]) + repmat(-num_values:1:num_values, [sum(cutEv) 1])', -1));
                
                ML_Methode.GaussianProb = Gaussian_MLM(repmat(PEAK_AREA_MLM(cutEv), [2*num_values+1 1]), repmat(sphes(cutEv), [2*num_values+1 1]), max(ML_Methode.PEAK_AREA_CutEv, 0), repmat(lrf_iq.sigma_zero.^2,  [2*num_values+1 sum(cutEv)]));
                ML_Methode.GaussianProb(isinf(ML_Methode.GaussianProb) | isnan(ML_Methode.GaussianProb)) = 0;
                ML_Methode.GaussianProb(ML_Methode.PEAK_AREA_CutEv<= 0) = 0;
                
                ML_Methode.GaussianProb(repmat(PEAK_AREA_MLM(cutEv), [2*num_values+1 1]) == 0) = 0;
                ML_Methode.GaussianProb(repmat(PEAK_AREA_MLM(cutEv), [2*num_values+1 1]) == 0 & ML_Methode.PEAK_AREA_CutEv == 0) = 1;
                                
                % This is the point of comparison. It needs some discussion but I am not sure  how this works
                
                clear Factor_Normalizacao
                Factor_Normalizacao(1:numel(PEAK_AREA_MLM(cutEv))) = 1;
                
                ML_Methode.Normalization_Factor = 1;
                ML_Methode.InitalPoint = sum(2*log(sum(poisspdf(ML_Methode.PEAK_AREA_CutEv,  repmat(PEAK_AREA_MLM(cutEv), [2*num_values+1 1])).*ML_Methode.GaussianProb, 1)./Factor_Normalizacao))/ML_Methode.Normalization_Factor;
                [p_min, fval, exitflag, output] = Min_Handle(@(p_min) MLM_Handle(p_min, PEAK_AREA, ML_Methode, PULSE_AREA, cutEv, lrfmat),double(In_Guess), minim_opt);
                
            %% Method III Photon Counting
            elseif PULSE_AREA <= rec_set.PhotonCount_ML & sum(cutEv) > 2 

                Method_Of_Reconstruction = 'PhC';       
                dp.reconstructed(p, e) = 3; %% Meaning the pulse was  reconstructed  using the Maximum Likelihood method with the number of spikes method
                spike_num_top = double(dp.spike_count(p, topchs,e));
                if max(spike_num_top) < 2
                    In_Guess = mean([PMT_r(find(dp.spike_count(p, :,e)==max(double(dp.spike_count(p, topchs,e))) & PMT_r(:,3)' > 0) , 1) PMT_r(find(dp.spike_count(p, :,e)==max(double(dp.spike_count(p, topchs,e))) & PMT_r(:,3)' > 0), 2)], 1) + normrnd(0,0.1,1, 2);
                end
                spike_num = squeeze(spike_num_top(cutEv)); %% spike numbers
                peak_areas = PEAK_AREA(cutEv);
                SPIKE_SUM = sum(dp.spike_count(p,:,e));

                lrfmat.Width = pulse_width(p,e); %% pile-up cut
                lrfmat.Cut_Pile_Up = PEAK_AREA(cutEv)< min((lrfmat.Width/(10*rec_set.resolution_phe)), 20); %% Selecionar of PMTs onde se usa a contagem vs  contagem normal                   
               
                if SPIKE_SUM > 1.35*PULSE_AREA
                    if rec_set.verbose > 1, disp(sprintf('\n Problems with ringing spikes: %.2f and areas %.2f\n', SPIKE_SUM, PULSE_AREA)); end
                    dp.reconstructed(p, e) = 4;
                    Method_Of_Reconstruction = 'PhC-Rin';       
                    lrfmat.Cut_Pile_Up = (lrfmat.Cut_Pile_Up*0)==0;
                end
                peak_areasc = peak_areas.*lrf_iq.QE.QE_Values(cutEv);
                lrfmat.max_prob = sum(-2*spike_num(lrfmat.Cut_Pile_Up).*(1-log(max(spike_num(lrfmat.Cut_Pile_Up),1)))) + sum(-2*peak_areasc(~lrfmat.Cut_Pile_Up).*(1-log(max(peak_areasc(~lrfmat.Cut_Pile_Up),0.125))));

                [p_min, fval, exitflag, output] = Min_Handle(@(p_min) PC_Handle(p_min, peak_areas, spike_num, [PULSE_AREA PULSE_AREA], cutEv, lrfmat),double(In_Guess), minim_opt);
                %% Error Mode  
            else
                Method_Of_Reconstruction = 'None';
                dp.reconstructed(p, e) = 5;
                if rec_set.verbose > 1 & sum(cutEv)<3
                    disp('****************************')
                    disp(sprintf('Not sufficient PMTs to be Able to reconstruct the pulse'));
                    disp('****************************')
                end
                p_min = In_Guess.*0-400;
                fval = -400;
                exitflag = -400;
                output.iterations = -400;
            end

            %% Step 
            %% Value Attribution
            %%
            
            
            % Estimative of the areas
            if rec_set.get_estimated_peak_areas == 1
                
                if PULSE_AREA > rec_set.min_num_of_phe 
                    
                    Chi2_divisor_N = (max(abs(data_not_eq'), val_minimum_data).*(1+lrf_iq.sphe([1:60 121]).^2)');

                    [pulse_dir pulse_ref ENERGY_MINIMIZED] = peakval_Matricial(p_min, PEAK_AREA, Chi2_divisor_N, PULSE_AREA, 1:61, lrfmat);
                    pulse_dir(122) = 0; pulse_ref(122) = 0;
                    pulse_dir(121) = pulse_dir(61); pulse_ref(121) = pulse_ref(61);
                    pulse_dir(61) = 0; pulse_ref(61) = 0;

                else 
                    pulse_dir = ones(1, 122) -101; pulse_ref = ones(1, 122) -101; ENERGY_MINIMIZED = -101;
                end                            
                dp.peak_area_rec(p,:, e) = ENERGY_MINIMIZED*(pulse_dir + pulse_ref).*gains;
            end

            dp.x_cm(p,e) = p_min(1);
            dp.y_cm(p,e) = p_min(2);
            dp.chi2(p,e) = fval;
            dp.rec_dof(p,e) = sum(cutEv) - numel(p_min);
            dp.s2_rec(p,e) = ENERGY_MINIMIZED;
            dp.minimization_flag(p,e) = exitflag;
            
            %% Step 
            %% Compute chi2 map to study the minimization
            %%
            
            if isfield(rec_set, 'compute_map')
                if rec_set.compute_map == 1
                    chi2_map = zeros(161, 161);
                    xxx = zeros(161, 161);
                    yyy = zeros(161, 161);
                    
                    for X_n = 1:161
                        for Y_n = 1:161
                            xxx(X_n, Y_n) = (round(p_min(1))-2)+0.025*X_n;
                            yyy(X_n, Y_n) = (round(p_min(2))-2)+0.025*Y_n;
                            if dp.reconstructed(p,e)==1;
                                chi2_map(X_n,Y_n) = chi2_Handle([xxx(X_n, Y_n) yyy(X_n, Y_n)], PEAK_AREA, Chi2_divisor, PULSE_AREA, cutEv, lrfmat);
                            elseif dp.reconstructed(p,e)==2;
                                chi2_map(X_n,Y_n) = MLM_Handle([xxx(X_n, Y_n) yyy(X_n, Y_n)], PEAK_AREA, ML_Methode, PULSE_AREA, cutEv, lrfmat);
                            elseif dp.reconstructed(p,e)==3;
                                chi2_map(X_n,Y_n) = PC_Handle([xxx(X_n, Y_n) yyy(X_n, Y_n)], peak_areas, spike_num, [PULSE_AREA PULSE_AREA], cutEv, lrfmat);
                            elseif dp.reconstructed(p,e)==4;
                                chi2_map(X_n,Y_n) = PC_Handle([xxx(X_n, Y_n) yyy(X_n, Y_n)], peak_areas, spike_num, [PULSE_AREA PULSE_AREA], cutEv, lrfmat);
                            end
                        end
                    end
                    chi2_map(imag(chi2_map)~=0) = real(max(chi2_map(imag(chi2_map)==0)));
                    
                    if 0%dp.event_number(e) == 17291
                        
                        
                        
                        xxx(X_n, Y_n) = (round(p_min(1))-2)+0.025*X_n;
                        yyy(X_n, Y_n) = (round(p_min(2))-2)+0.025*Y_n;
                        
                        Raja = sqrt(p_min(1).^2 + p_min(2).^2);
                    
                        
                        phi = atan(abs(p_min(1))./abs(p_min(2)));
                        if isnan(phi)
                            phi=0;
                        end
                        WALL = lrfmat.disphi( round(phi*10000/pi+1));
                        dis_wall = (WALL - sqrt(p_min(1).^2+p_min(2).^2))*24.5/WALL;
                        
                        
                        m = 1;
                        for Tanja = (-1.170):0.001:(-0.910)
                            [chis(m) chita(m,1:61) Dir(m, 1:61) Ref(m, 1:61)] = PC_Handle([sin(Tanja)*Raja cos(Tanja)*Raja], peak_areas, spike_num, [PULSE_AREA SPIKE_SUM], cutEv, lrfmat);
                            phis(m) = Tanja*180/pi;
                            m = m + 1;
                        end
                        
                        chita_tot = 0
                        for PMT = 1:60
                            plot(phis, chita_tot, 'ko')
                            PMT
                            pause
                        end
                        chis
                    end

                    
                end
            end 
            
            %% Step 
            %% Search for Double Scatters
            %% This is a new module that is under developement           
            
            if isfield(rec_set, 'search_for_double_scatterer')
                if rec_set.search_for_double_scatterer == 1
                    PEAKQ = squeeze(dp.peak_area_eq(p,:,e));
                    PEAKQ([61:120 122]) = 0;
                    PEAK_AREA_C = PEAKQ-ENERGY_MINIMIZED.*(pulse_dir+pulse_ref);
                    %disp(sprintf('--------------------------------------------------------------------------------'));
                    %disp(sprintf('%s %2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%4.2f', 'I', dp.x_cm_1st(e), dp.y_cm_1st(e), dp.x_cm_2nd(e), dp.y_cm_2nd(e),1, fval));
                    fval = 10000000000000000000;
                    
                    PEAK_AREA_C(find(PEAKQ==max(PEAKQ))) = 0;
                    for s = 1:4
                        In_Guess = [PMT_r(find(PEAKQ==max(PEAKQ)), ...
                                          1) PMT_r(find(PEAKQ==max(PEAKQ)), 2) PMT_r(find(PEAK_AREA_C==max(PEAK_AREA_C)), 1) PMT_r(find(PEAK_AREA_C==max(PEAK_AREA_C)), 2) min(max(PEAKQ)./(max(PEAKQ)+max(PEAK_AREA_C)), 0.75)];
                        p_min_s = In_Guess;
                        Method_Of_Reconstruction = 'MLM';
                        PEAK_AREA_MLM = PEAK_AREA.*lrf_iq.QE.QE_Values(topchs);
                        % Number of points calculated = num_values*2 +1
                        num_values = 12;
                        cutEv = ones(1,61);
                        cutEv(~rec_set.PMTS_To_Use(topchs)) = 0;
                        cutEv = cutEv==1;
                        ML_Methode.PEAK_AREA_CutEv = round(max(repmat(PEAK_AREA_MLM(cutEv), [2*num_values+1 1]) + repmat(-num_values:1:num_values, [sum(cutEv) 1])', -1));
                        ML_Methode.GaussianProb = Gaussian_MLM(repmat(PEAK_AREA_MLM(cutEv), [2*num_values+1 1]), repmat(sphes(cutEv), [2*num_values+1 1]), max(ML_Methode.PEAK_AREA_CutEv, 0), repmat(lrf_iq.sigma_zero.^2,  [2*num_values+1 sum(cutEv)]));
                        ML_Methode.GaussianProb(isinf(ML_Methode.GaussianProb) | isnan(ML_Methode.GaussianProb)) = 0;
                        ML_Methode.GaussianProb(ML_Methode.PEAK_AREA_CutEv<= 0) = 0;
                        
                        if rec_set.cutzerophotonsdata==0
                            ML_Methode.GaussianProb(repmat(PEAK_AREA_MLM(cutEv), [2*num_values+1 1]) == 0) = 0;
                            ML_Methode.GaussianProb(repmat(PEAK_AREA_MLM(cutEv), [2*num_values+1 1]) == 0 & ML_Methode.PEAK_AREA_CutEv == 0) = 1;
                        end
                        
                        clear Factor_Normalizacao
                        Factor_Normalizacao(1:numel(PEAK_AREA_MLM(cutEv))) = 1;
                        
                        % this is provisory. We have to introduce the pedestal
                        
                        %ML_Methode.Normalization_Factor = max(sum(cutEv)-2, 0.0000001);
                        ML_Methode.Normalization_Factor = 1;
                        ML_Methode.InitalPoint = sum(2*log(sum(poisspdf(ML_Methode.PEAK_AREA_CutEv,  repmat(PEAK_AREA_MLM(cutEv), [2*num_values+1 1])).*ML_Methode.GaussianProb, 1)./Factor_Normalizacao))/ML_Methode.Normalization_Factor;
                        
                        [p_min_s, fval_A, exitflag, output] = Min_Handle(@(p_min_s) MLM_Doubles(p_min_s, PEAK_AREA, ML_Methode, PULSE_AREA, cutEv, lrfmat),In_Guess, minim_opt);
                        if fval_A<fval
                            p_min = p_min_s;
                            fval = fval_A;
                        end
                        
                        PEAK_AREA_C(find(PEAK_AREA_C==max(PEAK_AREA_C))) = 0;
                        disp(sprintf(['%d %d %2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%4.2f'], e, s, p_min_s(1), p_min_s(2), p_min_s(3), p_min_s(4), p_min_s(5), fval_A));
                        
                    end
                    disp(sprintf(['%s %d %2.2f\t%2.2f\t%2.2f\t%2.2f\t%2.2f\t%4.2f'], 'F', s, p_min(1), p_min(2), p_min(3), p_min(4), p_min(5), fval));
                    disp(sprintf('--------------------------------------------------------------------------------'));
                    dp.x_mul(p, 1, e) = p_min(1);
                    dp.y_mul(p, 1, e) = p_min(2);
                    dp.x_mul(p, 2, e) = p_min(3);
                    dp.y_mul(p, 2, e) = p_min(4);
                    dp.chi2_mul(p, e) = p_min(5);
                        
                    A(1, 1) = sqrt((p_min(1)-dp.x_cm_1st(e)).^2 + (p_min(2)-dp.y_cm_1st(e)).^2);
                    A(2, 2) = sqrt((p_min(3)-dp.x_cm_2nd(e)).^2 + (p_min(4)-dp.y_cm_2nd(e)).^2);
                    A(2, 1) = sqrt((p_min(3)-dp.x_cm_1st(e)).^2 + (p_min(4)-dp.y_cm_1st(e)).^2);
                    A(1, 2) = sqrt((p_min(1)-dp.x_cm_2nd(e)).^2 + (p_min(2)-dp.y_cm_2nd(e)).^2);
                    %disp(sprintf('--------------------------------------------------------------------------------'));
                    %disp(A)
                    %disp(min(A(1, 1)+A(2, 2), A(2, 1)+A(1, 2)))
                    
                    if 0
                        figure(5);
                        clf
                        myfigview, grid on, grid minor, hold on
                        title('Understand the doubles');                
                        PEAKQ = squeeze(dp.peak_area_eq(p,:,e));
                        PlotTopArray(0.7)
                        plot(dp.x_cm_1st(e), dp.y_cm_1st(e), 'bo', 'MarkerFaceColor', 'b');
                        %plot(p_min_alte(1), p_min_alte(2), 'mv', 'MarkerFaceColor', 'm');
                        plot(dp.x_cm_2nd(e), dp.y_cm_2nd(e), 'ro', 'MarkerFaceColor', 'r');
                        plot(PMT_r(find(PEAKQ==max(PEAKQ)), 1), PMT_r(find(PEAKQ==max(PEAKQ)), 2), 'b+', 'MarkerFaceColor', 'b');
                        plot(PMT_r(find(PEAK_AREA_C==max(PEAK_AREA_C)), 1), PMT_r(find(PEAK_AREA_C==max(PEAK_AREA_C)), 2), 'r+', 'MarkerFaceColor', 'r');
                        plot(p_min(1), p_min(2), 'bs', 'MarkerFaceColor', 'b');
                        plot(p_min(3), p_min(4), 'rs', 'MarkerFaceColor', 'r');
                        for PMT = [1:60 121]
                            if PEAKQ(PMT)~=0
                                text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)+1.00, sprintf('%.1f', PEAKQ(PMT)), 'color', 'b');
                                text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)-1.00, sprintf('%.1f', PEAK_AREA_C(PMT)), 'color', 'g');
                                text(20,22,sprintf('CHI 2 %4.0f', fval));
                                text(20,20,sprintf('PMT A %d', find(PEAK_AREA_C==max(PEAK_AREA_C))));
                                text(20,18,sprintf('PMT B %d', find(PEAKQ==max(PEAKQ))));
                            end
                        end
                        pause
                    end
                
                end
            end 
            
            %%
            %% Compute standard deviations
            %%

            if rec_set.compute_sd == 1 & dp.pulse_area_phe(p,e) > 33 & dp.reconstructed(p,e)< 5

                Phi_O = atan2(p_min(2), p_min(1));
                R_O = sqrt(p_min(2).^2+ p_min(1).^2);
                if isfield(rec_set, 'uncertainty_step'), increment_r = rec_set.uncertainty_step./sqrt(dp.pulse_area_phe(p,e));
                else, increment_r = 0.5/sqrt(dp.pulse_area_phe(p,e)); end
                increment_phir = increment_r;
                chi_max = 2.30;

                R_S = R_O; chi_r = fval;
                
                while chi_r < fval+chi_max & R_S < 26,
                    R_S = R_S + increment_r;
                    if dp.reconstructed(p,e)==1;
                        chi_r =  chi2_Handle([R_S*cos(Phi_O) R_S*sin(Phi_O)], PEAK_AREA, Chi2_divisor, PULSE_AREA, cutEv, lrfmat);
                    elseif dp.reconstructed(p,e)==2;
                        chi_r = MLM_Handle([R_S*cos(Phi_O) R_S*sin(Phi_O)], PEAK_AREA, ML_Methode, PULSE_AREA, cutEv, lrfmat);
                    elseif dp.reconstructed(p,e)==3;
                        chi_r = PC_Handle([R_S*cos(Phi_O) R_S*sin(Phi_O)], peak_areas,spike_num, [PULSE_AREA PULSE_AREA], cutEv, lrfmat);
                    elseif dp.reconstructed(p,e)==4;
                        chi_r = PC_Handle([R_S*cos(Phi_O) R_S*sin(Phi_O)], peak_areas,spike_num, [PULSE_AREA PULSE_AREA], cutEv, lrfmat);
                    end
                end
                
                R_I = R_O; chi_r = fval;
                while chi_r < fval+chi_max & R_I > 0,
                    R_I = R_I - increment_r;
                    if dp.reconstructed(p,e)<2;
                        chi_r =  chi2_Handle([R_I*cos(Phi_O) R_I*sin(Phi_O)], PEAK_AREA, Chi2_divisor, PULSE_AREA, cutEv, lrfmat);
                    elseif dp.reconstructed(p,e)==2;
                        chi_r = MLM_Handle([R_I*cos(Phi_O) R_I*sin(Phi_O)], PEAK_AREA, ML_Methode, PULSE_AREA, cutEv, lrfmat);
                    elseif dp.reconstructed(p,e)==3;
                        chi_r = PC_Handle([R_I*cos(Phi_O) R_I*sin(Phi_O)], peak_areas,spike_num, [PULSE_AREA PULSE_AREA], cutEv, lrfmat);
                    elseif dp.reconstructed(p,e)==4;
                        chi_r = PC_Handle([R_I*cos(Phi_O) R_I*sin(Phi_O)], peak_areas,spike_num, [PULSE_AREA PULSE_AREA], cutEv, lrfmat);
                    end
                end
                
                if R_O > 1,
                    Phi_S = Phi_O*R_O; chi_r = fval;
                    while chi_r < fval+chi_max  & Phi_S < 4*R_O,
                        Phi_S = Phi_S + increment_phir;
                        if dp.reconstructed(p,e)<2;
                            chi_r = chi2_Handle([R_O*cos(Phi_S/R_O) R_O*sin(Phi_S/R_O)], PEAK_AREA, Chi2_divisor, PULSE_AREA, cutEv, lrfmat);
                        elseif dp.reconstructed(p,e)==2;
                            chi_r = MLM_Handle([R_O*cos(Phi_S/R_O) R_O*sin(Phi_S/R_O)], PEAK_AREA, ML_Methode, PULSE_AREA, cutEv, lrfmat);
                        elseif dp.reconstructed(p,e)==3;
                            chi_r = PC_Handle([R_O*cos(Phi_S/R_O) R_O*sin(Phi_S/R_O)], peak_areas,spike_num, [PULSE_AREA PULSE_AREA], cutEv, lrfmat);
                        elseif dp.reconstructed(p,e)==4;
                            chi_r = PC_Handle([R_O*cos(Phi_S/R_O) R_O*sin(Phi_S/R_O)], peak_areas,spike_num, [PULSE_AREA PULSE_AREA], cutEv, lrfmat);
                        end
                    end
                    
                    Phi_I = Phi_O*R_O;
                    chi_r = fval;
                    while chi_r < fval+chi_max  & Phi_I > -4*R_O,
                        Phi_I = Phi_I - increment_phir;
                        if dp.reconstructed(p,e)<2;
                            chi_r = chi2_Handle([R_O*cos(Phi_I/R_O) R_O*sin(Phi_I/R_O)], PEAK_AREA, Chi2_divisor, PULSE_AREA, cutEv, lrfmat);
                        elseif dp.reconstructed(p,e)==2;
                            chi_r = MLM_Handle([R_O*cos(Phi_I/R_O) R_O*sin(Phi_I/R_O)], PEAK_AREA, ML_Methode, PULSE_AREA, cutEv, lrfmat);
                        elseif dp.reconstructed(p,e)==3;
                            chi_r = PC_Handle([R_O*cos(Phi_I/R_O) R_O*sin(Phi_I/R_O)], peak_areas,spike_num, [PULSE_AREA PULSE_AREA], cutEv, lrfmat);
                        elseif dp.reconstructed(p,e)==4;
                            chi_r = PC_Handle([R_O*cos(Phi_I/R_O) R_O*sin(Phi_I/R_O)], peak_areas,spike_num, [PULSE_AREA PULSE_AREA], cutEv, lrfmat);
                        end
                    end
                else
                    Phi_I = 0;
                    Phi_S = 0;
                end
                    
                dp.sd_radius_inf(p,e) = abs(R_I-R_O);
                dp.sd_radius_sup(p,e) = abs(R_S-R_O);
                dp.sd_phiXR(p,e) = max(abs(Phi_I-Phi_O*R_O), abs(Phi_S-Phi_O*R_O));
            else    
            	R_O = 1; Phi_O = 0;% they are not defined, this is patch
            	R_I = R_O;
                R_S = R_O;
                Phi_I = 0;
                Phi_S = 0;
            end
            if rec_set.verbose > 1 
                The_Phi = atan2(dp.y_cm_old(p,e), dp.x_cm_old(p,e));
                disp(sprintf(['%d\t%d\t%s\t%.2f\t%1.f\t%1.f\t%1.f\' ...
                              't%.2f\t%.2f\t%.2f\t%.2f\t%.2f\' ...
                              't%.2f\t%.2f\t%.2f\t%.3f\t%.3f\' ...
                              't%.3f\t%.3f\t%d\t%d\t%d\t%.2f\t%.2f'], ...
                             dp.event_number(e), p, Method_Of_Reconstruction,cputime-t_0, ...
                             squeeze(sum(dp.peak_area_phe(p,:, e), 2)), PULSE_AREA, ENERGY_MINIMIZED, ...
                             sqrt(dp.x_cm_old(p,e).^2+dp.y_cm_old(p,e).^2), sqrt(p_min(1)^2+ ...
                                                       p_min(2)^2), ...
                             The_Phi*180/(2*pi), atan2(p_min(2), ...
                                                       p_min(1))* ...
                             180/(2*pi), In_Guess(1), In_Guess(2), p_min(1), p_min(2),abs(R_I-R_O), abs(R_S-R_O), abs(Phi_I-Phi_O*R_O), abs(Phi_S-Phi_O*R_O), exitflag, output.iterations, sum(cutEv), chi2_old(p,e), fval));
            else

                %waitbar(e/NEVTS, 'wait');
            end
            
            if rec_set.tests

                my_fig = figure(122);
                clf
                Resul = get(my_fig, 'Position');
                %% This is Just Stupid
                %set(gcf,'Position',[Resul(1)*0.8 Resul(2) Resul(3)*2.25 Resul(4)*2.25]);
                if isfield(dp, 's1_xyz_phe')
                    suptitle(sprintf('File Name: %s -- Event: %d  -- Pulse Area %.2f phe -- S2S1Ratio: %.3f -- Drift time %.0f -- S2 Width %d\n\n ', dp.dataset(e), dp.event_number(e), dp.s1_xyz_phe(e), dp.log10s2bs1_xyz(e), dp.pulse_area_phe(p,e), pulse_width(p,e)));
                elseif ~isfield(dp, 'file_number') | ~isfield(dp, 'hft_t10l_samples') | ~isfield(dp, 's1_xyz_phe')
                    suptitle(sprintf('Event: %d Pulse: %d Pulse Area %.2f phe\n\n ', dp.event_number(e), p, dp.pulse_area_phe(p,e)));

                else
                    suptitle(sprintf(['File Number: %d -- Event: %d -- Pulse: %d -- Pulse Area %.2f phe -- Starting at: %d\n\n '], dp.file_number(e), dp.event_number(e), p, dp.pulse_area_phe(p,e), dp.hft_t10l_samples(p,e)));
                end
                    
                subplot(2,2,2);
                title('Predicted Vs. Original Result');
                myfigview, grid on, grid minor, hold on
                [peak_area_sorted pmts_sorted] =  sort(dp.peak_area_eq(p,topchs,e), 'descend');
                plot(1:61, ENERGY_MINIMIZED*(pulse_dir(topchs(pmts_sorted)) + pulse_ref(topchs(pmts_sorted))), 'rd', 'MarkerFaceColor', 'r')
                plot(1:61, dp.peak_area_eq(p,topchs(pmts_sorted),e), 'bs', 'MarkerFaceColor', 'b')
                plot(1:61, (2*dp.(pmt_sat_flag)(p,topchs(pmts_sorted),e)-1.0)*1.1*max(ENERGY_MINIMIZED*(pulse_dir(topchs(pmts_sorted)) + pulse_ref(topchs(pmts_sorted)))), 'go', 'MarkerFaceColor', 'g')
                plot(1:61, (2*(PEAK_AREA(pmts_sorted)' > rec_set.saturated_pmt_phe_limit(pmts_sorted)')-1.)*1.1*max(ENERGY_MINIMIZED*(pulse_dir(topchs(pmts_sorted)) + pulse_ref(topchs(pmts_sorted)))), 'mo', 'MarkerFaceColor', 'm')

                xlabel('From the Highest to the lowest #');
                ylabel(sprintf('Photons (phe)'));
                ylim([0 1.3*max(ENERGY_MINIMIZED*(pulse_dir(topchs) + pulse_ref(topchs)))])
                legend('Predicted Result', 'Original Result')
                
                subplot(2,2,1);
                myfigview, grid on, grid minor, hold on
                
                PEAKQ = squeeze(dp.peak_area_eq(p,:,e));
                EXPEC = ENERGY_MINIMIZED*(pulse_dir+pulse_ref);
                PEAKT = sum(PEAKQ);
                if isfield(dp, 'spike_count')
                    SPIKES = squeeze(dp.spike_count(p,:,e));
                end
                PlotTopArray(0.7)
                
                if rec_set.compute_sd
                    plot(dp.sd_radius_inf(p,e)*cos(0:(pi*0.001):2*pi)+p_min(1), dp.sd_radius_inf(p,e)*sin(0:(pi*0.001):2*pi)+p_min(2), 'g-', 'LineWidth',2);
                end
                plot(p_min(1), p_min(2), 'mo', 'MarkerFaceColor', 'm');
                plot(dp.x_cm_old(p,e), dp.y_cm_old(p,e), 'rs', 'MarkerFaceColor', 'r');

                for PMT = [1:60 121]
                    if ~isfield(dp, 'spike_count')
                        if PEAKT<20000
                            if PEAKQ(PMT)~=0
                                text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)+1, sprintf('%.1f', PEAKQ(PMT)), 'color', 'r');
                            end
                            text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)-1, sprintf('%.1f', EXPEC(PMT)), 'color', 'b');
                        else
                            if PEAKQ(PMT)~=0
                                text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)+1, sprintf('%.1f', PEAKQ(PMT)/1000.), 'color', 'r');
                            end
                            text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)-1, sprintf('%.1f', EXPEC(PMT)/1000.), 'color', 'b');
                        end
                    else
                        if SPIKES(PMT)~=0
                            text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)+1.25, sprintf('%.d', SPIKES(PMT)), 'color', 'k');
                            text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)-0.25, sprintf('%.1f', PEAKQ(PMT)), 'color', 'r');
                        end
                        if cutEv(min(PMT, 61))==1
                            text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)-2.25, sprintf('%.1f', EXPEC(PMT)), 'color', 'b');
                        else
                            text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)-2.25, sprintf('%.1f', EXPEC(PMT)), 'color', 'g');
                        end
                    end
                end
                
                %text(14,25, sprintf('Ev: %.0f', dp.event_number(e)))
                %text(14,22, sprintf('Phe: %.1f', sum(PEAKQ)))
                plot(20*cos(0:(pi*0.001):2*pi), 20*sin(0:(pi*0.001):2*pi), 'b-');
                %plot(17*cos(0:(pi*0.001):2*pi),
                %17*sin(0:(pi*0.001):2*pi), 'r-');
                
                if rec_set.compute_map == 1
                    %% This will plot the chi2 curve / Minimization curve
                    subplot(2,2,4);
                    title(sprintf('%s Map', '\chi^2_{red}'))
                    PlotTopArray(0.25, 1);
                    
                    myfigview
                    MaxValueOfChi2Map = min([min(chi2_map(1,:)) min(chi2_map(end,:)) min(chi2_map(:,1)) min(chi2_map(:,end))]);
                    C = contour3(xxx,yyy,min(chi2_map, MaxValueOfChi2Map),40);
                    %C = contour3(xxx,yyy,chi2_map,20);
                    colormap(jet);
                    axis([min(xxx(:)) max(xxx(:)) min(yyy(:)) max(yyy(:))])
                    plot(p_min(1), p_min(2), 'mo',  'MarkerFaceColor', 'm')
                    try
                        clabel(C);
                    end
                    xlabel('x (cm)')
                    ylabel('y (cm)')
                    
                end
                

                
                comparisonpeaks(1:122) = 0;
                comparisonpeaks(topchs) = dp.peak_area_phe(p,topchs,e);
                comparisonpeaks([61:120 122]) = ENERGY_MINIMIZED* ...
                    (pulse_dir(topchs) + pulse_ref(topchs));
                options.fig_off = 1;
                subplot(2,2,3);
                LUXHitPattern(comparisonpeaks, options);
                myfigview
                title(['Observed (top) versus expected (bottom) hit pattern\n\n'])
                %save_graphic(sprintf('%s_%d', dp.filename{e}, dp.event_number(e)), [15 10], 'eps')
                           
                
                if isfield(rec_set, 'compute_map_doubles')
                    if rec_set.compute_map_doubles == 1
                        %% This will plot the chi2 curve / Minimization curve
                        figure(1);
                        title(sprintf('%s Map', '\chi^2_{red}'))
                        PlotTopArray(0.25);

                        C = contour3(xxx_b, yyy_b, chi2_map_doubles_100,12);
                        
                        colormap(jet);
                        plot(p_min(1), p_min(2), 'mo',  'MarkerFaceColor', 'm')
                        try
                            clabel(C);
                        end
                        xlabel('x (cm)')
                        ylabel('y (cm)')
                    end
                end
                pause
            end
            
            %close(my_fig);
            
            if rec_set.verbose == 0
                PMTs = topchs;
                pmts=PMTs(cutEv);
            end
                        
                        %%
            %% Step  Compute Doubles
            %%
            if isfield(rec_set, 'test_double')
                if rec_set.test_double
                    %% Initial Guess Value
                    pos_1 = p_min(1:2);
                    rho_reject = sqrt((pos_1(1)-PMT_r(topchs,1)).^2+(pos_1(2)-PMT_r(topchs,2)).^2)<9 | (PMT_r(topchs,1).^2+PMT_r(topchs,2).^2)>400;
                    PEAK_AREA_C = PEAK_AREA-pulse_dir(topchs)-pulse_ref(topchs);
                    pos_2 = PMT_r(find(PEAK_AREA_C==max(PEAK_AREA(~rho_reject)-pulse_dir(~rho_reject)-pulse_ref(~rho_reject))), 1:2);
                    In_Guess = [pos_1(1) pos_1(2) pos_2(1) pos_2(2) 0.75]; %% Initial Guess 5 parameters are used
                    
                    
                    if dp.reconstructed(p, e) == 3; %% Meaning the pulse was  reconstructed  using the Maximum Likelihood method with the number of spikes method
                        cutEv(1:61) = rec_set.PMTS_To_Use(topchs)' == 1; %% ALL PMTS ARE BEING CONSIDERED
                        
                        spike_num = squeeze(spike_num_top(cutEv)); %% spike numbers
                        lrfmat.Width = pulse_width(p,e); %% pile-up cut
                        lrfmat.Cut_Pile_Up = PEAK_AREA(cutEv)< min((lrfmat.Width/(5*rec_set.resolution_phe)), 3); %% Selecionar of PMTs onde se usa a contagem vs normal contagem
                        lrfmat.degrees_of_freedom = (max(sum(cutEv)-6,1)); %% Degrees of freedom
                        
                        lrfmat.max_prob = sum(-2*spike_num(lrfmat.Cut_Pile_Up).*(1-log(max(spike_num(lrfmat.Cut_Pile_Up),1))))/lrfmat.degrees_of_freedom;
                        
                        data_not_eq = double(dp.peak_area_phe(p,topchs,e));
                        Chi2_divisor = max(sum(cutEv)-2, 1)*(max(abs(data_not_eq(cutEv)'),0.5)+(rec_set.chib*PULSE_AREA)^2);
                        lrfmat.Chi2_divisor = Chi2_divisor(~lrfmat.Cut_Pile_Up);
                        QE_Corr = lrfmat.lrf_iq.QE.QE_Values(cutEv);
                        lrfmat.QE_Corrections = QE_Corr(lrfmat.Cut_Pile_Up);
                        
                        [p_min_double, fval, exitflag, output] = Min_Handle(@(p_min_double) MLM_PhotonCount_Doubles(p_min_double, PEAK_AREA, spike_num, PULSE_AREA, cutEv, lrfmat), In_Guess, minim_opt);
                    end
                    predicted_pulses = peakval_Doubles(p_min_double, PEAK_AREA, Chi2_divisor, PULSE_AREA, 1:61, lrfmat);
                    
                    %% Now the Plots
                    
                    my_fig = figure(123);
                    Resul = get(my_fig, 'Position');
                    %% This is Just Stupid
                    set(gcf,'Position',[Resul(1)*0.8 Resul(2) Resul(3)*2.25 Resul(4)*2.25]);
                    %suptitle(sprintf('Checking Double Scatter Events ::: Event: %d ::: Pulse: %d ::: First %.2f phe ::: Second  %.2f phe \n\n ', dp.event_number(e), p, dp.pulse_area_phe(p,e)*p_min_double(5), dp.pulse_area_phe(p,e)*(1-p_min_double(5))));
                    
                    subplot(2,2,2);
                    title('Predicted Vs. Original Result')
                    myfigview, grid on, grid minor, hold on
                    [peak_area_sorted pmts_sorted] =  sort(dp.peak_area_eq(p,topchs,e), 'descend');
                    plot(1:61, predicted_pulses(topchs(pmts_sorted)), 'rd', 'MarkerFaceColor', 'r')
                    plot(1:61, dp.peak_area_eq(p,topchs(pmts_sorted),e), 'bs', 'MarkerFaceColor', 'b')
                    
                    xlabel('From the Highest to the lowest #');
                    ylabel(sprintf('Photons (phe)'));
                    ylim([0 1.3*max(ENERGY_MINIMIZED*(pulse_dir(topchs) + pulse_ref(topchs)))])
                    legend('Predicted Result', 'Original Result')
                    
                    subplot(2,2,1);
                    myfigview, grid on, grid minor, hold on
                    
                    PEAKQ = squeeze(dp.peak_area_eq(p,:,e));
                    EXPEC = predicted_pulses;
                    PEAKT = sum(PEAKQ);
                    SPIKES = squeeze(dp.spike_count(p,:,e));
                    
                    PlotTopArray(0.7)
                    
                    plot(p_min_double(1), p_min_double(2), 'mo', 'MarkerFaceColor', 'm');
                    plot(p_min_double(3), p_min_double(4), 'mo', 'MarkerFaceColor', 'm');
                    
                    for PMT = [1:60 121]
                        if dp.reconstructed(p,e)~=3
                            if PEAKT<20000
                                if PEAKQ(PMT)~=0
                                    text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)+1, sprintf('%.1f', PEAKQ(PMT)), 'color', 'r');
                                end
                                text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)-1, sprintf('%.1f', EXPEC(PMT)), 'color', 'b');
                            else
                                if PEAKQ(PMT)~=0
                                    text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)+1, sprintf('%.1f', PEAKQ(PMT)/1000.), 'color', 'r');
                                end
                                text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)-1, sprintf('%.1f', EXPEC(PMT)/1000.), 'color', 'b');
                            end
                        else
                            if SPIKES(PMT)~=0
                                text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)+1.25, sprintf('%.d', SPIKES(PMT)), 'color', 'k');
                                text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)-0.25, sprintf('%.1f', PEAKQ(PMT)), 'color', 'r');
                            end
                            text(PMT_r(PMT,1)-1.75, PMT_r(PMT,2)-2.25, sprintf('%.1f', EXPEC(PMT)), 'color', 'b');
                        end
                    end
                
                    %text(14,25, sprintf('Ev: %.0f', dp.event_number(e)))
                    %text(14,22, sprintf('Phe: %.1f', sum(PEAKQ)))
                    plot(20*cos(0:(pi*0.001):2*pi), 20*sin(0:(pi*0.001):2*pi), 'b-');
                    %plot(17*cos(0:(pi*0.001):2*pi),
                    %17*sin(0:(pi*0.001):2*pi), 'r-');
                    
                    comparisonpeaks(1:122) = 0;
                    comparisonpeaks(topchs) = dp.peak_area_phe(p,topchs,e);
                    comparisonpeaks([61:120 122]) = predicted_pulses(topchs);
                    options.fig_off = 1;
                    subplot(2,2,3);
                    LUXHitPattern(comparisonpeaks, options);
                    myfigview
                    title(['Observed (top) versus expected (bottom) hit pattern\n\n'])

                end
            end
        end
    end
end    
if isfield(dp, 'pmt_saturation_flag_prov')
    dp = rmfield(dp, 'pmt_saturation_flag_prov');
end

dtsetReconstruct = dp;

if rec_set.verbose > 0
    fprintf('Done!\n')
    sprintf('Position reconstruction of %d pulses in %.2f s\n', EVTS, toc);
end
