function dp = Corrections_PositionCorrection_Function(dp, table_corrections) 

% dp = Corrections_PositionCorrection_Function(dp, table_of_correction)
%
% Inputs:
% dp - The data structure with the following RQs (s1s2_pairing, z_drift_samples, pulse_classification, x_cm, y_cm)
% The following inputs are optional:
% factor_of_correction - The matrix with the factors of correction that is applied to the positions 
% In the case that these inputs are not used the function will look for the default settings which should be 
% in the same folder of this function.
%
%Outputs:
% dp - Output datastructure with two new RQs with the corrected positions (dp.x_corrected and dp.y_corrected)
% When the final position has an absurd output that output should
% be understood as an error code.
% -100 Initial position of -100
% -200 Are placed clearly outside the chamber
% -300 The drift time is not between 0 and maximum drift time
% -400 The S2 does not have the correspondent pair
% -500 It's not classified as an S2 (dp.pulse_classification~=2)
%
%EXAMPLE OF USAGE:    
%
% dp = PositionCorrection_Function(dp, factor_of_correction)    
%
%Versioning:
% v1.00 - % 20130725 Created by Claudio
% v1.10 - % 20130729 Implemented a Translation along XY and the code is adapted to different sizes of the factor_of_correction table.
% v1.10.1 - % 20130731 Some Corrections related with the name of the function. Now the user does not need to provide the IQs for the corrections.
% v1.20 - % 20130903 Modified to include the reconstruction from
% TAXY. The calling of the function was modified
% v1.20.1 % 20130903 Patch 
% v1.20.2 % 20130910 Removed some printing information
% v1.20.3 % 20130912 Corrections to the sigma variable introduced
% v1.30.1 % 20140107 Small correction to allow more sections along
% the drift time
% error codes modified PAT 171228
% error code modified to exclude only s1_like_class5 and keep s2like 180205
%   180227 PAT added correction for single event files (not sure why needed but it is.    
%   

if nargin < 2  %probably don't need this? - this checks that there are 2 input args and finds correction structs (both of the structs at .mat files) 
    myname = 'Corrections_PositionCorrection';
    position_correction_path = which(myname);
    IQs = dir([position_correction_path(1:(end-numel(myname)-2)) '*Mercury.mat']);  %IQs holds the names of both correction structs
    if numel(IQs)<1
        disp(sprintf('***\nCorrections_PositionCorrection: MAJOR ERROR - No IQs for this function.\n Please check if settings (a mat file with the factor of corrections and translations) is in this folder %s\n***', position_correction_path(1:(end-numel(myname)-2))));
        return
    else
        table_corrections = load(IQs(end).name); %gets name of last element in IQs
        disp(sprintf('***\nCorrections_PositionCorrection: You have not indicated the IQs for this module.\n Loading the default IQ with the name %s and within the following path\n %s\n***', IQs(end).name, position_correction_path(1:(end-numel(myname)-2))));
    end
end

if ~isfield(table_corrections.cor_xy_iq, 'Translation') %checks properies of the table corrections struct, which was either loaded in the above if then statement or provided as input
    table_corrections.cor_xy_iq.Translation = zeros(2, size(table_corrections.cor_xy_iq.factor_of_correction, 1)); % if 'Translation' isn't there, it adds it
end

factor_of_correction = table_corrections.cor_xy_iq.factor_of_correction; %be careful, this is huge
Translation = table_corrections.cor_xy_iq.Translation;

if isfield(table_corrections.cor_xy_iq.global, 'algorithm_reconstruction') 
    algt = table_corrections.cor_xy_iq.global.algorithm_reconstruction;
else
    algt = 'Mercury'; %defaults Mercury instead of TAXY
end
    
%dp = which_class5_type(dp);  %180205 PAT added 

if strcmp(algt, 'TAXY') %probably don't need this since we don't use TAXY, and can just use what is in the 'else' part
    
    X_Unc = double(dp.x_cm_tmplt);
    Y_Unc = double(dp.y_cm_tmplt);
    Sigma_Unc = double(dp.xy_sigma_cm);
else
    X_Unc = double(dp.x_cm);
    Y_Unc = double(dp.y_cm);
    Sigma_Unc = 0; %% It is not working for Mercury right now
end
    
if isfield(table_corrections.cor_xy_iq.global, 'max_drift_time') % set max drift time to 33000 if not already there
    max_drift_time = table_corrections.cor_xy_iq.global.max_drift_time;
else
    max_drift_time = 33000;
end

if isfield(table_corrections.cor_xy_iq.global, 'min_drift_time') % sets min drift to 400 unless already there
    min_drift_time = table_corrections.cor_xy_iq.global.min_drift_time;
else
    min_drift_time = 400;
end


if isfield(table_corrections.cor_xy_iq.global, 'Contraction') %checks for 'contraction' - I don't know what that is
    Contraction = table_corrections.cor_xy_iq.global.Contraction;
else
    Contraction = 0;
end

[N_DT N_PHI N_R2] = size(factor_of_correction); %gets dimensions of correction struct - drift time (segments of 1000 samples), angle, rad^2

[NPULS, NEVTS] = size(X_Unc); % determins number of pulses and events in the dataset
X_Cor = ones(NPULS, NEVTS).*100;
Y_Cor = ones(NPULS, NEVTS).*100;
Sigma_Cor = ones(NPULS, NEVTS).*100;

if isfield(dp, 'dtime_us')
    dp.z_drift_samples = dp.dtime_us.*100;
end 
if ~isfield(dp, 's1s2_pairing') & NPULS == 2 %structure change if no s1d2 pairing and only 2 pulses, probably don't need
    dp.s1s2_pairing = ones(NPULS, NEVTS);
end



if ~isfield(dp,'s1s2_pairing') | ~isfield(dp,'z_drift_samples') | ~isfield(dp,'pulse_classification') | ~isfield(dp,'x_cm') |  ~isfield(dp,'y_cm') %check that required fields are present
    disp(sprintf('Corrections_PositionCorrection: You need the following RQs: s1s2_pairing, z_drift_samples, pulse_classification, x_cm, y_cm'))
end


%%
%% Test if we have the RQs needed for this correction
%%
  
% Drift time

DT = ceil(double(dp.z_drift_samples)./round(max_drift_time/N_DT)); %looks at drift time, creates drift time bins where 1=1000 samples, rounds that up. 
DT(~inrange(DT, 0.5, N_DT+0.5)) = 1; % checks that these in the range 0.5 to 32.5, otherwise assigns it value of 1, this will then be used as index for translation. I think that the translation matrix in r3 is zero so this might not be important
STEP_DT = round(max_drift_time/(N_DT*100));
%% 
%% perform the Translation in XY
%%

X_T = squeeze(Translation(1,:));
x_cm_trans = X_Unc+X_T(DT);

Y_T = squeeze(Translation(2,:));
y_cm_trans = Y_Unc+Y_T(DT);

if length(dp.luxstamp_samples) == 1 % ie there is one event in the file
    x_cm_trans = x_cm_trans(:,1);
    y_cm_trans = y_cm_trans(:,1); 
    %for reasons beyond my conception of logic, the lines for y_cm_trans
    %and x_cm_trans ADD nPulse by 1 array to another nPulse by 1 array and
    %make it an nPulse by nPulse array. PAT 190227
end

%%
%% create the azimuthal variables
%%

phi = atan2(y_cm_trans, x_cm_trans).*180/pi;
radius2 = x_cm_trans.^2+y_cm_trans.^2;

%%
%% Get the corrective elements
%%


% Azimuthal angle
PHI = round((phi+180)./round(360/N_PHI));
PHI(PHI==0) = 360; %% This is an error catch for TAXY, thought it is useful also for Mercury

%Radius - this is a bit confusing and not easy to be mantained
Phi_Rec_pr = atan(y_cm_trans./x_cm_trans); %polar conversion
Phi_Rec_pr(isnan(Phi_Rec_pr)) = 0; % Modified for TAXY  - check to see that there are no infs that came in the polar conversion
PhiCut_pr = inrange(abs(Phi_Rec_pr), pi*1./6, pi/2);  % range cut
WALL = (24.5*cos(pi/12)./cos(abs(abs(abs(Phi_Rec_pr)-PhiCut_pr*pi/3)-pi/12)));  %I really don't know, wall positioning? but it is based off of angle
RATIO2 = (24.5./WALL).^2; 
if Contraction>0   %I think the standard contraction is 0.5, so we need to keep this
    SHRINK = 600./(800-Contraction.*double(STEP_DT*DT));
else
    SHRINK = 600./(720-double(DT*10));
end

RADI = ceil((double(radius2).*double(RATIO2).*double(SHRINK))); 
RAD = RADI;
RAD(RADI>600) = 1;
RAD(RADI==0) = 1;

cor_factor = ones(NPULS, NEVTS);

for p = 1:NPULS
    for e = 1:NEVTS
        cor_factor(p,e) = factor_of_correction(DT(p,e), PHI(p,e), RAD(p,e)); % finds correction factor in 3 space of a postion
    end
end

X_Cor = x_cm_trans./cor_factor;
Y_Cor = y_cm_trans./cor_factor;
Sigma_Cor = Sigma_Unc./cor_factor;

%% Elements outside the chamber or zero
X_Cor(RADI>600 | RADI == 0) = -200;
Y_Cor(RADI>600 | RADI == 0) = -200;
Sigma_Cor(RADI>600 | RADI == 0) = -200;


%% Elements outside the drift time between min_drift_time and
%% max_drift_time as defined by the look-up table
X_Cor(~inrange(dp.z_drift_samples, min_drift_time, max_drift_time)) = -300;
Y_Cor(~inrange(dp.z_drift_samples, min_drift_time, max_drift_time)) = -300;
Sigma_Cor(~inrange(dp.z_drift_samples, 0, max_drift_time)) = -300;

%% Cuts - previously badly reconstructed elements
X_Cor(X_Unc==-100) = -100;
Y_Cor(Y_Unc==-100) = -100;
Sigma_Cor(Y_Unc==-100) = -100;
%{
this section removed by PAT 171228

%% S2 signals without an S1 Pairing
X_Cor(dp.s1s2_pairing~=1) = -400;
Y_Cor(dp.s1s2_pairing~=1) = -400;
Sigma_Cor(dp.s1s2_pairing~=1) = -400;

%% Elements with the wrong pulse classification pulse_classification
X_Cor(dp.pulse_classification~=2) = -500;
Y_Cor(dp.pulse_classification~=2) = -500;
Sigma_Cor(dp.pulse_classification~=2) = -500;
%}

%% Elements with the wrong pulse classification pulse_classification  % this modified version was added by PAT 171228
% this checks allows for class 2 and 4 (s2 and SE) to get the correction
% while eliminating others. Note that this allows for -100 error code for
% no pulse, which seemed to be the intent of the originial program but was
% implemented incorrectly
X_Cor(dp.pulse_classification==1 | dp.pulse_classification==3 | dp.s1_like_class5) = -500; %changed to s1_like_class5 condition 180205
Y_Cor(dp.pulse_classification==1 | dp.pulse_classification==3 | dp.s1_like_class5) = -500;
Sigma_Cor(dp.pulse_classification==1 | dp.pulse_classification==3 | dp.s1_like_class5) = -500;

if strcmp(algt, 'TAXY')  %probably don't need taxy part
    dp.x_tmplt_corrected = X_Cor;
    dp.y_tmplt_corrected = Y_Cor;
    dp.xy_sigma_corrected = Sigma_Cor;
else
    dp.x_corrected = X_Cor;
    dp.y_corrected = Y_Cor;
end
