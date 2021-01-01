function status = dp_framework_split_sim_orig(filename_evt , num_in_each)

addpath(genpath(['/home/paul/LUXCode_original/Stable_Releases/DataProcessingStableReleases/v2.0/DataProcessing/']));
pathbase= '/home/paul/LUXdata/'

filename_prefix=filename_evt(1:19)

data_path_evt = [pathbase filesep filename_prefix filesep];
filename_rq = strrep(filename_evt,'evt','rq');
data_path_rq = [pathbase filesep];

 lug_iqs_xml_file = '~/MATLAB/lug_iqs_new4.xml'
iq_xml_path=lug_iqs_xml_file;

%************check dp_settings_new for max num pulses ******* also check
%rubiks cube

data_processing_settings_path = '~/MATLAB/dp_settings_new.xml'
data_processing_xml_path =data_processing_settings_path 
%might want to edit /home/paul/Matlab/dp_settings_new.xml for max num pulses

% Don't use % status = InitializeRQFile_Default(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_settings_path,lug_iqs_xml_file);

%% Here we start with some cpp modules
%in order to get the cpp modules to work, they must first be complied. Go
%to the dp framework folder for the relevant cpp files version and type
%'make' . That should have all the modules in the bin folder to run a
%module, just use its full path a give it the same 6 parameters. if you’re
%inside the CppModules folder, binocessing_settings_path/PhotonTiming path_to_evt_file . I have added the path to .bashrc. these
%can only be done in the linux terminal.

%use instead: 

system(['/home/paul/LUXcode/Stable_Releases/v2.0/DataProcessing/CppModules/bin/InitializeRQFile_Initialize ' filename_evt ' ' data_path_evt ' ' filename_rq ' ' pathbase ' 1 ' data_processing_settings_path ' ' lug_iqs_xml_file ])

%status = PulseCalibration_BaselineZen(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_settings_path,lug_iqs_xml_file);

%n_divisions = make_smaller_sim(filename_evt,data_path_evt,filename_rq,data_path_rq , num_in_each);
n_divisions = split_sim_early(filename_evt,data_path_evt,filename_rq,data_path_rq, num_in_each);

for i = 1:n_divisions
    filename_evt = [filename_evt(1:end-6) num2str(i) '0.evt']
    filename_rq = strrep(filename_evt,'evt','rq');
  
status = PulseCalibration_BaselineZen_split_sim(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_settings_path,lug_iqs_xml_file);

status = PODSummer_LUXSumPOD(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_settings_path,lug_iqs_xml_file);

%{ 
I might need to compile perusePeeks.c before
PulseFinder_TransparentRubiksCube will work. This code should do it. 
cd ~/LUXcode/Stable_Releases/DataProcessingStableReleases/v2.0/DataProcessing
cd MatlabModules/
cd PulseFinder_TransparentRubiksCube/
mex perusePeeks.c 
%}


status = PulseFinder_TransparentRubiksCube(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_settings_path,lug_iqs_xml_file);

system(['/home/paul/LUXcode/Stable_Releases/v2.0/DataProcessing/CppModules/bin/PulseTiming_HeightTiming ' filename_evt ' ' data_path_evt ' ' filename_rq ' ' pathbase ' 5 ' data_processing_settings_path ' ' lug_iqs_xml_file ])

%not used status = PulseTiming_PerusePeeksMatlab(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
%don't use: status = PulseTiming_BasicSet(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path); % make sure the the xml has the snipet added for this module

status = PulseQuantities_MinimumSet_split_sim(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);


status = PulseQuantities_PhotonCounting_split_sim(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
status = PulseClassifier_MultiDimensional(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);

status = S1S2Pairing_Naive(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);

status = Event_Classification(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%for my purposes I might not need Event_Classification

%don't need% PositionReconstructizon_CorrCentroid lux10_20130506T2323_f000502_eb00020.evt /home/paul/LUXdata/lux10_20130506T2323/ lux10_20130506T2323_f000502_eb00020.rq /home/paul/LUXdata 13  /home/paul/Matlab/data_processing_settings.xml /home/paul/Matlab/iqs.xml
% hitmap module% don't need

PositionReconstruction_MercuryI (filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path);
Corrections_PositionCorrection_split_sim(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
Corrections_ApplyCorrections(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)

% 170608 tomaz's talk recommended removal of EnergyReconstruction_Naive(
%EnergyReconstruction_Naive(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
PulseQuantities_TimeSince(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path) %deals with etrains, and looks at the time since the last 'big' event

%PulseQuantities_WaterPmtRQs lux10_20130506T2323_f000502_eb00020.evt /home/paul/LUXdata/lux10_20130506T2323/ lux10_20130506T2323_f000502_eb00020.rq /home/paul/LUXdata 20  /home/paul/Matlab/dp_settings_new.xml /home/paul/Matlab/lug_iqs_new4.xml
%TriggerRQModule(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%%

%status =AdditionalFileFormat_SaveMatFile(filename_evt,data_path_evt,filename_rq,data_path_rq,data_processing_xml_path,iq_xml_path)
%this is to save the .rq file as a .mat file saved /home/paul/LUXdata/matfiles//lux10_20130506T2323_f000502_eb00020.rq.mat
end
for ii = 1:n_divisions
    filename_evt = [filename_evt(1:28) num2str(ii) '0.evt']
    filename_rq = strrep(filename_evt,'evt','rq');
    dp = LUXLoadRQ1s_framework(filename_rq, data_path_rq);

    rq_list = fieldnames(dp);
    if ii == 1
        dp_total = dp;
    else
        for jj=1:length(rq_list)
            current_rq = rq_list{jj};
            if strcmp(current_rq, 'admin') | strcmp(current_rq, 'source_filename') %...
                  %  | strcmp(current_rq, 'luxstamp_samples') | strcmp(current_rq, 'time_since_livetime_start_samples') ...
                   % | strcmp(current_rq, 'time_until_livetime_end_samples')
                %do nothing. hopefully this will be ok. 
            else
                dp_total.(current_rq) = cat(ndims(dp.(current_rq)), dp_total.(current_rq), dp.(current_rq));
            end
        end
    end
end

filename_rq = [filename_rq(1:19) '_merged.rq'];
LUXWriteRQMFile_framework(filename_rq, [data_path_rq '/matfiles/'], dp_total, dp_total.admin);
fprintf('saving mat file %s.mat\n',filename_rq);
       
    
%exit
end
