function [dp rec_set] = MercurySelectEventsForReconstruction(dp, rec_set)

% This function selects the events for the reconstruction. It can
% select everything or only some kind of pulses
%
% [dp rec_set] = MercurySelectEventsForReconstruction(dp, rec_set)
%
% Inputs:
%         dp     - the data array
%         rec_set - the reconstruction settings
% Outputs:
%         dp     - the data array with the field dp.Reconstructed
%         rec_set - the reconstruction settings
% Example usage:
% 
%
% 20130212 CFPS - Created
% 20130319 CFPS - Modified to work without the field dp.Reconstructed
% 20130325 CHF  - Minor fix near line 52, where a stray
%                 pulse_classification was giving errors
% 20130326 CFPS - A new option for the selection of the events: select_all
% 20130701 CFPS - Some code cleaning. The option to qualify the 
% 20140507 CFPS - File Cleaned
%%
%% Input Check
%%


