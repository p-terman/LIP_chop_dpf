function [num_pC_main APR tus ymean tot_sum range APratio] = LUXGetAP(pathbase,filename_prefix,num_ch,bias,amp)
%{

LUXGetAP(pathbase,filename_prefix,ch,bias,amp)

Converts an event to a trace and calculates the afterpulsing and size of main pulse.

Inputs:

           pathbase - path to dataset
    filename_prefix - name of dataset
                 ch - channel
               bias - bias voltage in V
                amp - amplification factor


Outputs:

        num_pC_main - Size of main pulse in picocoulombs
                APR - afterpulsing ratio, afterpulses/mainpulse
                tus - time data in us
              ymean - ydata in mV
            tot_sum - struct. containing the sum of the main pulse, and the sum of the afterpulses

20090525 CHF
20090928 JRV
20091029 CHF - Removed event list, added amp as an input
20100114 JRV - Now one call gets AP info. for all channels

%}

%% Load settings

xml_settings = LUXLoadSettings(filename_prefix,pathbase);
        
%% Load events
event_struct = LUXPMTEventLoader(xml_settings);
event_traces = LUXEvent2Trace_old(xml_settings,[],event_struct);
pulses_all = squeeze(-event_traces.yy_mV)/amp;

Wn = 0.5;
[b a] = butter(1,Wn); % apply butterworth filter to trace

for ii = 1:num_ch
    pulses{ii} = filter(b,a,squeeze(pulses_all(:,ii,:)));
end

% Use no sphe suppresion at this point
% pulses(pulses<20) = 0;

%% Get mean of pulses
for ii = 1:num_ch
    ymean{ii} = mean(pulses{ii},2);
end

for ii = 1:num_ch
    [tpk{ii} pk{ii}] = max(ymean{ii});
end

% Hard coded right now, must calculate ranges from bias
mainbegins= -0.05; %us
afterbegins = 0.10; % us, when main pulse ends and AP begins
afterends = 5;

for ii = 1:num_ch
    dt = 10/1.05; %ns
    t = (0:(length(ymean{ii})-1)) * dt;
    tus{ii} = (t - t(pk{ii}))/1000 + 0.02;


    maincut{ii} = inrange(tus{ii},mainbegins,afterbegins);
    pulse_main{ii} = ymean{ii}(maincut{ii});
    aftercut{ii} = inrange(tus{ii},afterbegins,afterends);
    afterpulses{ii} = ymean{ii}(aftercut{ii});

    tot_sum{ii}.mainsum = sum(pulse_main{ii});
    tot_sum{ii}.aftersum = sum(afterpulses{ii});

    APR{ii} = tot_sum{ii}.aftersum / tot_sum{ii}.mainsum;
    num_pC_main{ii} = ceil(tot_sum{ii}.mainsum*dt/25); %in pC
    num_pC_ap{ii} = ceil(tot_sum{ii}.aftersum*dt/25); %in pC
end
%% get individual element afterpulsing

q = -1.609e-19; % electron charge
A = 6.02e23; % Avogadro's number
d_ions = 6.10e-2; % approximate distance where ions are formed inside PMT 

ion_list  = ['H+   ';'He+  ';'N+,O+';'Ar+  ';'Xe+  '];
mass_list = [1; 4; 15; 40; 131] ./ (1e3*A); % kg

for ii = 1:num_ch
    t_list{ii} = (d_ions).*sqrt(2.*(mass_list)/(q*bias{ii})) ./ 1e-6; % ion transit time in us

    dt = 10/1.05/1000;

    % 20091030 CHF - added hard-coded AP ranges, cleaned-up code
    range{ii}(1).cut_range = inrange(tus{ii},[0.20 0.35]); %H, empirical CHF
    range{ii}(2).cut_range = inrange(tus{ii},[0.39 0.65]); %He, empirical CHF
    range{ii}(3).cut_range = inrange(tus{ii},[0.77 1.20]); %N,O, empirical CHF
    range{ii}(4).cut_range = inrange(tus{ii},[1.30 1.75]); %Ar, empirical CHF
    range{ii}(5).cut_range = inrange(tus{ii},[2.50 3.00]); %He, empirical CHF

    for jj = 1:length(t_list{ii})
        range{ii}(jj).time_range = tus{ii}(range{ii}(jj).cut_range);
        afterpulses{ii}(jj) =  sum(ymean{ii}(range{ii}(jj).cut_range));
    end

    APratio{ii}.individual = afterpulses{ii}./tot_sum{ii}.mainsum.*100;
    APratio{ii}.main_total = tot_sum{ii}.aftersum/tot_sum{ii}.mainsum.*100;
end

