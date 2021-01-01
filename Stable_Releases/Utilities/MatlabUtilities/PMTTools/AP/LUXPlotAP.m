function LUXPlotAP(tus,ymean,PMTserial,num_pC_main,APratio,filename,bias,range,Ch)
%
% Plots the afterpulsing analysis in PMTHealthMonitor.m
%
% Inputs: 
%                tus - time data in us
%              ymean - ydata in mV
%          PMTserial - the pmt serial number
%        num_pC_main - size of main pulse in picoCoulombs
%            APratio - area afterpulses / area main pulse for each ion population and total
%           filename - dataset name
%               bias - pmt bias voltage in Volts
%              range - gives range of times for afterpulses of each element
%                 ch - channel number
% 
% 20091028 JRV - created
% 20091030 CHF - changed APR input to APratio, which contains all ion populations as well 


colorset = 'brkm';
plot(tus{Ch},zeros(length(tus{Ch}),1),'k','HandleVisibility','off'); hold on % plot zero baseline
APR = APratio.main_total;

name = filename;
[maxvalue max_index] = max(ymean);
plot(tus{Ch}-tus{Ch}(max_index),ymean,colorset(1)); hold on
leg = sprintf('Dataset: %s\n%2.1f%% APR (%d pC signal) \nPMT Bias: %d V \nChannel: %d ',name(7:end),APR,num_pC_main,bias,Ch);


axis([-0.2 4 -0.5 maxvalue/20]);
ax2 = axis;

%% Place afterpulsing ion locators
q = -1.609e-19; % electron charge
A = 6.02e23; % Avogadro's number
d_ions = 6.10e-2; % approximate distance where ions are formed inside PMT
%bias = 1500; %V

ion_list  = ['H+   ';'He+  ';'N+,O+';'Ar+  ';'Xe+  '];
mass_list = [1; 4; 15; 40; 131] ./ (1e3*A); % kg

t_list = (d_ions).*sqrt(2.*(mass_list)/(q*bias)) ./ 1e-6; % ion transit time in us

for ii = 1:numel(mass_list)
    plot(t_list(ii)*[1 1],ax2(3:4),'m--');
    text(t_list(ii)+0.02,ax2(3)*0.4,ion_list(ii,:),'Color','m');
    text(t_list(ii)+0.02,ax2(3)*0.6,sprintf('%3.1f%%',APratio.individual(ii)),'Color','k');
end

max_plot = max(ymean);

for kk = 1:length(range)
    g(kk) = patch([(range(kk).time_range(kk)) (range(kk).time_range(length(range(kk).time_range))) (range(kk).time_range(length(range(kk).time_range))) (range(kk).time_range(kk))],[-5 -5 max_plot max_plot],[0 0 0 0],'y');
    set(g(kk),'linestyle','none')
end
alpha(g,0.4);

%%

legend(leg,'Location','NE');
title(sprintf('Afterpulsing for %s',PMTserial));
xlabel('Time (us)');
ylabel('mV');