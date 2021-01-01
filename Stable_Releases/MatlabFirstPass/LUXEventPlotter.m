function LUXEventPlotter(filename_prefix, data_path, event_list, event_struct, options)
% LUXEventPlotter(filename_prefix, data_path, event_list, event_struct,
% options)
%
% Inputs:
% filename_prefix is name of dataset to load
% data_path is path to data
% event_list is the list of event numbers that you wish to plot.
% event_struct is the structure of events that has been loaded from
% LUXEventLoader.m - optional, will load it if empty
%
% options:
% flag_axis_auto 1 - axis auto after plotting
% flag_cal_pulses 1 - plot in phe/sample (call LUXCalibratePulses)
% flag_sum_pulses 1 - sum and plot in top plot (call LUXSumPOD)
% flag_plot_ns 1 - change scale to nanoseconds from samples
% flag_plot_rq1s 1 - show choice RQ1 variables on the plot (and tell value
% in output) NOTE: this needs rq1s_path - a string with the path to load
% the RQ1s from.
%
% 2008-07-02 jjc
% 2008-07-15 LdV - renamed it from plot_event to LUXEventPlotter
% 2008-07-15 LdV - it plots with the baseline subtracted
% 2009-05-31 DCM - if event_list isn't specified or is [], LUXEventPlotter
%  loops through all events in event_struct
% 2010-02-02 JJC - fixing so it doesn't suck so much. calls LUXEventLoader
% if you want, has some plot options too.
%2010-03-24 JJC - adding the option to plot RQ1s on the sum.
%2010-05-27 JJC - v5 

%figure(100);


if ~exist('data_path','var') || isempty(data_path)
	data_path=[];
end

%% Create options variables
if nargin > 4
    field_names = fieldnames(options);
    for ii=1:length(field_names)
        eval([field_names{ii} ' = options.' field_names{ii} ';']);
    end
end

if ~exist('flag_axis_auto','var')
    flag_axis_auto = false;
end
if ~exist('flag_cal_pulses','var')
    flag_cal_pulses = false;
end
if ~exist('flag_sum_pulses','var')
    flag_sum_pulses = false;
end
if ~exist('flag_plot_ns','var')
    flag_plot_ns = false;
end
if ~exist('flag_plot_rq1s','var')
    flag_plot_rq1s = false;
elseif ~exist('rq1s_path', 'var')
    fprintf('You asked to plot RQ1s, but did not provide the RQ1s path, so I wont do it.\n');
    flag_plot_rq1s = false;
elseif flag_sum_pulses == false
    fprintf('You need to sum the pulses first. Dont worry, I will do that for you...\n');
    flag_sum_pulses = true;
else
    RQ1s_path = rq1s_path;
end
    

%if ~exist('flag_axis_auto','var') || isempty(flag_axis_auto)
%	flag_axis_auto=0;
%end
%if ~exist('flag_cal_pulses','var') || isempty(flag_cal_pulses)
%	flag_cal_pulses=0;
%    flag_sum_pulses=0;
%end
%if ~exist('flag_sum_pulses','var') || isempty(flag_sum_pulses)
%    flag_sum_pulses=0;
%end
%keyboard

if ~exist('filename_prefix','var')
    filename_prefix=[];
end

if ~exist('event_struct','var')
    event_struct = [];
end

if isempty(filename_prefix) && isempty(event_struct)
    if isempty(filename_prefix)
    fprintf('Looking for the latest dataset in %s ...\n', data_path);
    latest_filename_list = {};
    latest_filename_n = [];
    
    % Making assumption here -- all files / datasets are to start with 'lux' prefix
    % Note -- inclusive with LUX 0.1 files, so this is OK
    filename_list = dir([data_path,filesep,'lux*']);
    if ~isempty(filename_list)
        % Check whether the folder listed in filename_list fits the right format
        ii_file = 1;
        while ii_file < length(filename_list)
            if filename_list(ii_file).isdir
                % test to see if bit after underscore is a well-formed date
                % string
                [prefix, rem] = strtok(filename_list(ii_file).name, '_');
                datestring = strtok(rem, '_');
                try
                    datenum = datestr(datestring, 'yyyymmddTHHMM');
                    ii_file = ii_file+1;
                catch
                    filename_list(ii_file) = [];
                end
            else
                filename_list(ii_file) = [];
            end
        end
    end
    
    if isempty(filename_list) % no results found or left after name checking
        fprintf('There are no datasets in %s\n', data_path);
        return;
    else
        dates = [filename_list.datenum];
        [tmp, ii_latest] = max(dates);
        filename_prefix = filename_list(ii_latest).name;
        fprintf('Loading settings from dataset %s\n', filename_prefix);
    end
end

end

if isempty(event_struct)
   %s = LUXSuperLoader_v5(filename_prefix, data_path);
   %if ~exist('event_list','var')
   % event_list = LUXNumberEvents(s.daq_settings);
   %end
   
   [event_struct settings] = LUXEventLoader_v5(filename_prefix,data_path,event_list);
elseif isempty(event_list)
    event_list = length(event_struct);
end

if flag_cal_pulses && ~isfield(event_struct(end).ch(1).pulse(1),'pulse_data_phe')
    fprintf('calibrating pulses...\n');
    options.query_lug=1;
    %s=LUXSuperLoader_v5(filename_prefix, data_path,options);
    s.lug_settings = LUX01LoadState(filename_prefix);
    event_struct = LUXCalibratePulses(event_struct, s);
end

if flag_sum_pulses && ~isfield(event_struct(end),'chsum')
    fprintf('summing pulses...\n');
    event_struct = LUXSumPOD(event_struct);
end

if flag_plot_rq1s % load the rq1s for this data to plot them
    dp = LUXLoadRQ1s(filename_prefix, RQ1s_path, event_list);
end

scrsz = get(0,'ScreenSize');
h=figure(100);
set(100,'Position',[1 scrsz(4)/2 scrsz(3)/2 scrsz(4)/2]);

xlab = 'samples';
if flag_plot_ns
    xlab = 'ns';
end
ylab = 'mV';
if flag_cal_pulses
    ylab = 'phe/sample';
    if flag_plot_ns
        ylab = 'phe/10ns';
    end
end

for ii_evt=event_list
    clf
    %title(['Event ' num2str(ii_evt)], 'fontsize',13,'fontweight','b');
    if flag_sum_pulses
        subplot(2,1,2);
        %title('individual channels');
    end
    legstr = {}; % Initialize legend matrix
    %nextcol({'b', 'r', 'g', 'k', 'c', 'm'}); % Plot colors
    col{1} = [0 0 1];
    col{2} = [1 0 0];
    col{3} = [0 .6 0];
    col{4} = [0 0 0];
    if event_struct(ii_evt).empty ==1
        disp('There are no pulses in this event!');
    end
    if event_struct(ii_evt).empty==0

        %xlabel('Samples');
        %ylabel('ADC bins');
       % title(['Event ' num2str(ii_evt)], 'fontsize',13,'fontweight','b');


        nb_chs = length(event_struct(ii_evt).ch);

        for ch=1:nb_chs
            plot(0,0,'color',col{ch}, 'LineStyle','.')
            hold on
            %legstr{end+1} = sprintf('ch %d',ch);
        end

        for ch=1:nb_chs
            if isfield(event_struct(ii_evt).ch(1), 'pulse')
                for ps=1:length(event_struct(ii_evt).ch(ch).pulse)
                    if flag_cal_pulses==1
                    	tp = double(1:1:length(event_struct(ii_evt).ch(ch).pulse(ps).pulse_data_phe));
                    	xx = double(tp + double(event_struct(ii_evt).ch(ch).pulse(ps).start));
                        if flag_plot_ns
                            xx = xx.*9.52;
                        end
                    	yy = -event_struct(ii_evt).ch(ch).pulse(ps).pulse_data_phe;
                    else
                    	tp = double(1:1:(event_struct(ii_evt).ch(ch).pulse(ps).length-1));
                    	xx = double(tp + double(event_struct(ii_evt).ch(ch).pulse(ps).start));
                    	if flag_plot_ns
                            xx = xx.*9.52;
                        end
                        yy_data = double(event_struct(ii_evt).ch(ch).pulse(ps).pulse_data_mV(2:end));
                    	yy_baseline = event_struct(ii_evt).ch(ch).pulse(ps).baseline_mV;
                    	yy = -(yy_data - yy_baseline);
                    end
                    plot(xx, yy, 'color', col{ch}, 'LineStyle','-');
                        
                        
                    hold on

                end
            end
          
            
        end
        xlabel(xlab);
        ylabel(ylab);
        %legend(legstr,'Location','SouthEast');
        if flag_sum_pulses
             if flag_plot_rq1s
                clear t0 t2 area;
                if length(event_list)>1
                t0 = dp{1}.t0(:,ii_evt) - settings.evt_settings.pretrigger;
                t2 = dp{1}.t2(:,ii_evt) - settings.evt_settings.pretrigger;
                t1 = dp{1}.t1(:,ii_evt) - settings.evt_settings.pretrigger;
                area = dp{1}.peak_area_per_ch(:,:,ii_evt);
                whichpeak = dp{1}.whichpeak(:,ii_evt);
                whichpulse = dp{1}.whichpulse(:,ii_evt);
                else
                    t0 = dp{1}.t0 - settings.evt_settings.pretrigger;
                    t2 = dp{1}.t2 - settings.evt_settings.pretrigger;
                    t1 = dp{1}.t1 - settings.evt_settings.pretrigger;
                    area = sum(dp{1}.peak_area_per_ch,2);
                    whichpeak = dp{1}.whichpeak;
                    whichpulse = dp{1}.whichpulse;
                end
            end
            subplot(2,1,2);
            title('Individual Channels','fontsize',13,'fontweight','b');
            subplot(2,1,1);
            if isfield(event_struct(ii_evt).chsum,'pulse')
                for ps=1:length(event_struct(ii_evt).chsum.pulse)
                    xx = event_struct(ii_evt).chsum.pulse(ps).tt;
                    if flag_plot_ns
                        xx = xx.*9.52;
                    end
                    yy = -event_struct(ii_evt).chsum.pulse(ps).data_phe;
                    plot(xx, yy, 'Color', [0.3 0.3 0.7]);
                    hold on;
                    if flag_plot_rq1s
                        clear Cl Cr IAl IBl IAr IBr;
                        [Cl,IAl,IBl] = intersect(t0,xx);
                        for jj=1:length(Cl)
                           plot(t0(IAl(jj)),yy(IBl(jj)),'o','MarkerSize',(2*(length(Cl)-jj+1)),'MarkerEdgeColor','m'); 
                           
                           
                        end
                        [Ca,IAa,IBa] = intersect(t1,xx);
                        subplot(2,1,1);
                        for jj=1:length(Ca)
                           h=text(t1(IAa(jj)),max((yy(IBa(jj))*1.05),(yy(IBa(jj))+150)),sprintf('%4.0f',sum(area(IAa(jj)),2)));
                           set(h,'FontSize',8);
                           if whichpulse(IAa(jj))==1
                              marker='mo'; 
                           end
                           if whichpulse(IAa(jj))==2
                               marker='co';
                           end
                           plot(t1(IAa(jj)),yy(IBa(jj)),marker);
                           
                        end
                        
                        
                        for jj=1:length(Ca)
                           subplot(2,1,2);
                           h=text(t1(IAa(jj)),300+(yy(IBa(jj))*1.05),sprintf('%4.0f %4.0f %4.0f %4.0f',area(IAa(jj),1)));
                           set(h,'FontSize',8,'Color',col{1});
                           h=text(t1(IAa(jj)),200+(yy(IBa(jj))*1.05),sprintf('%4.0f %4.0f %4.0f %4.0f',area(IAa(jj),2)));
                           set(h,'FontSize',8,'Color',col{2});
                           h=text(t1(IAa(jj)),100+(yy(IBa(jj))*1.05),sprintf('%4.0f %4.0f %4.0f %4.0f',area(IAa(jj),3)));
                           set(h,'FontSize',8,'Color',col{3});
                           h=text(t1(IAa(jj)),(yy(IBa(jj))*1.05),sprintf('%4.0f %4.0f %4.0f %4.0f',area(IAa(jj),4)));
                           set(h,'FontSize',8,'Color',col{4});
                           ax=axis; axis([ax(1) ax(2) ax(3) max(ax(4), 200+yy(IBa(jj)))]);
                        end
                        subplot(2,1,1);
                        [Cr,IAr,IBr] = intersect(t2,xx);
                        for jj=1:length(Cr)
                           plot(t2(IAr(jj)),yy(IBr(jj)),'o','MarkerSize',(2*(length(Cr)-jj+1)),'MarkerEdgeColor','m');
                           
                        end
                        
                        %plot(t10l(IA),yy(IB),'x','MarkerSize',6,'MarkerEdgeColor','m');
                        %hold on;
                        %[C,IA,IB] = intersect(t10r,xx);
                        %plot(t10r(IA),yy(IB),'x','MarkerSize',6,'MarkerEdgeColor','c');
               % keyboard
                    end
                end
            end
           
            title(['Event ' num2str(ii_evt) ' - Sum'], 'fontsize',13,'fontweight','b');
        else
            title(['Event ' num2str(ii_evt)], 'fontsize',13,'fontweight','b');
        end
        
        x2 = settings.evt_settings.posttrigger;
        x1 = - settings.evt_settings.pretrigger;
        ax=axis; axis([x1 x2 ax(3:4)]);
        if flag_axis_auto==1
        	axis auto;
            ax = axis;
            axis([ax(1) ax(2) ax(3) ax(4)*1.1]);
        end
        xlabel(xlab);
        ylabel(ylab);
        if flag_plot_rq1s
           
        end
        
        
        %legend('sum','Location','SouthEast');

        disp(sprintf(['Event #' num2str(ii_evt) '. Click in figure window to advance to the next event. Ctrl-C to escape.']))
        waitforbuttonpress;
    end
    hold off
    clf
end