function APGlobalStatus(serials)
%
% Use serials = getPMTserials to generate full R8778 list
%
% This function will bar plot latest measured AP levels for all LUX PMTs
%
% 20100218 CHF - Created
%
%% Initialize

xx = 1:numel(serials);

APmat = zeros(numel(serials),3);

%% Get latest AP for each PMT

for ii = xx
    
    try
        [AP dates] = getAPfromPMTdb(serials{ii},1); % Use this function to get it from LUG
    catch
        AP.H = 0;
        AP.He = 0;
        AP.NO = 0;
        AP.Ar = 0;
        AP.Xe = 0;
        % If it can't be found in LUG, zero it
    end
    
    APmat(ii,1:3) = [AP.H AP.He AP.NO]; % Build matrix with values for bar stack plotting

end

APmat(APmat < 0) = 0; %Zero all negative AP values

%% Plot

edge = xx(end) + 0.5; % right limit for x
col = [0.6 0.7 0.9]; % Fill color

% Separate advance order PMTs from LUX order
plot(12*[1 1] + 0.5,[0 50],'k--','handlevisibility','off'); hold on

% Make a fill rectangle for AP 'OK' threshold
a = fill([0 0 edge edge],[0 5 5 0],col,'LineStyle','none','handlevisibility','off'); hold on
alpha(a,0.2); % Make fill transparent

% Make bar plot of all data
bar(xx,APmat,0.8,'stacked');
colormap([0 0 0; 1 0 0 ;0 0 1]); %k r b

% Label it and pretty it up
ylabel('AP (%)')
title('The Current State of Affairs')
set(gca,'XTickLabel',serials,'XTick',xx);
box('on');
rotateticklabel(gca,90,8);
legend('H^+','He^+','N^+,O^+','location','NE');

axis([0 edge 0 40])
    