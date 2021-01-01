function APHist(PMTserial)

%% Query AP history about this PMT

bias = -1300;
limit = 10;

try
    [AP dates] = getAPfromPMTdb(PMTserial,limit,bias);
    different = 0;
catch
    [AP dates] = getAPfromPMTdb(PMTserial,limit);
    different = 1;
end

N = size(AP,2);

%% Combine into plot data
t = zeros(1,N);
APmat = zeros(N,3);

for ii = 1:N
    APmat(ii,1:3) = [AP(ii).H AP(ii).He AP(ii).NO];
    t(ii) = dates{ii}.datenum;
    tlabel{ii} = dates{ii}.label(1:10);
end

%% Plot history

% figure(98); clf

plot(t,APmat(:,1),'ko-'); hold on
plot(t,APmat(:,2),'ro-')
plot(t,APmat(:,3),'bo-')

ylabel('AP (%)')

if different
    tit(sprintf('Afterpulsing History for PMT %s, low bias',PMT));
else
    tit(sprintf('Afterpulsing History for PMT %s',PMTserial));
end

legend('H^+','He^+','N^+,O^+','location','NW');

if max(max(APmat)) < 5
    ylim([0 5])
end

[tu Iu] = unique(t);
tlabelu = {tlabel{Iu}};
[ti I] = sort(tu);
tlabeli = tlabelu{I};

set(gca,'XTickLabel',tlabelu,'XTick',ti);
rotateticklabel(gca,30);
