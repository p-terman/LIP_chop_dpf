function PlotTopArray(White, TDs)
% This function will plot a 3D view of the LUX detector with numbers, on the current
% figure.
%
% LUXPlotTopPMTs
%
% No inputs or outputs.
%
% Versioning
%   20130104 CHF - Created using LUXPlotTopPMTs
%   20130128 CFPS - Simplified to show only the top array
%% Initialize

if nargin < 2
    TDs = 0;
end

if nargin < 1
    White = 0;
end

load pmt_pos_map

pmt_width = 5.7; % cm
tt = 0:pi/64:2*pi;
xx = pmt_width./2 .* cos(tt);
yy = pmt_width./2 .* sin(tt);

%% Plot

for ii_pmt=[1:60 121]
    hold on
    ii_pmt_name = ii_pmt;
    if White == 0,
        plot(pmt_pos_cm(ii_pmt,1)+xx, pmt_pos_cm(ii_pmt,2)+yy,'color',[1 1 1]*0);
    else
        plot(pmt_pos_cm(ii_pmt,1)+xx, pmt_pos_cm(ii_pmt,2)+yy,'color',[1 1 1]*White);
    end    
        
end

ww = 18.632 * 2.54; % ptfe panel face-face width, in->cm
ll = ww/2;
alpha_d = 30; % degrees
beta_d = alpha_d/2;
xy = [-ll*tand(beta_d) ll*tand(beta_d);ll ll];
xy = [cosd(15) -sind(15);sind(15) cosd(15)]*xy;
for ii = 1:12
  hh = line(xy(1,:),xy(2,:),[0 0]);
  if White==0,
      set(hh,'color',[1 1 1]*0);
  else
      set(hh,'color',[1 1 1]*White);
  end
  Ralpha = [cosd(alpha_d) -sind(alpha_d);sind(alpha_d) cosd(alpha_d)];
  xy_new(:,1) = Ralpha*xy(:,1);
  xy_new(:,2) = Ralpha*xy(:,2);
  xy = xy_new;
end;
if ~TDs
    axis([-25 25 -25 25])
    axis square
end

%% Pretty it up 

% axis equal;
