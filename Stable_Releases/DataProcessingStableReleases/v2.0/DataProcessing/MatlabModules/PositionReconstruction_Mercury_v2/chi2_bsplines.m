function chi2 = chi2_BSplines(p_min, data, data_not_eq,pulse_area_eq, cutEv, lrfmat)

%%    
%% Definitions of the variables that are used
%%    


XX = p_min(1); % p_min is in polar coordinates
YY = p_min(2);
rho = sqrt((XX-lrfmat.PMT_r(lrfmat.topchs,1)).^2+(YY-lrfmat.PMT_r(lrfmat.topchs,2)).^2); % Variable 1 - Distance PMT-Event
phi = atan(abs(XX)./abs(YY)); % Variable 3 - the phi
phi_cut =  inrange(abs(phi), pi*1./6, pi/2);
dis_event_wall = max(24.5*cos(pi/12)./cos(abs(abs(abs(phi)-phi_cut*pi/3)-pi/12))-sqrt(XX.^2+YY.^2), -1.0); 

%% Get the variables that are necessary for the spline calculation
phi_B = atan2(XX, YY);
Delta_Phi_pr  = abs(lrfmat.PMT_phi - repmat(atan2(YY, XX), [122 1]));
Delta_Phi_Cut = Delta_Phi_pr > pi;
Delta_Phi_A = (2*pi-Delta_Phi_pr).*Delta_Phi_Cut + Delta_Phi_pr.*~Delta_Phi_Cut;        
Delta_Phi = Delta_Phi_A(lrfmat.topchs);
%%
%% Definition of the cuts in rho (distance PMT-event) and the photomultiplier amplitude
%%

rho_Evcut = rho(cutEv);
NumPMTs = sum(cutEv);
the_Direct_Component = lrfmat.Rho(round(10*(60-min(rho_Evcut,59)).^2))';

listpmts = 1:61;
mpta = 1;
for mpt = listpmts(cutEv)
    UpVal = lrfmat.spval(mpt, min(max(ceil((dis_event_wall+1/16)*8), 1), 160), min(round(Delta_Phi(mpt).*3600/pi+1), 1800));
    MedVal = lrfmat.spval(mpt, min(max(ceil((dis_event_wall)*8), 1), 160), min(round(Delta_Phi(mpt).*3600/pi+1), 1800));
    DownVal = lrfmat.spval(mpt, min(max(ceil((dis_event_wall-1/16)*8), 1), 160), min(round(Delta_Phi(mpt).*3600/pi+1), 1800));
    the_Reflected_Component(mpta) = mean([UpVal MedVal DownVal]);
    mpta = mpta+1;
end 

DOF = max(NumPMTs-2, 1); % The degrees of freedom

if dis_event_wall > 0
    SystError = 1+0.001*(XX.^2+YY.^2);
else
    SystError = 0.0001;
end    
   
chi2 = sum((pulse_area_eq*(the_Direct_Component+the_Reflected_Component')-data(cutEv)').^2./(DOF*SystError*(abs(data_not_eq(cutEv)'))));
