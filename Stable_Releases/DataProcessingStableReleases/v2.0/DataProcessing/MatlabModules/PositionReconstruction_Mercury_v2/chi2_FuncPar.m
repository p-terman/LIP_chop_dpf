function chi2 = chi2_FuncPar(p_min, data, data_not_eq, pulse_area_eq, cutEv, lrf_iq)
%%
%% That's probably should defined as global.
%%

[pmt_pos sextant_arrangement] = LUXPMTArray; 
for i=1:122,
    PMT_r(i,1:3)=pmt_pos(find(sextant_arrangement==i),:);
    the_sextant(i,:) = Find_Sextant(i, pmt_pos, sextant_arrangement);
end
[PMT_Radius PMTG] = LUXLRF_GroupDefinition;

topchs = [1:60 121];

%%    
%% Definitions of the variables that are used
%%    

XX = p_min(1); % (XX YY are the estimated positions of the event in question)
YY = p_min(2);
rho = sqrt((XX-PMT_r(topchs,1)).^2+(YY-PMT_r(topchs,2)).^2); % Variable 1 - Distance PMT-Event
radius = sqrt(XX.^2+YY.^2); % Variable 2 - the radius of the event
phi = atan(XX./YY); % Variable 3 - the phi
phi_cut =  inrange(abs(phi), pi*1./6, pi/2);
phi_sextant = abs(abs(phi)-phi_cut*pi/3);
dis_event_wall = 24.*cos(pi/12)./cos(abs(abs(abs(phi)-phi_cut*pi/3)-pi/12))-radius; 

%%
%%

rho_Evcut = rho(cutEv);
NumPMTs = sum(cutEv);

%%
%% The for along the PMT groups according the radius
%%

the_Direct_Component = lrf_iq.rad.CallFunction_In(rho(cutEv), lrf_iq.rad.Parameters);

for PMTGroup = 1:9,
    Tau_Decay = lrf_iq.ref.Tau_Decay_Function_In(dis_event_wall, lrf_iq.ref.Tau_Decay(PMTGroup,1), lrf_iq.ref.Tau_Decay(PMTGroup,2), lrf_iq.ref.Tau_Decay(PMTGroup,3));
    Amplitude = lrf_iq.ref.Amplitude_Function_In(dis_event_wall, lrf_iq.ref.Exp_Amplitude(PMTGroup,1), lrf_iq.ref.Exp_Amplitude(PMTGroup,2), lrf_iq.ref.Exp_Amplitude(PMTGroup,3));     
    cutPMTs = PMT_Radius == PMTG(PMTGroup);   
    Reflected_Component(cutPMTs(topchs)) = Amplitude*exp(-rho(cutPMTs(topchs)).*Tau_Decay);
    %% This accounts for the symmetry of the array of 22 cm.
    if PMTGroup == 8 & inrange(phi_sextant, 0.2426, pi/6),
        PMTGroup = 18;
        Tau_Decay = lrf_iq.ref.Tau_Decay_Function_In(dis_event_wall, lrf_iq.ref.Tau_Decay(PMTGroup,1),lrf_iq.ref.Tau_Decay(PMTGroup,2), lrf_iq.ref.Tau_Decay(PMTGroup,3));
        Amplitude = lrf_iq.ref.Amplitude_Function_In(dis_event_wall, lrf_iq.ref.Exp_Amplitude(PMTGroup,1), lrf_iq.ref.Exp_Amplitude(PMTGroup,2), lrf_iq.ref.Exp_Amplitude(PMTGroup,3));
        Reflected_Component(cutPMTs(topchs)) = Amplitude*exp(-rho(cutPMTs(topchs)).*Tau_Decay);
    end
end
the_Reflected_Component = Reflected_Component(cutEv)';

%% The version of chi-square used right now is the reduced chi-square. That's the reason why we are dividing it by the number of PMTs

DOF = max(NumPMTs-2, 1); % The degrees of freedom
chi2 = sum((pulse_area_eq*(the_Reflected_Component+the_Direct_Component)-data(cutEv)').^2./(DOF*(abs(data_not_eq(cutEv)')+1)));
