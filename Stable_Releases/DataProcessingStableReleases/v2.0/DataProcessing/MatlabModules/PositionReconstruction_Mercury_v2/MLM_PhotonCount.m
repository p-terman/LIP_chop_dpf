function [chi2 chiE Dir Ref] = MLM_PhotonCount(p_min, data, spike_num, energy_est, cutEv, lrfmat)

phi = atan(abs(p_min(1))./abs(p_min(2)));
if isnan(phi)
    phi=0;
end
WALL = lrfmat.disphi( round(phi*10000/pi+1));
dis_wall = (WALL - sqrt(p_min(1).^2+p_min(2).^2))*24.5/WALL;
dis_event_wall = max(dis_wall, -1.0); 
rho = (sqrt((p_min(1)-lrfmat.PMT_r(lrfmat.topchs(cutEv),1)).^2 + (p_min(2)-lrfmat.PMT_r(lrfmat.topchs(cutEv),2)).^2)); % Variable 1 - Distance PMT-Event

para_dis_event_wall = max(min(round(100*(24.5-dis_event_wall)^2), 65000), 1);
para_rho_wall = max(min(round(100*(p_min(1).^2+p_min(2).^2)), 65000), 1);
the_Direct_Component = diag(lrfmat.Rho(cutEv, round(10*(60-min(rho,59)).^2)));

raiooudist = lrfmat.variabledependence(cutEv) == 1;
the_Reflected_Component = raiooudist'.*0;
the_Reflected_Component(~raiooudist) = (lrfmat.Amp(cutEv & ~lrfmat.variabledependence',para_dis_event_wall).*exp(-rho( ~raiooudist').*lrfmat.Tau(cutEv & ~lrfmat.variabledependence',para_dis_event_wall)));
the_Reflected_Component(raiooudist) = (lrfmat.Amp(cutEv & lrfmat.variabledependence',para_rho_wall).*exp(-rho(raiooudist').*lrfmat.Tau(cutEv & lrfmat.variabledependence',para_rho_wall)));
the_Reflected_Component = the_Reflected_Component + lrfmat.Amp2(lrfmat.topchs(cutEv)).* exp(-rho.*lrfmat.Gam2(lrfmat.topchs(cutEv)));

if numel(p_min)>2
    TotalEnergy_Spike = p_min(3);
    TotalEnergy_Areas = p_min(3);
else
    rho = [sqrt((p_min(1)-lrfmat.PMT_r(5,1)).^2 + (p_min(2)-lrfmat.PMT_r(5,2)).^2) sqrt((p_min(1)-lrfmat.PMT_r(32,1)).^2 + (p_min(2)-lrfmat.PMT_r(32,2)).^2)]';
    the_Direct_Component_NW = [lrfmat.Rho(5, round(10*(60-min(rho(2),59)).^2))' lrfmat.Rho(32, round(10*(60-min(rho(2),59)).^2))']';
    the_Reflected_Component_NW = [lrfmat.Amp(5,para_dis_event_wall).*exp(-rho(1).*lrfmat.Tau(5,para_dis_event_wall))+ lrfmat.Amp2(5).* exp(-rho(1).*lrfmat.Gam2(5))  lrfmat.Amp(32,para_dis_event_wall).*exp(-rho(2).*lrfmat.Tau(32,para_dis_event_wall))+ lrfmat.Amp2(32).* exp(-rho(2).*lrfmat.Gam2(32))]';
    TotalEnergy_Spike = energy_est(1)/(1-sum(the_Direct_Component_NW +  the_Reflected_Component_NW));
    TotalEnergy_Areas = energy_est(2)/(1-sum(the_Direct_Component_NW +  the_Reflected_Component_NW));
end

%% The Poisson Distribution of the observed spikes

Corrections = lrfmat.lrf_iq.QE.QE_Values(cutEv);

mu_val = max((TotalEnergy_Spike*(the_Direct_Component(lrfmat.Cut_Pile_Up)+the_Reflected_Component(lrfmat.Cut_Pile_Up)))'.*Corrections(lrfmat.Cut_Pile_Up),0.000001);
chi2_Pois = sum(2*(mu_val - spike_num(lrfmat.Cut_Pile_Up).*log(mu_val)));

mu_valG = max((TotalEnergy_Areas*(the_Direct_Component(~lrfmat.Cut_Pile_Up)+the_Reflected_Component(~lrfmat.Cut_Pile_Up)))'.*Corrections(~lrfmat.Cut_Pile_Up),0.000001); 
chi2_Gau = sum(2*(mu_valG - data(~lrfmat.Cut_Pile_Up).*log(mu_valG).*(Corrections(~lrfmat.Cut_Pile_Up))));


if dis_event_wall>0
    chi2 = chi2_Gau+chi2_Pois + lrfmat.max_prob;
else
    chi2 = chi2_Gau+chi2_Pois + lrfmat.max_prob + abs(50*dis_event_wall);
end

%%% For Debugging only

chisP = zeros(1, numel(lrfmat.Cut_Pile_Up));
chisP(lrfmat.Cut_Pile_Up) = 2*(mu_val - spike_num(lrfmat.Cut_Pile_Up).*log(mu_val));
chisP(~lrfmat.Cut_Pile_Up) = 2*(mu_valG - data(~lrfmat.Cut_Pile_Up).*log(mu_valG).*(Corrections(~lrfmat.Cut_Pile_Up)));

chiE = zeros(1, 61);
chiE(cutEv) = chisP;
Dir = zeros(1, 61); Ref = zeros(1, 61);
Dir(cutEv) = the_Direct_Component;
Ref(cutEv) = the_Reflected_Component;