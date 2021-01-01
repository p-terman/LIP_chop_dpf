function chi2 = MLM_Matricial(p_min, data, ML_Methode,pulse_area_eq, cutEv, lrfmat)

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
    TotalEnergy = p_min(3);
else
    TotalEnergy = pulse_area_eq;
    rho = [sqrt((p_min(1)-lrfmat.PMT_r(5,1)).^2 + (p_min(2)-lrfmat.PMT_r(5,2)).^2) sqrt((p_min(1)-lrfmat.PMT_r(32,1)).^2 + (p_min(2)-lrfmat.PMT_r(32,2)).^2)]';
    the_Direct_Component_NW = [lrfmat.Rho(5, round(10*(60-min(rho(2),59)).^2))' lrfmat.Rho(32, round(10*(60-min(rho(2),59)).^2))']';
    %the_Reflected_Component_NW = lrfmat.Amp(lrfmat.notwor',para_dis_event_wall).*exp(-rho.*lrfmat.Tau(lrfmat.notwor',para_dis_event_wall))+ lrfmat.Amp2(lrfmat.topchs(lrfmat.notwor')).* exp(-rho.*lrfmat.Gam2(lrfmat.topchs(lrfmat.notwor')))
    the_Reflected_Component_NW = [lrfmat.Amp(5,para_dis_event_wall).*exp(-rho(1).*lrfmat.Tau(5,para_dis_event_wall))+ lrfmat.Amp2(5).* exp(-rho(1).*lrfmat.Gam2(5))  lrfmat.Amp(32,para_dis_event_wall).*exp(-rho(2).*lrfmat.Tau(32,para_dis_event_wall))+ lrfmat.Amp2(32).* exp(-rho(2).*lrfmat.Gam2(32))]';
    TotalEnergy = TotalEnergy/(1-sum(the_Direct_Component_NW +  the_Reflected_Component_NW));
end

%% Maximum Likelihood method starts here

mu_val = max(repmat((TotalEnergy*(the_Direct_Component+the_Reflected_Component))'.*lrfmat.lrf_iq.QE.QE_Values(cutEv), [25 1]), 0);
data_val = repmat(data(cutEv).*lrfmat.lrf_iq.QE.QE_Values(cutEv), [25 1]);

Val_Prob = -2*log(sum(poisspdf(ML_Methode.PEAK_AREA_CutEv,mu_val).*ML_Methode.GaussianProb, 1));

chi2 = sum(Val_Prob)/(ML_Methode.Normalization_Factor)+ML_Methode.InitalPoint;

if dis_wall>0
    chi2 = chi2;
else
    chi2 = chi2+abs(50*dis_wall);
end




% The following code is for debugging
%disp(sprintf('%.2f\t%.2f\t%.2f\t%.2f\t', p_min(1), p_min(2),TotalEnergy, chi2));
%g = 0;
%for k = 1:61
%    if(cutEv(k)==1)
%        g = g +1;
%        disp(sprintf('%d\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\t%.2f\n', k, (TotalEnergy*(the_Direct_Component(g)+the_Reflected_Component(g))-data(k)').^2./(SD(g).*SystError), TotalEnergy*the_Direct_Component(g), TotalEnergy*the_Reflected_Component(g),data(k), SD(g)));
%    end
%end