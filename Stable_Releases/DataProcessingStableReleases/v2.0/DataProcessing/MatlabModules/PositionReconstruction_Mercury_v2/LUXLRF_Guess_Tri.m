function [XXs YYs Status] = LUXLRF_Guess_Tri(data, PMT_r, the_sextant, mode, stretch)
    
%This is Vladimir's method to get an initial estimative of the
%event's position. It worked well for ZEPLIN, but it seems there
%are some problems for the LUX detector.
%
% [xy, Status] = LUXLRF_Guess_Tri(data, PMT_r, the_sextant, mode, stretch)
% 
% 
% Inputs:
%  data - the peak area for each PMT
%  PMT_r - the positions for each PMT
%  the_sextant - the sextant arrangement
%  mode - the number of sextants used for the estimative of the
%  position  (the only choices are 1 or 3) 
%  stretch - this variable is algorithm related. Just choose 1 if
%  you don't know how this works. You can play with it if you don't
%  have nothing to do.
%    
% Outputs:
%  
% sp - the new dataset with the important quantities
%
% 2013-01-25 CFPS v1.0
% 
    
    
function Point_XY = Point_UV(DT)
    Point_XY(:,1) = log(DT(:,2)./DT(:,5)) + 0.5*log(DT(:,3)./DT(:,4)) + 0.5*log(DT(:,7)./DT(:,6));
    Point_XY(:,2) = (log(DT(:,3)./DT(:,6)) + log(DT(:,4)./DT(:,7)))*sqrt(3.)/2.;
end

function chisq = XYComplet(Point, MLine, MColumn)
    chisq = sum((MLine(:,1).*Point(1)+MLine(:,2).*Point(2)-MColumn(:)).^2);
end

PMT_Radius = round(sqrt(PMT_r(1:122,2).^2 + PMT_r(1:122,1).^2)-PMT_r(1:122,3)+24.4000);
topchsinner = PMT_Radius < 20;

PMTsAll = 1:122;
PMTs = PMTsAll(topchsinner);
for i=PMTs,
    SextantAmplitude(i) = sum(data(the_sextant(i,:)));
end

[SextantAmplitudeArr, SextantAmplitudeid]  = sort(SextantAmplitude,'descend');


try
    minimum = min(min(data(the_sextant(SextantAmplitudeid(1:3),:))));
    if minimum == 0,
        mode = 1;
    end
catch err  
    mode = 1;
end

if mode == 1,

    
    id0 = SextantAmplitudeid(1);
    DT = data(the_sextant(SextantAmplitudeid(1),:));
    if min(DT) > 0,
        Status = id0;
        tt = Point_UV(DT);
        xy = stretch*0.5*3*tt(1:2) + PMT_r(id0,1:2);
    else,
        Status = 0;
        xy = PMT_r(id0,1:2);
    end

elseif mode == 3,    
     tt = Point_UV(data(the_sextant(SextantAmplitudeid(1:3),:)));
     XY = stretch*0.5*3*tt + PMT_r(SextantAmplitudeid(1:3),1:2);
     PMTCPos = PMT_r(SextantAmplitudeid(1:3),:);
     MLine = [PMTCPos(:,2)-XY(:,2) XY(:,1)-PMTCPos(:,1)];
     MColumn = PMTCPos(:,2).*XY(:,1)-PMTCPos(:,1).*XY(:,2);
     Initial_Guess = [mean(XY(:,1)) mean(XY(:,2))];
     xy = fminsearch(@(p_min) XYComplet(p_min, MLine, MColumn), Initial_Guess);
     Status = 3;

end

XXs = xy(1);
YYs = xy(2);

end