function [xy, Status] = Guess_xy(data, PMT_r, the_sextant, topchsinner, stretch)

%This is Vladimir's function xy
for i=topchsinner,
    SextantAmplitude(i) = sum(data(the_sextant(i,:)));
end

[SextantAmplitudeArr, SextantAmplitudeid]  = sort(SextantAmplitude,'descend');
id0 = SextantAmplitudeid(1)

DT = data(the_sextant(id0,:))
if min(DT) > 0,
    Status = id0;
    tt(1) = log(DT(2)/DT(5)) + 0.5*log(DT(3)/DT(4)) + 0.5*log(DT(7)/DT(6))
    tt(2) = (log(DT(3)/DT(6)) + log(DT(4)/DT(7)))*sqrt(3.)/2.
    PMT_r(id0,1:2)
    xy = stretch*0.5*tt(1:2) + PMT_r(id0,1:2)
else,
    Status = 0;
    xy = PMT_r(id0,1:2);
end
