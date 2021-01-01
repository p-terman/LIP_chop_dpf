function the_sextant = Find_Sextant(pmt, pmt_pos, sextant_arrangement)

Position_PMT_Sextant(1,:) = pmt_pos(find(sextant_arrangement == pmt),:);
Pmt_Enum =  0:1:5;

%Get The positions of the PMTs
%This should be replaced for something more efficient ATTENTION TO
%THE GODDAM 6
Position_PMT_Sextant(Pmt_Enum+2,1) = Position_PMT_Sextant(1,1)+6*cos(Pmt_Enum*pi/3); 
Position_PMT_Sextant(Pmt_Enum+2,2) = Position_PMT_Sextant(1,2)+6*sin(Pmt_Enum*pi/3);
Position_PMT_Sextant(Pmt_Enum+2,3) = Position_PMT_Sextant(1,3);

%There is any other logical way to do this. It works but not the
%way I like it

for j = 1:7,   
    try
        sextant_position(j)=find(int8(pmt_pos(:,1))== int8(Position_PMT_Sextant(j,1)) & int8(pmt_pos(:,2))==int8(Position_PMT_Sextant(j,2)) & int8(pmt_pos(:,3))==int8(Position_PMT_Sextant(j,3)));
    catch err
        sextant_position(j) = 140;
    end
end;
sextant_arrangement(140)=0;
the_sextant = sextant_arrangement(sextant_position(sextant_position>0));