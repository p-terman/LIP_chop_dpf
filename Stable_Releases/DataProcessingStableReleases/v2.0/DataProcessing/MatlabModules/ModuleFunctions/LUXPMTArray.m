function [pmt_pos_out sextant_arrangement] = LUXPMTArray(fignum)
% [pmt_pos sextant_arrangement] = LUXPMTArray(fignum)
% 
% Creates a 1x2 plot with top and bottom PMT arrays and corresponding
% position numbers (1-122)
% 
% Inputs:
%  fignum: if specified, will draw PMT array in figure #[fignum] [default:
%          do not draw array]
% 
% Outputs:
%              pmt_pos: 122x3 matrix of PMT centers in cm (LUXSim standard 
%                       coordinates)
%  sextant_arrangement: 122 vector to map to go from LUXFirstPass output
%                       into LUXSim PMT coordinates, i.e.
%                       pmt_hits_luxsim_coords = ...
%                          pmt_hits_from_LUXFirstPass(sextant_arrangement);
%                       
% 
% DCM 2010-03-10
% JRV 2010-10-10 - Updated to have new PMT sextant numbers
% 


%% Create PMT position array

sextant_arrangement_top = [ 1, 54, 53, 52, 51, 2, 5, 57, 56, 55, 44, ...
    3, 6, 8, 59, 58, 47, 43, 4, 7, 9, 10, 60, 49, 46, 42, 11, 15, 18, ...
    20, 121, 50, 48, 45, 41, 12, 16, 19, 30, 40, 39, 37, 34, 13, 17, ...
    28, 29, 38, 36, 33, 14, 25, 26, 27, 35, 32, 21, 22, 23, 24, 31 ];

sextant_arrangement_bottom = [ 61, 114, 113, 112, 111, 62, 65, 117, ...
    116, 115, 104, 63, 66, 68, 119, 118, 107, 103, 64, 67, 69, 70, ...
    120, 109, 106, 102, 71, 75, 78, 80, 122, 110, 108, 105, 101, 72, ...
    76, 79, 90, 100, 99, 97, 94, 73, 77, 88, 89, 98, 96, 93, 74, 85, ...
    86, 87, 95, 92, 81, 82, 83, 84, 91 ];

sextant_arrangement = [sextant_arrangement_top sextant_arrangement_bottom];

pmt_width = 5.7; % cm
pmt_spacing = 0.3; % cm

pmt_top_pos = 12.8+11.6; % cm
pmt_bot_pos = -42; % cm

NUM_ROWS = 9;
pmts_per_row = [ 5, 6, 7, 8, 9, 8, 7, 6, 5 ];
start_x = [ -2*(pmt_width+pmt_spacing),
    -2.5*(pmt_width+pmt_spacing),
    -3*(pmt_width+pmt_spacing),
    -3.5*(pmt_width+pmt_spacing),
    -4*(pmt_width+pmt_spacing),
    -3.5*(pmt_width+pmt_spacing),
    -3*(pmt_width+pmt_spacing),
    -2.5*(pmt_width+pmt_spacing),
    -2*(pmt_width+pmt_spacing) ];

pmt_pos = zeros(122,3);

pmt_ctr = 1;
for ii_j = length(pmts_per_row)-1:-1:0
    for ii_i = 0:pmts_per_row(ii_j+1)-1
        pmt_pos(pmt_ctr,1:2) = ...
            [start_x(ii_j+1)+ii_i*(pmt_width+pmt_spacing), ...
            (pmt_width+pmt_spacing)*(ii_j-(NUM_ROWS-1)/2)*sin(pi/3)];
        pmt_ctr = pmt_ctr+1;
    end
end

pmt_pos(62:end,1:2) = pmt_pos(1:61,1:2);

pmt_pos(1:61,3) = pmt_top_pos;
pmt_pos(62:end,3) = pmt_bot_pos;

if nargout >= 1
    pmt_pos_out = pmt_pos;
end



%% Draw

if nargin >= 1 && ~isempty(fignum)
    
    figure(fignum); clf;
    
    subplot(1,2,1);

    tt = 0:pi/64:2*pi;
    xx = pmt_width./2 .* cos(tt);
    yy = pmt_width./2 .* sin(tt);

    for ii_pmt=1:61

        fill(pmt_pos(ii_pmt,1)+xx, ...
            pmt_pos(ii_pmt,2)+yy, ...
            [1 1 1]);
        hold on
        text(pmt_pos(ii_pmt,1), pmt_pos(ii_pmt,2), num2str(sextant_arrangement(ii_pmt)),'HorizontalAlignment','center','FontSize',8);

    end
    hold off

    axis equal
    axis off

    title('Top PMT Array');

    subplot(1,2,2);

    for ii_pmt=62:122

        fill(pmt_pos(ii_pmt,1)+xx, ...
            pmt_pos(ii_pmt,2)+yy, ...
            [1 1 1]);
        hold on
        text(pmt_pos(ii_pmt,1), pmt_pos(ii_pmt,2), num2str(sextant_arrangement(ii_pmt)),'HorizontalAlignment','center','FontSize',8);

    end
    hold off

    axis equal
    axis off

    title('Bottom PMT Array');
    
end

