% mexfunction [area_diff s2area max_s1area] = ...
%     S2Filter(data, s2window, s1window)
% 
% S2Filter passes the s2 filter over the input waveform -- the output is
% used to identify S2 pulses in the input waveform.
%
% algorithm:  Two box-filters are passed over the waveform, one S2-sized
% (default 2 us), and one with S1-sized (default, 200 ns).  The 'area_diff'
% at a given point is the area within the S2-box centered on that point,
% minus the maximum of the areas of all S1-boxes contained within the
% S2-box centered on that point.  ie, the diff is the area of the pulse
% minus the area of it's highest 200 ns (for pulses up to 2 us long).  For
% pulses less than 200 ns long, this comes to zero, and for an S1 with some
% noise after, it will be very small, including in the area_diff only the
% noise after the main pulse.
%
% Inputs:
%       data            -   the waveform to be searched for S2's -- usually the
%                             sum of all pmt's (after baseline flattening)
%       s2window        -   the width, in samples, of the S2 box, default 200
%       s1window        -   the width, in samples, of the S1 box, default 20
%
% Outputs:
%       area_diff       -   the area_diff at all points in the waveform
%                             (can't be calculated at the edges, so zeros
%                             are filled in.)
%       s2area          -   the area in the S2 box centered at that point
%                             in the waveform
%       max_s1area      -   the maximum area of the S1 boxes contained in
%                             the S2 box centered at that point in the
%                             waveform
%
%  CED, 08/31/06, revised for LUX 02/12/08

disp('S2Filter is a mex file -- if you see this, you need to compile S2Filter.c');
