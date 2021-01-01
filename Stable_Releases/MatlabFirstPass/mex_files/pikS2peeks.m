% mexfunction [pulse_times pulse_areas pulse_diffs trace_diffs] = ...
%     pikS2peeks(data, s2window, s1window, threshold, rel_area_thresh, x50_cutpoint, num_to_find)
% 
% pikS2peeks finds the n most likely S2 pulses in the summed waveform
% passed to it, and returns their timing info, their areas within a given
% window, and the 'diffs', or the parameter used to determine their
% S2-ness, both the maximum for each pulse and the value over the summed
% waveform.
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
%       threshold       -   the threshold, in phe, placed on area_diff to
%                             qualify as a potential S2.  default .5 (ie,
%                             pick up any two phe separated by between 200
%                             and 2000 ns)
%       rel_area_thresh -   the point at which a pulse is truncated (if it 
%                             hasn't already been ended):  If the area in 
%                             the s2window beyond this point is less than
%                             rel_area_thresh * pulse_diff, where
%                             pulse_diff is the maximum area_diff for the
%                             pulse, the pulse is truncated.  default 0
%                             (ie, cut off the pulse if there is no more
%                             area in the next 2 us, otherwise wait for one
%                             of the other thresholds)
%       x50_cutpoint    -   the maximum interval to include in the pulse --
%                             on the left and right, the pulse will extend
%                             to at most (t1 - (t1-t50l)*x50_cutpoint) and
%                             (t1 + (t50r-t1)*x50_cutpoint) respectively,
%                             where t1 is the time of maximum area_diff,
%                             and t50l/r are the times of area_diff half
%                             that at t1 (the pulse always extends at least
%                             out to t50lr, unless the rel_area_thresh
%                             cutoff is encountered). default 2.5
%       num_to_find     -   the maximum number of S2's to look for.
%                             Actually, 1 more than this number are found,
%                             but the smallest is continuously overwritten,
%                             so that only the first num_to_find are
%                             guarenteed to be the n with the largest
%                             maximum area_diff
%
% Outputs:
%       pulse_times     -   a double array of size [7, num_to_find+1],
%                             giving the following 7 times for each pulse
%                             (the num_to_find+1th entry is to be ignored).
%                               t0:  the point at which integration
%                                 should begin -- determined either by
%                                 rel_area_thresh, x50_cutpoint, or by when
%                                 the area_diff drops below threshold (plus
%                                 some buffer)
%                               t10l:  the point at which the area_diff
%                                 drops below threshold (if it comes after
%                                 t0, otherwise it's t0+1)
%                               t50l:  the point closest to t1 on the left
%                                 at which the area_diff is half that of
%                                 the maximum (if this point comes after
%                                 t10l -- otherwise it's t10l+1)
%                               t1:  the point of maximum area_diff
%                               t50r:  same as t50l, but on the right
%                               t10r:  same as t10l, but on the right
%                               t2:  the point at which integration should
%                                 end, determined same way as t0.
%       pulse_areas     -   the area in the s2window centered at t1
%       pulse_diffs     -   the area_diff at t1 (maximum area_diff for the
%                             pulse -- the output is sorted in this
%                             parameter.)
%       trace_diffs     -   the area_diff at all points in the waveform
%                             (can't be calculated at the edges, so zeros
%                             are filled in.)
%
%  CED, 08/31/06

