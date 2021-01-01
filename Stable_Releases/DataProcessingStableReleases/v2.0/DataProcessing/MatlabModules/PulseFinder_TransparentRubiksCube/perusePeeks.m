%
% usage: 
%   >> [afTiming hfTiming pp_areas flags] = perusePeeks(eventSamples,eventRecord,parameters);
%
% raison d'etre:
%   We wish to know the temporal extent (the edges) of pulses. Further analysis 
%   is predicated on having a clear definition of the edges. At the same time, we
%   recognize that features may be present just beyond the edges (in the "wings")
%   and we wish to know about this as well, in order to fully characterize the
%   pulse.
%                 ___
% description: ---| |---
%   A sliding box filter is passed over the event record. The width should be
%   chosen to be larger than the largest expected (physical) pulse. BY computing 
%   the convolution of the box and the eventRecord (cumulative sum in box), the box 
%   aligns itself to enclose the maximum possible area. Then, pulse edges and
%   and pulse areas are returned, for features INSIDE the box. Additionally, the 
%   "wing" areas before and after the box are also returned.
%
%   algorithm update in svn Rev 3777, 26 June 2013:
%     afTiming edges are now determined by starting from the maximum point in the box
%     region, work backwards in time calculating a smoothed/averaged bin content 
%     (nLookAhead + nLookBehind samples, either side of current sample).
%     If this average value stays below noiseThreshold (phe) for
%     at least maximumGap (samples), then stop and count this as the leading edge
%     of the pulse.
%
% inputs:
%   eventSamples = nx1 vector containing sample values, which may contain gaps (due to pod)
%   eventRecord  = nx1 vector containing event data, summed accross channels, corresponding to 
%                  eventSamples. Vector time base is 10 ns samples.  Agnostic about y units, 
%                  but likely this is phe/sample. IMPORTANT: vector is nx1. As a safety, use 
%                  eventRecord(:) to force nx1 dimensionality
%
%   parameters   = [preBoxSamples fullBoxSamples postBoxSamples edgeFraction]
%      preBoxSamples    = pre pulse region to check, before fullBox, sugest 400
%      fullBoxSamples   = width of fullBox, suggest 50
%      postBoxSamples   = post pulse region to check, after fullBox, suggest 50
%      edgeFraction     = fraction of pulse area used to define pulse edge. For example,
%                         edgeFraction=0.05 returns timing values which mark 5% and 95% 
%                         edges, by area, suggest 0.01
%
%      skinnyBoxSamples = width of skinnyBox used to calculate amis1_fraction, suggest 10
%      maximumGap       = number sequential samples consistent with baseline noise required to end pulse, suggest 50
%      nLookAhead       = number samples to sum in forward direction wrt scan direction, suggest 2
%      nLookBehind      = number of samples to sum in backwards direction wrt scan direction, suggest 2
%      noiseThreshold   = threshold below which we consider the smoothed signal to be just baseline noise, suggest 0.15
%
%
% outputs:
%   afTiming     = [aft_pulse_start_samples aft_t0_samples aft_t1_samples aft_t2_samples aft_pulse_end_samples] 
%                   Units are samples, and times are based on area fractions (rather 
%                   than height fractions)
%   hfTiming     = [hft_t0_samples hft_t10l_samples hft_t50l_samples hft_t1_samples hft_t50r_samples hft_t10r_samples hft_t2_samples]
%                   Units are samples, times are the points at which the pulse returns to baseline
%   pp_areas     = [preFullBoxArea fullBoxArea postFullBoxArea] 
%                   Area corresponds to input units (probably phe/sample*samples = phe)
%   flags        = [(fullBoxSamples > eventRecord)
%                   amis1    
%                  ]
%                     
%                  amis1, the 2nd flag value returned, is a fraction whos value =1 for a
%                  perfect s1, and 0 for a perfect s2.
%        
%         
% versioning:
% 120927 pfs - distilled from previous algorithms, for LUX
% 121110 pfs - the intermittent random stack crash issue was traced to the pointers for the 
%              output. now using dynamic memory allocation as an attempted fix. have
%              verified stable operation, details in LUX Matlab Bug Tracker incident 0000050
% 121203 pfs - stream-lined and stable. also, boxStart is now the right-most sample which 
%              satisfies the algorithm (rather than the 1st sample which satisfies it), 
%              as per CHF suggestion.
% 130322 pfs - new input structure to accomodate podsum
% 130328 pfs - now returns hfTiming as well. Also, employs a variable box width to avoid
%              chopping off the edges of pulses.
% 130626 JD & pfs - algorithm update (require additional 4 input params)
% 130702 pfs - enforce minimum fullBoxSamples=3
% 130705 JD  - fix issue with indexing when calculating skinny box area and deal with case 
%              where skinny box width greater than pulse width
% 130706 JD & KOS - minor fix to ensure min fullBoxSamples=3 never takes pulse beyond waveform
% 130708 pfs & KOS - minimum box width set = skinnyBoxSamples, and also now setting left edge
%                    of pulse -1 sample earlier (based on by-eye assesment)
% 130811 pfs & JD - additional aft_tlx_samples and aft_trx_samples timing variables 
% 131011 JD - made tlx and trx timing fractions configurable using xml settings
%
disp('*** if you are reading this message, you need to compile perusePeeks.c ***');
disp('(thus spoke perusePeeks.m, the help file)');

