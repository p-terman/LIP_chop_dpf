function height = gain_to_height(gain)
% height = gain_to_height(gain)
% 
% Gives the mean pulse height (mV) for a given PMT gain

height = 4.9e-7*gain - 0.084; % Empirical formula

