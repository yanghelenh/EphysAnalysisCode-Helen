% filtEphysNoSpikes.m
%
% Median filters ephys trace to remove spikes. 
%
% INPUTS:
%   voltage - voltage trace
%   t - time at each sample, in seconds
%   filtOrder - order of filter, in seconds; i.e. number of seconds to take
%       median over
%
% OUTPUTS:
%   filtVolt - filtered voltage trace
%
% CREATED: 9/14/20 - HHY
%
% UPDATED:
%   9/14/20 - HHY
%
function filtVolt = filtEphysNoSpikes(voltage, t, filtOrder)
    
    % find inter-sample interval (i.e. inverse of sample rate), in seconds
    isi = median(diff(t));
    
    % filter order in samples, rounded to whole number
    filtOrderSamp = round(filtOrder / isi);
    
    filtVolt = medfilt1(voltage, filtOrderSamp);

end