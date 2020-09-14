% computeSpikeRate.m
%
% Function that takes in vector of indicies for when spikes occured, a
%  timing vector, and a window size to compute the spike rate.
%
% INPUTS:
%   spikeInd - vector of indicies for when spikes occured
%   t - time at each index, in seconds
%   window - size of moving window, in seconds, to compute spike rate over
%
% OUTPUTS:
%   spikeRate - vector of spike rate (same size as t)
%
% CREATED: 9/14/20 - HHY
%
% UPDATED:
%   9/14/20 - HHY
%
function spikeRate = computeSpikeRate(spikeInd, t, winSize)

    % find inter-sample interval (i.e. inverse of sample rate), in seconds
    isi = median(diff(t));
    
    % get size of window, in samples
    winSizeSamp = round(winSize / isi);
    winSizeSec = winSizeSamp * isi; % window size in seconds, post rounding
    
    % vector of size t where 0 for no spike, 1 for spike
    spikeLog = zeros(size(t));
    spikeLog(spikeInd) = 1;
    
    % compute moving average over spike logical vector, using specified
    %  window size; divide by window size in seconds to get spike rate
    spikeRate = movmean(spikeLog, winSizeSamp);
    
    
    spikeTimes = t(spikeInd);
    intSpikeInt = diff(spikeTimes);
    
    spikeRate = ones(size(intSpikeInt))./intSpikeInt;
    
end