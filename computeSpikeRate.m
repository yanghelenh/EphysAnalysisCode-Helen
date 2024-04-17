% computeSpikeRate.m
%
% Function that takes in vector of indicies for when spikes occured and a
%  timing vector to compute spike rate as the inverse of the
%  inter-spike-interval
%
% INPUTS:
%   spikeInd - vector of indicies for when spikes occured
%   t - time at each index, in seconds
%
% OUTPUTS:
%   spikeRate - vector of spike rate (same size as t)
%
% CREATED: 9/14/20 - HHY
%
% UPDATED:
%   9/14/20 - HHY
%   1/30/24 - HHY - spike rate empty vector if no spikes
%
function spikeRate = computeSpikeRate(spikeInd, t)
    
    % times at which spikes occured
    spikeTimes = t(spikeInd);
    
    % inter-spike-interval
    intSpikeInt = gradient(spikeTimes);
    
    % inverse of the inter-spike-interval
    invIntSpikeInt = ones(size(intSpikeInt))./intSpikeInt;
    
    % interpolate vector to same size as time vector, instead of same size
    %  as spikeInd vector
    try
        spikeRate = interp1(spikeTimes, invIntSpikeInt, t);
    catch ME
        spikeRate = [];
    end
    
end