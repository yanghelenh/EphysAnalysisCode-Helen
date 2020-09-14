% detectSpikes.m
%
% Function that takes in voltage trace and returns indicies of when action
%  potentials occurred. Action potentials defined by threshold crossings in
%  first derivative of voltage (dv/dt). Adacent threshold crossings must be
%  at least refractoryPeriod length of time from each other.
% 
% Returns indicies of spike start.
%
% When two threshold crossings are too close together, assumes the latter
%  is wrong.
%
% INPUTS:
%   voltage - voltage trace
%   t - time trace
%   dvdtThresh - threshold to cross from below for start of spike
%   refractoryPeriod - length of time in seconds where another spike can't
%       occur
%
% OUTPUTS:
%   startInd - indicies for start of each action potential
%
% CREATED: 9/13/20 - HHY
%
% UPDATED:
%   9/13/20 - HHY
%   9/14/20 - HHY - change to only report spike start indicies, instead of
%       start, end, and peak (auto correction too hard for peak)
%
function startInd = detectSpikes(voltage, t, dvdtThresh, refractoryPeriod)

    % find inter-sample interval (i.e. inverse of sample rate), in seconds
    isi = median(diff(t));
    
    % find refractory period in samples
    refrPrdSamp = refractoryPeriod / isi;
    
    % find first derivative of voltage
    dvdt = diff(voltage) ./ diff(t);
    
    % find when first derivative of voltage exceeds threshold
    abvThresh = dvdt > dvdtThresh;
    
	% find indicies for when these threshold crossings occur
    spikeStartInd = find(diff(abvThresh) == 1); % 1st to exceed thresh
    
    % correct for taking difference
    spikeStartInd = spikeStartInd + 1;
    
    % correct for if there are threshold crossings too close together
    %  (within refractory period)
    spikeStartIndRefrCorr = rmRefrPrdThreshCross(spikeStartInd, ...
        refrPrdSamp);
    
    % output
    startInd = spikeStartIndRefrCorr;
end
