% getSpikeRateSteps.m
%
% Function that returns spike rate during each leg step. Can specify delay
%  between spikes and behavior (spike rate x sec before or after step).
% Spike rate computed as number of spikes during step divided by duration
%  of step in seconds.
%
% INPUTS:
%   spikeTimes - times (in sec) when all spikes occured, as column vector
%   stepTimes - start and end times (in sec) of all steps, as matrix of
%       number of steps x 2 (start is col 1, end is col 2)
%   tDelay - time (in sec) to shift spike and step times relative to each
%       other. Negative times are ephys before behavior. Scalar
%
% OUTPUTS:
%   spikeRate - spike rate (in Hz) for each step, as column vector of
%       length number of steps
%
% CREATED: 10/6/22 - HHY
%
% UPDATED:
%   10/6/22 - HHY
%
function spikeRate = getSpikeRateSteps(spikeTimes, stepTimes, tDelay)

    % incorporate time offset b/w ephys and behavior
    spikeTimesDelay = spikeTimes - tDelay;

    % preallocate spikeRate
    spikeRate = zeros(size(stepTimes,1),1);
    
    % loop through all steps
    for i = 1:size(stepTimes,1)
        % get step start, mid, end times
        stepStartT = stepTimes(i,1);
        stepEndT = stepTimes(i,2);

        % for step, figure out how many spikes 
        numSpikes = sum((spikeTimesDelay >= stepStartT) &...
            (spikeTimesDelay < stepEndT));

        % convert to spike rate
        spikeRate(i) = numSpikes / (stepEndT-stepStartT);
    end
end