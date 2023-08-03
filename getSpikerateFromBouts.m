% getSpikerateFromBouts.m
%
% Helper function for saveSpikerate_bouts() that takes in indices of
%  ephys spiking data and returns the spike rate aligned to the velocity peak
% To allow averaging over different bouts and the inconsistent frame rate,
%  interpolate to specified frame rate
%
% INPUTS:
%   ephysSpikes - struct of spiking data
%   spikerateParams - struct of parameters, directly from
%     saveSpikerate_bouts()
%       maxDuration - time in seconds to consider on each side 
%       interpFrameRate - frame rate to interpolate to, in Hz
%   peakTime - start times of bout peaks
%   boutStartTime - start times of bout starts
%   boutEndTime - end times of bout ends
%   delay - time offset between behavior and ephys, neg for ephys before
%       behavior
%
% OUTPUTS:
%   spikerate - numTimePts x numBouts matrix for spikerate for bouts, 
%       aligned to yaw velocity peak
%   t - time vector for spike rate, interpolated
%
% CREATED: 7/17/23 - HHY
%
% UPDATED:
%   7/17/23 - HHY
%   8/2/23 - HHY - fix comments
%
function [spikerate, t] = getSpikerateFromBouts(ephysSpikes, ...
    spikerateParams, peakTime, boutStartTime, boutEndTime, delay)

    % get number time pts to either side of peak, after interpolation
    maxNumFrames = floor(spikerateParams.maxDuration * ...
        spikerateParams.interpFrameRate);
    % number of bouts
    numBouts = length(peakTime);
    % interframe interval, for interpolation
    ifi = 1/spikerateParams.interpFrameRate;

    % preallocate
    spikerate = nan(maxNumFrames * 2 + 1,numBouts);


    % loop through all bouts
    for i = 1:numBouts
        % zero time vector, so that yaw velocity peaks are at t = 0
        % account for delay
        tOrig = ephysSpikes.t - peakTime(i) - delay;

        % get time for bout start and end, zeroed
        boutStartT = boutStartTime(i) - peakTime(i);
        boutEndT = boutEndTime(i) - peakTime(i);

        % get time vector for interpolation - only consider +/- maxDuration
        %  around peak
        % doing it this way keeps 0 at 0
        newTDur = (maxNumFrames / spikerateParams.interpFrameRate);


        newTHalf1 = 0:ifi:newTDur;
        newTHalf1 = fliplr(newTHalf1) * -1;

        newTHalf2 = 0:ifi:newTDur;

        newT = [newTHalf1 newTHalf2(2:end)];

        % get interpolated leg X and Y position
        interpSpikerate = interp1(tOrig, ...
            ephysSpikes.spikeRate, newT,'spline');

        % filter for bout start and end
%         boutStartLog = newT < boutStartT;
%         interpSpikerate(boutStartLog) = nan;
% 
%         boutEndLog = newT > boutEndT;
%         interpSpikerate(boutEndLog) = nan;

        % add to output matrix
        spikerate(:,i) = interpSpikerate';
    end

    % get time vector for output matrix
    tHalf1 = fliplr((0:maxNumFrames) * ifi * -1);
    tHalf2 = (1:maxNumFrames) * ifi;

    t = [tHalf1 tHalf2]';    
end