% extractTrialsIinj.m
%
% Helper function for computeLegFictracEphysIinj_1Fly() that takes in one
%  trial's worth of data for one parameter and divides it into reps based 
%  on current injection
% Output appends trials into cell array given as input 
% Has option to align position data by subtracting value at stimulation
%  start
% Adaptation of extractTrialsOpto()
%
% INPUTS:
%   reps - cell array of size # NDs x # durations, to append individual 
%       trials into
%   repsIinjTimes - cell array of size # NDs x # durations, to keep
%       track of opto stim start times for each trial, to append into
%   repsPDatNames - cell array of size # NDs x # durations, to keep track of
%       pData name for each rep, to append into
%   var - input data variable values
%   t - time for input data parameter
%   iInj - struct of current injection data
%   amps - vector of all amplitudes
%   durs - vector of all durations of stimulation
%   bwStimDur - scalar value of time between stimulations to consider, in
%       seconds, rep will be stimulation plus this time before and after
%   pDataName - name of pData file for this trial
%   norm2StimStart - logical for whether to normalize to value at
%       stimulation start
%
% OUTPUTS:
%   reps - same cell array as inputs, with additional reps added
%   repsIinjTimes - same cell array as inputs, with additional current 
%       injection start times added
%   repsPdatName - same cell array as inputs, with additional pData names
%       for each trial added
%   durTs - cell array of same size as durs, with time at each sample point
%       for each duration; stimulation starts at 0 (in sec)
%
% CREATED: 6/29/22 - HHY
%
% UPDATED:
%   6/29/22 - HHY
%   9/22/22 - HHY - check to make sure stim start time is not later than
%       length of trial
%
function [reps, repsIinjTimes, repsPDatNames, durTs] = extractTrialsIinj(...
    reps, repsIinjTimes, repsPDatNames, var, t, iInj, amps, durs, ...
    bwStimDur, pDataName, norm2StimStart)

    % check that var is column vector; if not, make it so
    if (isrow(var))
        var = var';
    end

    % get ifi for behavior var
    ifi = mean(diff(t));

    % convert stimulation durations and between stimulation duration into 
    %  number of behavior samples
    dursNumSamps = floor(durs / ifi);
    bwStimNumSamps = floor(bwStimDur / ifi);

    % get times for each duration
    durTs = cell(size(durs));

    for i = 1:length(durs)
        thisDurSamps = 1:(dursNumSamps(i) + bwStimNumSamps * 2);

        durTs{i} = ((thisDurSamps - 1) * ifi) - (bwStimNumSamps * ifi);
    end

    % loop through all stimulations
    for i = 1:length(iInj.amps)

        % get amplitude for this rep
        thisAmp = iInj.amps(i);
        % get duration for this rep
        thisDur = iInj.durs(i);

        % get index corresponding to this amplitude
        thisAmpInd = find(amps == thisAmp);
        % get index corresponding to this duration
        thisDurInd = find(durs == thisDur);

        % if this duration or amplitude isn't one of the specified ones, 
        %  go to next stimulation
        if (isempty(thisDurInd) || isempty(thisAmpInd))
            continue;
        end

        % get index into var corresponding to stimulation start
        stimStartInd = find(t >= iInj.startTimes(i),1,'first');

        % if startTime is too late (not part of trial)
        if isempty(stimStartInd)
            continue;
        end

        % index for rep start (stim start with b/w stim period appended)
        repStartInd = stimStartInd - bwStimNumSamps;
        % index for rep end (stim start + stim duration + b/w stim period)
        repEndInd = stimStartInd + dursNumSamps(thisDurInd) + bwStimNumSamps - 1;

        % check that start index is valid (1 or greater), otherwise, skip
        %  this rep
        if (repStartInd < 1)
            continue;
        % check that end index is valid (not greater than length of whole
        %  trial; otherwise, skip this rep
        elseif (repEndInd > length(t))
            continue;
        end

        % get this rep var data
        thisRepBeh = var(repStartInd:repEndInd);
        thisRepStartVal = var(stimStartInd);

        % check if var data should be normalized to value at start of
        % stimulation (subtract out that value)
        if (norm2StimStart)
            thisRepBeh = thisRepBeh - thisRepStartVal;
        end

        % append this rep data to data cell array
        % each rep is a row
        reps{thisAmpInd,thisDurInd} = [reps{thisAmpInd,thisDurInd}; ...
            thisRepBeh'];

        % append this rep iInj start time to appropriate cell array
        repsIinjTimes{thisAmpInd,thisDurInd} = [...
            repsIinjTimes{thisAmpInd,thisDurInd}; ...
            iInj.startTimes(i)];

        % append this rep pData name to appropriate cell array
        repsPDatNames{thisAmpInd,thisDurInd}{end + 1} = pDataName;
    end
end