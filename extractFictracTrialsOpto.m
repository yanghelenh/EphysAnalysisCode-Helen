% extractFictracTrialsOpto.m
%
% Helper function for computeAvgFictracOpto_1Fly() that takes in one
%  trial's worth of data for one FicTrac parameter and divides it into
%  reps based on the opto stimulation
% Output appends trials into cell array given as input 
%
% INPUTS:
%   reps - cell array of size # NDs x # durations, to append into
%   fictracVar - FicTrac variable values
%   t - time for fictracVar
%   opto - struct of opto stimulation data
%   NDs - vector of all NDs 
%   durs - vector of all durations of stimulation
%   bwStimDur - scalar value of time between stimulations to consider, in
%       seconds, rep will be stimulation plus this time before and after
%
% OUTPUTS:
%   reps - same cell array as inputs, with additional reps added
%   durTs - cell array of same size as durs, with time at each sample point
%       for each duration; stimulation starts at 0 (in sec)
%
% CREATED: 3/29/22 - HHY
%
% UPDATED:
%   3/29/22 - HHY
%
function [reps, durTs] = extractFictracTrialsOpto(reps, fictracVar, t, ...
    opto, NDs, durs, bwStimDur)

    % get ND for this trial
    thisND = opto.stimParams.ndFilter;

    % get index corresponding to this ND
    thisNDInd = find(NDs == thisND);

    % if this ND isn't one of the specified ones, return
    if (isempty(thisNDInd))
        return;
    end

    % get ifi for FicTrac var
    ifi = mean(diff(t));

    % convert stimulation durations and between stimulation duration into 
    %  number of FicTrac samples
    dursNumSamps = floor(durs / ifi);
    bwStimNumSamps = floor(bwStimDur / ifi);

    % get times for each duration
    durTs = cell(size(durs));

    for i = 1:length(durs)
        thisDurSamps = 1:(dursNumSamps(i) + bwStimNumSamps * 2);

        durTs{i} = ((thisDurSamps - 1) * ifi) - (bwStimNumSamps * ifi);
    end

    % loop through all stimulations
    for i = 1:length(opto.stimCmdDurs)
        % get duration for this trial
        thisDur = opto.stimCmdDurs(i);

        % get index corresponding to this duration
        thisDurInd = find(durs == thisDur);

        % if this duration isn't one of the specified ones, go to next
        %  stimulation
        if (isempty(thisDurInd))
            continue;
        end

        % get FicTrac index corresponding to stimulation start
        stimStartInd = find(t >= opto.stimStartTimes(i),1,'first');

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

        % get this rep FicTrac data
        thisRepFictrac = fictracVar(repStartInd:repEndInd);

        % check for NaNs in FicTrac data; if present, skip this rep
        if (any(isnan(thisRepFictrac)))
            continue;
        end

        % append this rep FicTrac data to data cell array
        % each rep is a row
        reps{thisNDInd,thisDurInd} = [reps{thisNDInd,thisDurInd}; ...
            thisRepFictrac'];
    end
end