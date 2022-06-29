% preprocessIInj.m
%
% Function for processing current injection from iInj field of ephysData
% Extracts start and end times of each iInj step as well as its amplitude.
%  Also, logical for iInj step or not. And original iInj vector.
% For steps of 0 pA, uses known duration of time between steps to estimate
%  when this 0 pA step was presented. Saves parameters of iInj that were
%  originally saved in inputParams.
% Does not use ephysData.current because of contamination with action
%  potentals
% Returns iInj struct
%
% INPUTS:
%   daqData - struct of data from experimental DAQ, processed by
%       preprocessUserDaq()
%   daqOutput - struct of output signals sent by experimental DAQ,
%       processed by preprocessUserDaq() 
%   daqTime - vector of times corresponding to each sample point of daqData
%   inputParams - struct of input parameters, specific to particular
%       experiment type, saved in trial.mat file when experiment run
%
% OUTPUTS:
%   iInj - output struct with the following fields:
%       startTimes - start times of steps
%       endTimes - end times of steps
%       amps - amplitude of current steps, in pA
%       durs - duration of current steps, in seconds
%       onLogical - logical for whether iInj step was presented (1) or
%           not (0)
%       cmd - full commanded current sequence, in pA; ephysData.iInj
%       params - parameters, from inputParams
%
% CREATED: 6/28/22 - HHY
%
% UPDATED:
%   6/28/22 - HHY
%
function iInj = preprocessIInj(ephysData, inputParams)

    % get indices where iInj changes
    changeInds = find(diff(ephysData.iInj) ~= 0);

    % changeInd with addition of index 1 and end to capture last steps
    changeIndsW1End = [1; changeInds; length(ephysData.t) + 1];

    % times when iInj changes
    changeTimes = ephysData.t(changeInds + 1);

    ifi = median(diff(ephysData.t));

    % durations of each step between changes (add index 1 and end to
    %  capture first and last steps)
    bwChangeDurs = diff([ephysData.t(1); ...
        changeTimes; ephysData.t(end)+ifi]);

    % if 0 pA is a possible step amplitude, look for durations greater than
    %  the space duration, back-fill 0 pA steps
    if any(inputParams.iInjParams.allStepAmps == 0)
        % if steps are of different durations, extract dur of 0 pA step
        if (length(inputParams.iInjParams.stepDur) > 1)
            zeroStepInd = find(inputParams.iInjParams.allStepAmps == 0);
            zeroStepDur = inputParams.iInjParams.stepDur(zeroStepInd);
        else
            zeroStepDur = inputParams.iInjParams.stepDur;
        end

        % find all b/w changes that are too long
        longStepInd = find(bwChangeDurs > inputParams.iInjParams.spaceDur);

        % loop through all b/w changes that are too long
        for i = 1:length(longStepInd)
            thisBwChangeDur = bwChangeDurs(longStepInd(i));

            % number of 0 pA steps in this b/w change
            numSteps = round((thisBwChangeDur - inputParams.iInjParams.spaceDur) / ...
                (inputParams.iInjParams.spaceDur + zeroStepDur));

            % fill in step transitions (append to changeInds)
            thisBwChangeStartInd = changeIndsW1End(longStepInd(i));

            % add start and end indices of steps
            for j = 1:numSteps
                % start ind is space duration after last transition
                newStartInd = thisBwChangeStartInd + ...
                    round(inputParams.iInjParams.spaceDur / ifi) * j + ...
                    round(zeroStepDur / ifi) * (j - 1);
                % end ind is zero step duration after start ind
                newEndInd = newStartInd + round(zeroStepDur / ifi); 

                changeInds = [changeInds; newStartInd; newEndInd];
            end
        end

        % sort changeInds in ascending order (since appending doesn't
        %  account for sorting)
        changeInds = sort(changeInds);

        % update changeTimes with new changeInds
        changeTimes = ephysData.t(changeInds + 1);

        % durations of each step between changes (add index 1 and end to
        %  capture first and last steps)
        bwChangeDurs = diff([ephysData.t(1); ...
            changeTimes; ephysData.t(end)+ifi]);

        % new changeInd with addition of index 1 and end to capture last steps
        changeIndsW1End = [1; changeInds; length(ephysData.t) + 1];
    end

    % convert changeTimes to start and end times
    % note that step injections always start with space (in iInj code)
    firstSpaceInd = find(bwChangeDurs == inputParams.iInjParams.spaceDur,...
        1,'first');

    iInjStartTimes = changeTimes(firstSpaceInd:2:length(changeTimes));
    iInjEndTimes = changeTimes((firstSpaceInd + 1):2:length(changeTimes));

    iInjStartInds = changeInds(firstSpaceInd:2:length(changeInds));
    iInjEndInds = changeInds((firstSpaceInd + 1):2:length(changeInds));

    % if acquistion stops in middle of step, so more start times than end
    % times
    if (length(iInjStartTimes) > length(iInjEndTimes))
        iInjStartTimes = iInjStartTimes(1:length(iInjEndTimes));
        iInjStartInds = iInjStartInds(1:length(iInjEndInds));
    end

    % duration of steps
    iInjDurs = iInjEndTimes - iInjStartTimes;

    % amplitude of steps
    % preallocate
    iInjAmps = zeros(size(iInjDurs));


    %  loop through all steps, find amplitude as value at midpoint
    for i = 1:(length(iInjStartInds))
        stepStartInd = iInjStartInds(i);
        stepEndInd = iInjEndInds(i);

        stepMidInd = floor((stepEndInd - stepStartInd)/2 + stepStartInd);

        iInjAmps(i) = ephysData.iInj(stepMidInd);
    end

    % generate iInjOnLogical
    % preallocate
    iInjOnLogical = zeros(size(ephysData.iInj));

    % loop through all steps, flip to 1
    for i = 1:length(iInjStartInds)
        iInjOnLogical(iInjStartInds(i):iInjEndInds(i)) = 1;
    end

    % convert to logical array
    iInjOnLogical = logical(iInjOnLogical);

    % save into struct
    iInj.startTimes = iInjStartTimes;
    iInj.endTimes = iInjEndTimes;
    iInj.amps = iInjAmps;
    iInj.durs = iInjDurs;
    iInj.onLogical = iInjOnLogical;
    iInj.cmd = ephysData.iInj;
    iInj.params = inputParams.iInjParams;
end