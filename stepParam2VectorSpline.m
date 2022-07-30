% stepParam2VectorSpline.m
%
% Function that takes in legSteps struct and returns a vector of values for
%  a single field, during either stance or swing, at a specified times.
% Assigns parameter to midpoint of half step to which it belongs. Then uses
%  spline interpolation to fill in values for specified times. NaNs when
%  fly not walking.
% Alternative to stepParam2Vector()
%
% Helper function for: stepParamFiltlm and stepParamStepwiselm
%
% INPUTS:
%   legSteps - struct of leg parameter values, output of
%       extractLegStepsFromPData()
%   sampTime - times at which to return values
%   moveNotMove - struct of when fly is moving/not moving, output of
%       extractLegStepsFromPData()
%   whichLeg - 'R1', 'R2', 'R3', 'L1', 'L2', 'L3' string specifying which
%       leg
%   whichParam - name of field to get values for, as string
%   whichPhase - 'swing' or 'stance', for which phase of leg movement to
%       consider
%
% OUTPUTS:
%   outVec - column vector, output values of specified step parameter at
%       specified times. NaNs when fly not walking
%
% CREATED: 7/13/22 - HHY
%
% UPDATED:
%   7/19/22 - HHY
%
function outVec = stepParam2VectorSpline(legSteps, sampTime, moveNotMove, ...
    whichLeg, whichParam, whichPhase)

    % check that whichParam specifies a field of legSteps
    % fields of legSteps
    legStepsFields = fieldnames(legSteps);

    if ~any(strcmp(whichParam,legStepsFields))
        fprintf('%s is not a valid field of legSteps. Ending\n',whichParam);
        outVec = [];
        return;
    end

    % convert whichPhase to 1 for stance, -1 for swing, to match calls in
    %  legSteps
    if (strcmpi(whichPhase,'stance'))
        thisPhase = 1;
    elseif (strcmpi(whichPhase,'swing'))
        thisPhase = -1;
    else
        disp('Invalid value for whichPhase. Ending.');
        return;
    end

    % convert whichLeg to leg index
    whichLegInd = find(strcmpi(legSteps.legIDs.names, whichLeg));
    if isempty(whichLegInd)
        disp('Invalid value for whichLeg. Ending.');
        return;
    end
    thisLegInd = legSteps.legIDs.ind(whichLegInd);


    % loop through all steps for specified variable

    % get var, as steps, only for specified phase and leg - for correct
    % extraction of swing/stance, need to transpose matrix
    thisVarTemp = legSteps.(whichParam)';
    swingStanceLog = legSteps.stepSwingStance == thisPhase;
    swingStanceLog = swingStanceLog';
    thisVarAllLegs = thisVarTemp(swingStanceLog);
    % get specified leg only
    thisVar = thisVarAllLegs(legSteps.stepWhichLeg == thisLegInd);

    % get whether specified phase is first or second half step
    [swingStanceIndAllLegs, ~] = ind2sub(size(swingStanceLog), ...
        find(swingStanceLog));
    swingStanceInd = ...
        swingStanceIndAllLegs(legSteps.stepWhichLeg == thisLegInd);

    % get times for this leg, step start, mid, end times
    stepStartTimes = legSteps.stepT(legSteps.stepWhichLeg == thisLegInd,1);
    stepMidTimes = legSteps.stepT(legSteps.stepWhichLeg == thisLegInd,2);
    stepEndTimes = legSteps.stepT(legSteps.stepWhichLeg == thisLegInd,3);

    % loop through all steps, times for each step - midpoint of appropriate
    %  half step
    % initialize
    tStepTimes = zeros(size(stepStartTimes));

    for i = 1:length(stepStartTimes)
        % time assigned to this step parameter variable is at midpoint of
        %  relevant half step
        switch swingStanceInd(i)
            case 1 % first half step
                thisTime = (stepMidTimes(i) - stepStartTimes(i))/2 + ...
                    stepStartTimes(i);
            case 2 % second half step
                thisTime = (stepEndTimes(i) - stepMidTimes(i))/2 + ...
                    stepMidTimes(i);
        end
        tStepTimes(i) = thisTime;
    end

    % correct for multiple values at same time point
    if (length(tStepTimes) ~= length(unique(tStepTimes)))
        [tStepTimes, uniInd,~] = unique(tStepTimes);
        thisVar = thisVar(uniInd);
    end

    % using var vals at step times, interpolate to requested times
    outVec = interp1(tStepTimes, thisVar, sampTime);

    % when fly isn't moving, replace values in outVec with NaN

    % get not moving times
    notMoveStartTimes = moveNotMove.legT(moveNotMove.legNotMoveBout(:,1));
    notMoveEndTimes = moveNotMove.legT(moveNotMove.legNotMoveBout(:,2));

    outVecNotMoveInd = [];

    % loop through all not moving bouts, get indices of sampTime where fly
    %  not moving
    for i = 1:length(notMoveStartTimes)
        % include edges
        thisNotMoveInd = find((sampTime >= notMoveStartTimes(i)) & ...
            sampTime <= notMoveEndTimes(i));
        outVecNotMoveInd = [outVecNotMoveInd; thisNotMoveInd];
    end

    % replace outVec values with NaN
    outVec(outVecNotMoveInd) = nan;
end
