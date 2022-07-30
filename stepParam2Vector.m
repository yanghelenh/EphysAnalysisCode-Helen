% stepParam2Vector.m
%
% Function that takes in legSteps struct and returns a vector of values for
%  a single field, during either stance or swing, at a specified times.
% Holds the parameter value for the entire duration of the step (swing +
%  stance). Since steps are defined max-min-max position (forward walking, 
%  swing then stance), swing will generally fill forward in time while
%  stance will generally fill back.
%
% Helper function for: stepParamFitlm.m and stepParamStepwiselm.m
%
% INPUTS:
%   legSteps - struct of leg parameter values, output of
%       extractLegStepsFromPData()
%   sampTime - times at which to return values
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
% CREATED: 7/8/22 - HHY
%
% UPDATED:
%   7/8/22 - HHY
%
function outVec = stepParam2Vector(legSteps, sampTime, ...
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

    % get times for this leg, step start times, step end times (no need for
    %  mid times)
    stepStartTimes = legSteps.stepT(legSteps.stepWhichLeg == thisLegInd,1);
    stepEndTimes = legSteps.stepT(legSteps.stepWhichLeg == thisLegInd,3);

    % loop through all steps, get values at each time in the input time
    %  vector
    % when there are gaps in time, NaN
    % initialize
    outVec = nan(size(sampTime));

    for i = 1:length(stepStartTimes)
        % set output value to value for this step for indices of sampTime
        %  corresponding to this start and end time
        % greater than or equal to for start; less than only for end
        outVec((sampTime>=stepStartTimes(i)) & ...
            (sampTime < stepEndTimes(i))) = thisVar(i);
    end
end
