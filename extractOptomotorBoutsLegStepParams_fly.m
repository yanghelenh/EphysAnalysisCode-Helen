% extractOptomotorBoutsLegStepParams_fly.m
%
% Function that extracts turning bouts in response to optomotor stimulus
%  and saves leg step parameters during these bouts, relative to yaw
%  velocity peak.
% Adaptation of saveBallLegStepParamCond_bouts(), which extracts turning
%  bouts in flies walking on the ball, with conditioning, but not relative
%  to any visual stimulus
% Meant to operate on all pData files for 1 fly 
%
% INPUTS:
%   cond - struct of conditions for yaw velocity peak, if multiple 
%     conditions, treats it as AND
%       whichParam - cell array (even if 1 element) on which fictracSmo
%           field to condition on, one for each condition
%       cond - cell array of strings to condition on, for eval(); same size
%           as whichParam
%       legs - which leg (R1, R2, R3, L1, L2, L3) for stepFwdBool, [] for 
%           FicTrac
%       turnDur - 2 element vector [minTurnDuration maxTurnDuration] to
%           specify the min and max duration of the turning bout for it to
%           be included
%       minYawThresh - minimum yaw velocity to define start and end of bout
%   fwdVelCond - struct of conditions on forward velocity to apply to
%     turning bout
%       initVel - 2 element vector defining [min max] range of acceptable
%           initial forward velocities (at bout start)
%       change - 2 element vector defining [min max] range of acceptable
%           changes in forward velocity, peak - start 
%   pkTimeRelOpto - 2 element vector defining time relative to optomotor
%       stimulus yaw velocity peak has to fall into 
%       [after opto start, after opto end] (negative is before, positive is
%       after)
%   maxNumSteps - number of steps to each side of peak to consider as part
%       of bout (max bout length is this x2 + 1)
%   pDataFNames - cell array of pData file names or [] if select through
%       GUI
%   pDataPath - full path to pData directory
%   saveFilePath - directory in which to save output file
%   saveFileSuffix - suffix to append to name of first pData file, to
%       generate output file name
%
% OUTPUTS:
%   none, but saves output file in saveFilePath
%       selLegSteps - struct of aligned step parameters, where each one is 
%           numSteps x numLegs x 2 x numBouts matrix 
%       selStanceParams - struct of aligned step parameters during stance
%           only, where each one is numSteps x numLegs x numBouts matrix
%       selSwingParams - same as selStanceParams, but for swing
%       selNDs - numBouts x 1 vector of ND for each bout
%       selOptoVel - numBouts x 1 vector of optomotor stim velocity for
%           each bout
%       pkTimeOpto - numBouts x 1 vector for time of peak, relative to
%           optomotor stimulus (0 is start of optomotor stimulus)
%       pkSwingStance - numLegs x numBouts matrix for whether peak is
%           during swing or stance for each leg and bout
%       numBouts - total number of bouts (max n for stepParam (if no NaNs))
%       boutPeakVel - struct for velocity values at each bout peak
%           yaw - vector of length numBouts for peak yaw velocity
%           fwd - vector of length numBouts for peak forward velocity
%           slide - vector of length numBouts for peak lateral velocity
%       boutStartVel - struct for velocity values at each bout start
%           yaw - vector of length numBouts for start yaw velocity
%           fwd - vector of length numBouts for start forward velocity
%           slide - vector of length numBouts for start lateral velocity
%       boutEndVel - struct for velocity values at each bout end
%           yaw - vector of length numBouts for end yaw velocity
%           fwd - vector of length numBouts for end forward velocity
%           slide - vector of length numBouts for end lateral velocity
%       boutPeakVelSmo - as boutPeakVel, but for fictracSmo values
%       boutStartVelSmo - as boutStartVel, but for fictracSmo values
%       boutEndVelSmo - as boutEndVel, but for fictracSmo values
%       pkTimeRelOpto - same as INPUT
%       cond - same as INPUT
%       fwdVelCond - same as INPUT
%       maxNumSteps - same as INPUT
%
% CREATED: 4/17/24 - HHY
%
% UPDATED:
%   4/17/24 - HHY
%
function extractOptomotorBoutsLegStepParams_fly(cond, fwdVelCond, ...
    pkTimeRelOpto, maxNumSteps, pDataFNames, pDataPath, saveFilePath, ...
    saveFileSuffix)

    % names of all step parameters to save
    stepParamNames = {'stepLengths', 'stepXLengths',...
        'stepYLengths', 'stepDirections', 'stepDurations', 'stepSpeeds',...
        'stepVelX', 'stepVelY', 'stepAEPX', 'stepAEPY', 'stepPEPX', ...
        'stepPEPY', 'stepFtFwd', 'stepFtLat', 'stepFtYaw'};

    % prompt user to select pData files
    if isempty(pDataFNames)
        [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
            'Select pData files for 1 fly', ...
            pDataPath, 'MultiSelect', 'on');
    else
        pDataDirPath = pDataPath;
    end
    
    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end

    % preallocate outputs
    for i = 1:length(stepParamNames)
        selLegSteps.(stepParamNames{i}) = [];
        selStanceParams.(stepParamNames{i}) = [];
        selSwingParams.(stepParamNames{i}) = [];
    end

    boutPeakVelSmo.yaw = [];
    boutPeakVelSmo.fwd = [];
    boutPeakVelSmo.slide = [];

    boutStartVelSmo.yaw = [];
    boutStartVelSmo.fwd = [];
    boutStartVelSmo.slide = [];

    boutEndVelSmo.yaw = [];
    boutEndVelSmo.fwd = [];
    boutEndVelSmo.slide = [];

    boutPeakVel.yaw = [];
    boutPeakVel.fwd = [];
    boutPeakVel.slide = [];

    boutStartVel.yaw = [];
    boutStartVel.fwd = [];
    boutStartVel.slide = [];

    boutEndVel.yaw = [];
    boutEndVel.fwd = [];
    boutEndVel.slide = [];

    pkSwingStance = [];

    pkTimeOpto = [];

    selNDs = [];

    selOptoVels = [];

    countNumBouts = 0; % counter for number of bouts

    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end

        pDataFullPath = [pDataDirPath filesep pDataName];

        % save fly name as first pDataName's date, fly, cell (19 characters)
        if (i == 1)
            flyName = pDataName(1:19);
        end

        % get variables saved in pData file
        pDatVars = whos('-file', pDataFullPath);
    
        pDatVarsNames = cell(size(pDatVars));
        
        % convert pDatVars into cell array of just names
        for j = 1:length(pDatVars)
            pDatVarsNames{j} = pDatVars(j).name;
        end

        % check if this pData file has legSteps, fictracSmo, moveNotMove
        %  opto, visstim structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legSteps')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracSmo')) || ...
                ~any(strcmpi(pDatVarsNames, 'opto')) || ...
                ~any(strcmpi(pDatVarsNames, 'visstim')) || ...
                ~any(strcmpi(pDatVarsNames, 'moveNotMove')))
            continue;
        end

        % load variables from pData
        load(pDataFullPath, 'legTrack', 'moveNotMove', 'fictracSmo', ...
                'legSteps', 'stanceStepParams', 'swingStepParams', ...
                'opto', 'fictracProc', 'visstim')

        % convert legTrack.refPts to legIDs
        legIDs.ind = legSteps.legIDs.ind;
        legIDs.name = legSteps.legIDs.names;



        % get yaw velocity peaks for right turns
        [rightPeakTimes, rightStartTimes, rightEndTimes, ...
            rightPeakInd, rightStartInd, rightEndInd, rightOptoVels, ...
            rightTimesRelOpto] = findOptomotorCondYawVelPeaksFT(...
            fictracSmo, legSteps, cond, fwdVelCond, moveNotMove, ...
            visstim, pkTimeRelOpto, true);


        % get yaw velocity peaks for left turns
        [leftPeakTimes, leftStartTimes, leftEndTimes, ...
            leftPeakInd, leftStartInd, leftEndInd, leftOptoVels, ...
            leftTimesRelOpto] = findOptomotorCondYawVelPeaksFT(...
            fictracSmo, legSteps, cond, fwdVelCond, moveNotMove, ...
            visstim, pkTimeRelOpto, false);
        

        % check if this pData file contributes any turns
        % if not, skip and move to next pData file
        if (isempty(rightPeakInd) && isempty(leftPeakInd))
            continue;
        end

        % number of bouts for this trial
        thisNumBouts = length(rightPeakInd) + length(leftPeakInd);

        % update counter
        countNumBouts = countNumBouts + thisNumBouts;


        % get indices for steps aligned to bouts
        % right turns
        [rightStepInd, rightPkSwingStance, rightRmInd] = ...
            findBoutStepIndsBall_full(legSteps, rightPeakTimes, ...
            legTrack.t, maxNumSteps, legIDs);

        % left turns
        [leftStepInd, leftPkSwingStance, leftRmInd] = ...
            findBoutStepIndsBall_full(legSteps, leftPeakTimes, ...
            legTrack.t, maxNumSteps, legIDs);



        % remove bouts
        rightPeakInd(rightRmInd) = [];
        leftPeakInd(leftRmInd) = [];
        rightStartInd(rightRmInd) = [];
        leftStartInd(leftRmInd) = [];
        rightEndInd(rightRmInd) = [];
        leftEndInd(leftRmInd) = [];
        rightOptoVels(rightRmInd) = [];
        leftOptoVels(leftRmInd) = [];
        rightTimesRelOpto(rightRmInd) = [];
        leftTimesRelOpto(leftRmInd) = [];
        rightPeakTimes(rightRmInd) = [];
        leftPeakTimes(leftRmInd) = [];

        % get FicTrac and FicTracSmo values at bout peaks
        % get FicTracSmo values at bout peaks
        boutPeakVelSmo.yaw = [boutPeakVelSmo.yaw; ...
            fictracSmo.yawAngVel(rightPeakInd)];
        boutPeakVelSmo.yaw = [boutPeakVelSmo.yaw; ...
            fictracSmo.yawAngVel(leftPeakInd)];
        boutPeakVelSmo.fwd = [boutPeakVelSmo.fwd; ...
            fictracSmo.fwdVel(rightPeakInd)];
        boutPeakVelSmo.fwd = [boutPeakVelSmo.fwd; ...
            fictracSmo.fwdVel(leftPeakInd)];
        boutPeakVelSmo.slide = [boutPeakVelSmo.slide; ...
            fictracSmo.slideVel(rightPeakInd)];
        boutPeakVelSmo.slide = [boutPeakVelSmo.slide; ...
            fictracSmo.slideVel(leftPeakInd)];

        % get FicTracSmo values at bout starts
        boutStartVelSmo.yaw = [boutStartVelSmo.yaw; ...
            fictracSmo.yawAngVel(rightStartInd)];
        boutStartVelSmo.yaw = [boutStartVelSmo.yaw; ...
            fictracSmo.yawAngVel(leftStartInd)];
        boutStartVelSmo.fwd = [boutStartVelSmo.fwd; ...
            fictracSmo.fwdVel(rightStartInd)];
        boutStartVelSmo.fwd = [boutStartVelSmo.fwd; ...
            fictracSmo.fwdVel(leftStartInd)];
        boutStartVelSmo.slide = [boutStartVelSmo.slide; ...
            fictracSmo.slideVel(rightStartInd)];
        boutStartVelSmo.slide = [boutStartVelSmo.slide; ...
            fictracSmo.slideVel(leftStartInd)];

        % get FicTracSmo values at bout ends
        boutEndVelSmo.yaw = [boutEndVelSmo.yaw; ...
            fictracSmo.yawAngVel(rightEndInd)];
        boutEndVelSmo.yaw = [boutEndVelSmo.yaw; ...
            fictracSmo.yawAngVel(leftEndInd)];
        boutEndVelSmo.fwd = [boutEndVelSmo.fwd; ...
            fictracSmo.fwdVel(rightEndInd)];
        boutEndVelSmo.fwd = [boutEndVelSmo.fwd; ...
            fictracSmo.fwdVel(leftEndInd)];
        boutEndVelSmo.slide = [boutEndVelSmo.slide; ...
            fictracSmo.slideVel(rightEndInd)];
        boutEndVelSmo.slide = [boutEndVelSmo.slide; ...
            fictracSmo.slideVel(leftEndInd)];

        % get FicTrac values at bout peaks
        boutPeakVel.yaw = [boutPeakVel.yaw; ...
            fictracProc.yawAngVel(rightPeakInd)];
        boutPeakVel.yaw = [boutPeakVel.yaw; ...
            fictracProc.yawAngVel(leftPeakInd)];
        boutPeakVel.fwd = [boutPeakVel.fwd; ...
            fictracProc.fwdVel(rightPeakInd)];
        boutPeakVel.fwd = [boutPeakVel.fwd; ...
            fictracProc.fwdVel(leftPeakInd)];
        boutPeakVel.slide = [boutPeakVel.slide; ...
            fictracProc.slideVel(rightPeakInd)];
        boutPeakVel.slide = [boutPeakVel.slide; ...
            fictracProc.slideVel(leftPeakInd)];

        % get FicTrac values at bout starts
        boutStartVel.yaw = [boutStartVel.yaw; ...
            fictracProc.yawAngVel(rightStartInd)];
        boutStartVel.yaw = [boutStartVel.yaw; ...
            fictracProc.yawAngVel(leftStartInd)];
        boutStartVel.fwd = [boutStartVel.fwd; ...
            fictracProc.fwdVel(rightStartInd)];
        boutStartVel.fwd = [boutStartVel.fwd; ...
            fictracProc.fwdVel(leftStartInd)];
        boutStartVel.slide = [boutStartVel.slide; ...
            fictracProc.slideVel(rightStartInd)];
        boutStartVel.slide = [boutStartVel.slide; ...
            fictracProc.slideVel(leftStartInd)];

        % get FicTrac values at bout ends
        boutEndVel.yaw = [boutEndVel.yaw; ...
            fictracProc.yawAngVel(rightEndInd)];
        boutEndVel.yaw = [boutEndVel.yaw; ...
            fictracProc.yawAngVel(leftEndInd)];
        boutEndVel.fwd = [boutEndVel.fwd; ...
            fictracProc.fwdVel(rightEndInd)];
        boutEndVel.fwd = [boutEndVel.fwd; ...
            fictracProc.fwdVel(leftEndInd)];
        boutEndVel.slide = [boutEndVel.slide; ...
            fictracProc.slideVel(rightEndInd)];
        boutEndVel.slide = [boutEndVel.slide; ...
            fictracProc.slideVel(leftEndInd)];



        % peak swing/stance call, merge across legs
        thisPkSwingStance = cat(2, rightPkSwingStance, leftPkSwingStance);
        % add to output matrix
        pkSwingStance = cat(2, pkSwingStance, thisPkSwingStance);


        % peak time relative to opto, add to output matrix
        pkTimeOpto = [pkTimeOpto; rightTimesRelOpto; leftTimesRelOpto];

        % optomotor velocity, add to output matrix
        selOptoVels = [selOptoVels; rightOptoVels; leftOptoVels];


        % get NDs for each peak
        rightNDs = zeros(size(rightPeakTimes));
        for j = 1:length(rightPeakTimes)
            startOptoInd = find(rightPeakTimes(j) >= ...
                opto.stimStartTimes,1,'last');
            endOptoInd = find(rightPeakTimes(j) <= ...
                opto.stimEndTimes,1,'first');

            % if peak falls during opto stim
            if (startOptoInd == endOptoInd)
                rightNDs(j) = opto.stimParams.ndFilter;
            else % no opto
                rightNDs(j) = -1;
            end
        end

        leftNDs = zeros(size(leftPeakTimes));
        for j = 1:length(leftPeakTimes)
            startOptoInd = find(leftPeakTimes(j) >= ...
                opto.stimStartTimes,1,'last');
            endOptoInd = find(leftPeakTimes(j) <= ...
                opto.stimEndTimes,1,'first');

            % if peak falls during opto stim
            if (startOptoInd == endOptoInd)
                leftNDs(j) = opto.stimParams.ndFilter;
            else % no opto
                leftNDs(j) = -1;
            end
        end

        % add NDs to output
        selNDs = [selNDs; rightNDs; leftNDs];



        % get legStep parameters from aligned step indices

        % preallocate
        oneParamValsRight = nan(size(rightStepInd,1), ...
            size(rightStepInd,2), 2, size(rightStepInd,3));
        oneParamStValsRight = nan(size(rightStepInd));
        oneParamSwValsRight = nan(size(rightStepInd));

        oneParamValsLeft = nan(size(leftStepInd,1), ...
            size(leftStepInd,2), 2, size(leftStepInd,3));
        oneParamStValsLeft = nan(size(leftStepInd));
        oneParamSwValsLeft = nan(size(leftStepInd));

        for j = 1:length(stepParamNames)
            thisParam = legSteps.(stepParamNames{j});
            thisStParam = stanceStepParams.(stepParamNames{j});
            thisSwParam = swingStepParams.(stepParamNames{j});
            
            % for right turns
            % loop over all steps of bout
            for k = 1:size(rightStepInd, 1)
                % loop over all legs
                for l = 1:size(rightStepInd, 2)
                    % loop over all bouts
                    for m = 1:size(rightStepInd, 3)
                        thisInd = rightStepInd(k,l,m);
                        % if index is NaN, don't do anything (will keep
                        %  NaN), otherwise, put in appropriate value
                        if ~isnan(thisInd)
                            % legSteps
                            thisVals = thisParam(thisInd,:);
                            oneParamValsRight(k,l,1,m) = thisVals(1);
                            oneParamValsRight(k,l,2,m) = thisVals(2);

                            % stance only
                            oneParamStValsRight(k,l,m) = ...
                                thisStParam(thisInd);
                            % swing only
                            oneParamSwValsRight(k,l,m) = ...
                                thisSwParam(thisInd);
                        end
                    end
                end
            end

            % for left turns
            % loop over all steps of bout
            for k = 1:size(leftStepInd, 1)
                % loop over all legs
                for l = 1:size(leftStepInd, 2)
                    % loop over all bouts
                    for m = 1:size(leftStepInd, 3)
                        thisInd = leftStepInd(k,l,m);
                        % if index is NaN, don't do anything (will keep
                        %  NaN), otherwise, put in appropriate value
                        if ~isnan(thisInd)
                            % legSteps
                            thisVals = thisParam(thisInd,:);
                            oneParamValsLeft(k,l,1,m) = thisVals(1);
                            oneParamValsLeft(k,l,2,m) = thisVals(2);

                            % stance only
                            oneParamStValsLeft(k,l,m) = ...
                                thisStParam(thisInd);
                            % swing only
                            oneParamSwValsLeft(k,l,m) = ...
                                thisSwParam(thisInd);
                        end
                    end
                end
            end


            % concatenate right and left
            oneParamVals = cat(4, oneParamValsRight, oneParamValsLeft);
            oneParamStVals = cat(3, oneParamStValsRight, oneParamStValsLeft);
            oneParamSwVals = cat(3, oneParamSwValsRight, oneParamSwValsLeft);

            % add to output matrices
            selLegSteps.(stepParamNames{j}) = cat(4, ...
                selLegSteps.(stepParamNames{j}), oneParamVals);

            selStanceParams.(stepParamNames{j}) = cat(3, ...
                selStanceParams.(stepParamNames{j}), oneParamStVals);
            selSwingParams.(stepParamNames{j}) = cat(3, ...
                selSwingParams.(stepParamNames{j}), oneParamSwVals);
        end
    end

    numBouts = countNumBouts;

    % save output file
    fullSavePath = [saveFilePath filesep flyName saveFileSuffix '.mat'];

    save(fullSavePath, 'selLegSteps', 'selStanceParams', ...
        'selSwingParams', 'selNDs', 'selOptoVels', 'pkTimeOpto', ...
        'pkSwingStance', 'numBouts', 'boutPeakVel', 'boutPeakVelSmo', ...
        'boutStartVel', 'boutStartVelSmo', 'boutEndVel', 'boutEndVelSmo', ...
        'pkTimeRelOpto', 'cond', 'fwdVelCond', 'maxNumSteps', '-v7.3');
end