% findOptomotorCondYawVelPeaksFT.m
%
% Helper function for extractOptomotorBoutsLegStepParams_fly() that takes 
%  in fictracSmo, cond, fwdVelCond, moveNotMove, visstim, boutTimeWin,
%  pkTimeRelOpto, and 
%  returns times for valid yaw velocity peaks as well as start/end times
%  for corresponding bouts.
% Boolean for whether to extract right or left turns
% 
% Adapted from findCondYawVelPeaksFT(), which doesn't incorporate vis
% Use times instead of indices b/c legStep indices and FicTrac indices
%  aren't the same, unlike in free-walking data. Also, indices just in case
%
% INPUTS:
%   fictracSmo - pData output struct, from computeSmoFictrac()
%   legSteps - pData output struct
%   cond - pass from extractOptomotorBoutsLegStepParams_fly()
%   fwdVelCond - pass from extractOptomotorBoutsLegStepParams_fly()
%   moveNotMove - pData output struct
%   visstim - pData output struct, of visual stimulus info
%   pkTimeRelOpto - pass from extractOptomotorBoutsLegStepParams_fly()
%   rightTurn - boolean for whether to extract right turns (false = left)
%
% OUTPUTS:
%   yawVelPeakTimes - vector of times corresponding to yaw
%       velocity peaks (that meet cond criteria), length n
%   boutStartTimes - vector of times corresponding to start of bouts,
%       length n
%   boutEndTimes - vector of times corresponding to end of bouts, length n
%   yawVelPeakInd - vector of indices of fictracSmo/Proc corresponding to
%       yaw velocity peaks (that meet cond criteria), length n
%   boutStartInd - vector of indices corresponding to start of bouts,
%       length n
%   boutEndInd - vector of indices corresponding to end of bouts, length n
%   peakOptoVels - vector of optomotor velocities for each peak, length n
%   peakTimesRelOpto - vector of times of peaks, relative to optomotor
%       stimlus, length n
%   
%
% CREATED: 4/17/24 - HHY
%
% UPDATED:
%   4/17/24 - HHY
%
function [yawVelPeakTimes, boutStartTimes, boutEndTimes, ...
    yawVelPeakInd, boutStartInd, boutEndInd, peakOptoVels, ...
    peakTimesRelOpto] = findOptomotorCondYawVelPeaksFT(fictracSmo, ...
    legSteps, cond, fwdVelCond, moveNotMove, visstim, pkTimeRelOpto, ...
    rightTurn)

    % check if we're extracting right or left turns
    if (rightTurn)
        angVelSmoS = fictracSmo.yawAngVel;
    else
        angVelSmoS = -1 * fictracSmo.yawAngVel;
    end
        
    % find all yaw velocity peaks, no conditioning yet
    % always operates on fictracSmo.yawAngVel
    [~, pkInds] = findpeaks(angVelSmoS);

    % remove peaks during not moving times
    pkInds = setdiff(pkInds, moveNotMove.ftNotMoveInd,'stable');

    % remove peaks during dropInd
    pkInds = setdiff(pkInds, fictracSmo.dropInd, 'stable');

    % convert pkInds to logical (for fictrac indices)
    pkLog = false(size(fictracSmo.t));
    pkLog(pkInds) = true;

    % find peaks that meet conditions of cond
    for i = 1:length(cond.whichParam)

        % condition on step walking forward/backwards
        if (strcmpi(cond.whichParam{i}, 'stepFwdBool'))
            thisTime = fictracSmo.t;
            thisFwdLog = getStepFwdLogical(legSteps, thisTime,...
                cond.legs{i});
            if ~(eval(cond.cond{i})) % if target is false, invert
                thisCondLog = ~thisFwdLog;
            else
                thisCondLog = thisFwdLog;
            end
        else
            % if condition is on yaw or on lateral, invert if left turns
            % don't touch if right turn
            if (rightTurn)
                thisCond = fictracSmo.(cond.whichParam{i});
            else % left turn
                if (contains(cond.whichParam{i}, 'angVel', 'IgnoreCase',true) || ...
                        contains(cond.whichParam{i}, 'slideVel', 'IgnoreCase',true))
                    thisCond = -1 * fictracSmo.(cond.whichParam{i});
                else
                    thisCond = fictracSmo.(cond.whichParam{i});
                end
            end
             
            % all time points that meet condition
            thisCondLog = eval(['thisCond' cond.cond{i}]);
        end

        % intersect all points that meet condition with peaks
        pkLog = pkLog & thisCondLog;
    end

    % convert logical back to indices
    pkInds = find(pkLog);


    % check that peaks occur within the valid time window relative to the
    %  optomotor stimulus
    pkTimes = fictracSmo.t(pkInds);

    % get start and end times of valid time windows relative to optomotor
    %  stimulus
    % only turns in correct direction
    if (rightTurn)
        valOptoTrials = visstim.rampVels > 0;
    else
        valOptoTrials = visstim.rampVels < 0;
    end
    valOptoStartTimes = visstim.rampStartTimes(valOptoTrials) + ...
        pkTimeRelOpto(1);
    valOptoEndTimes = visstim.rampEndTimes(valOptoTrials) + ...
        pkTimeRelOpto(2);

    rmvInd = [];
    % optomotor trial index for the peak
    pkOptoInd = zeros(size(pkInds));
    % for each peak, check if it falls within valid time
    for i = 1:length(pkTimes)
        % index of valid optomotor start and end times, relative to peak
        thisStartOpto = find(pkTimes(i) >= valOptoStartTimes, 1, 'last');
        thisEndOpto = find(pkTimes(i) <= valOptoEndTimes, 1, 'first');

        % if peak falls within window, start and end indices are the same,
        %  remove otherwise
        if (isempty(thisStartOpto))
            rmvInd = [rmvInd i];
        elseif (isempty(thisEndOpto))
            rmvInd = [rmvInd i];
        elseif (thisStartOpto ~= thisEndOpto)
            rmvInd = [rmvInd i];
        else
            pkOptoInd(i) = thisStartOpto;
        end
    end

    % remove peaks that don't fall within optomotor window
    pkInds(rmvInd) = [];
    pkOptoInd(rmvInd) = [];

    % all indices where yaw velocity is less than min
    yawMinInd = find(angVelSmoS < cond.minYawThresh);

    % preallocate start and end index vectors
    pkStartInd = zeros(size(pkInds));
    pkEndInd = zeros(size(pkInds));

    % find start and end for each peak
    % if start and end can't be found (too close to edge, remove peaks)
    rmvInd = [];
    for i = 1:length(pkInds)
        thisStart = find(yawMinInd < pkInds(i), 1, 'last');
        % if no index found, set peak start to beginning of trial
        if isempty(thisStart)
            rmvInd = [rmvInd i];
        else
            % +1 because presumably next index is greater than min
            pkStartInd(i) = yawMinInd(thisStart) + 1;
        end

        thisEnd = find(yawMinInd > pkInds(i), 1, 'first');
        % if no index found, set peak end to end of trial
        if isempty(thisEnd)
            rmvInd = [rmvInd i];
        else
            % -1 because this is first that is less than min
            pkEndInd(i) = yawMinInd(thisEnd) - 1;
        end
    end
    % remove peaks too close to edge
    pkInds(rmvInd) = [];
    pkStartInd(rmvInd) = [];
    pkEndInd(rmvInd) = [];
    pkOptoInd(rmvInd) = [];

    initSizeStart = length(pkStartInd);

    % remove any peaks whose edges fall within dropInd
    [pkStartInd, startNotRmvInd] = setdiff(pkStartInd, fictracSmo.dropInd, ...
        'stable');
    startRmvInd = setxor(1:initSizeStart, startNotRmvInd);
    pkInds(startRmvInd) = [];
    pkEndInd(startRmvInd) = [];
    pkOptoInd(startRmvInd) = [];

    initSizeEnd = length(pkEndInd);
    [pkEndInd, endNotRmvInd] = setdiff(pkEndInd, fictracSmo.dropInd, ...
        'stable');
    endRmvInd = setxor(1:initSizeEnd, endNotRmvInd);
    pkInds(endRmvInd) = [];
    pkStartInd(endRmvInd) = [];
    pkOptoInd(endRmvInd) = [];




    % check if any of the peaks share edges - look only at peak starts
    %  (will be the same as peak ends)

    % get indices into bodytraj, for not unique starts 
    [~, uniInd]  = unique(pkStartInd, 'stable');
    nonUniInd = setdiff(1:numel(pkStartInd), uniInd);
    startNonUniInd = pkStartInd(nonUniInd);

    % indices of peaks themselves, for non-unique peaks
    peakNonUniInd = pkInds(nonUniInd); 


    % if peaks share edges, keep only peak with greatest yaw velocity
    if (~isempty(startNonUniInd))
        notKeepInd = []; % initialize tracker of peaks to discard
        for i = 1:length(startNonUniInd)
            thisPeakStart = startNonUniInd(i);

            % find all peaks that share this start
            sharedInd = find(pkStartInd == thisPeakStart);

            % find peak with greatest yaw velocity
            bInd = pkInds(sharedInd);
            yawVals = angVelSmoS(bInd);

            [~, keepPkInd] = max(yawVals);
            keepPkBInd = bInd(keepPkInd);
            notKeepPkBInd = setdiff(bInd, keepPkBInd, 'stable');

            if (iscolumn(notKeepPkBInd))
                notKeepPkBInd = notKeepPkBInd';
            end

            % track all peaks to discard
            notKeepInd = [notKeepInd notKeepPkBInd];
        end

        % discard not keep peaks
        % discard for peaks, get indices of ones to keep
        [pkInds, keepInd] = setdiff(pkInds,notKeepInd,'stable');
        pkStartInd = pkStartInd(keepInd);
        pkEndInd = pkEndInd(keepInd);
        pkInds = pkInds(keepInd);
        pkOptoInd = pkOptoInd(keepInd);
    end

    % loop through all peaks and check that the bout meets the duration
    %  requirements
    % initialize to keep track of peaks to remove
    rmInd = []; 

    for i = 1:length(pkInds)

        thisBoutStartT = fictracSmo.t(pkStartInd(i));
        thisBoutEndT = fictracSmo.t(pkEndInd(i));

        thisBoutDur = thisBoutEndT - thisBoutStartT;

        % if this bout duration is too short or too long, flag this index
        %  for deletion
        if((thisBoutDur < cond.turnDur(1)) || ...
                (thisBoutDur > cond.turnDur(2)))
            rmInd = [rmInd i];
        end
    end
    % remove any bouts that don't meet the duration criteria
    pkInds(rmInd) = [];
    pkStartInd(rmInd) = [];
    pkEndInd(rmInd) = [];
    pkOptoInd(rmInd) = [];

    % loop through all peaks and check that the bout meets the forward
    %  velocity requirements
    % initialize to keep track of peaks to remove
    rmInd = [];
    for i = 1:length(pkInds)
        % forward velocity at start for this bout
        thisBoutInitFwdVel = fictracSmo.fwdVel(pkStartInd(i));
        % forward velocity at yaw peak for this bout
        thisBoutPeakFwdVel = fictracSmo.fwdVel(pkInds(i));

        % change in forward velocity
        thisBoutChangeFwdVel = thisBoutPeakFwdVel - thisBoutInitFwdVel;

        % if this bout doesn't meet forward velocity requirements, flag
        %  this index for deletion
        if ((thisBoutInitFwdVel < fwdVelCond.initVel(1)) || ...
                (thisBoutInitFwdVel > fwdVelCond.initVel(2)) || ...
                (thisBoutChangeFwdVel < fwdVelCond.change(1)) || ...
                (thisBoutChangeFwdVel > fwdVelCond.change(2)))
            rmInd = [rmInd i];
        end
    end
    % remove any bouts that don't meet the forward velocity criteria
    pkInds(rmInd) = [];
    pkStartInd(rmInd) = [];
    pkEndInd(rmInd) = [];
    pkOptoInd(rmInd) = [];

    % outputs, indices
    yawVelPeakInd = pkInds;
    boutStartInd = pkStartInd;
    boutEndInd = pkEndInd;

    % outputs, times
    yawVelPeakTimes = fictracSmo.t(pkInds);
    boutStartTimes = fictracSmo.t(pkInds);
    boutEndTimes = fictracSmo.t(pkInds);

    % outputs, optomotor velocity
    peakOptoVels = zeros(size(yawVelPeakTimes));
    peakTimesRelOpto = zeros(size(yawVelPeakTimes));

    valCmdVel = visstim.rampCmdVels(valOptoTrials);
    for i = 1:length(yawVelPeakTimes)
        % index of optomotor trial for this peak
        thisOptoInd = find(yawVelPeakTimes(i) >= valOptoStartTimes, ...
            1, 'last');
        % velocity of this optomotor stimulus
        peakOptoVels(i) = valCmdVel(thisOptoInd);

        % time of yaw velocity peak, relative to optomotor stimulus start
        peakTimesRelOpto(i) = yawVelPeakTimes(i) - ...
            valOptoStartTimes(thisOptoInd) + pkTimeRelOpto(1);
    end
end