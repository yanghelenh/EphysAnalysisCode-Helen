% saveSpikerate_bouts.m
%
% Function that saves spike rate for turning bouts, defined by yaw
%  velocity peaks. User can specify FicTrac conditions that yaw
%  velocity peaks must meet.
% Extracts bouts that meet FicTrac conditions and accounts for timing of 
%  step relative to yaw velocity peak
% Extracts both right and left turns (flips left turns over midline). Any
%  conditions on yaw or side velocity are considered relative to right
%  turns, and computed for both. Can remove bias in turning before computing
%  left and right together
% Removes peaks that fall during periods of stimulation (opto or I inj),
%  during periods of not moving, or when FicTrac lost
% User selects one or more pData files through GUI
%
% Adaptation of saveLegStepParamCond_bouts(), which is for free-walking
%  flies
% 
% INPUTS:
%   delay - time offset between behavior and ephys, neg for ephys before
%       behavior
%   minYawThresh - minimum yaw velocity to define start and end of bout
%   turnDur - 2 element vector [minTurnDuration maxTurnDuration] to
%       specify the min and max duration of the turning bout for it to
%       be included
%   spikerateParams - struct of parameters on spike rate
%       maxDuration - time in seconds to consider on each side 
%       interpFrameRate - frame rate to interpolate to, in Hz
%   postStimExclDur - additional time, in sec, after stimulation (opto or 
%       iInj to exclude turning bouts from overlapping with
%   pDataPath - full path to pData directory
%   saveFilePath - directory in which to save output file
%   saveFileName - name of output file, without .mat part
%
% OUTPUTS:
%   none, but saves output file with name saveFileName in saveFilePath
%       allSpikerate - matrix of numTimePts x numLegs x numBouts for 
%         aligned spike rate
%       allYaw - matrix of numTimePts x numLegs x numBouts for 
%         aligned yaw velocity
%       allFwd - matrix of numTimePts x numLegs x numBouts for 
%         aligned fwd velocity
%       numBouts - total number of bouts (max n for stepParam (if no NaNs))
%       pDataFiles - struct of info on pData files
%           names - name of each pData file with at least 1 valid step, as
%               cell array
%           inds - indices (corresponding to bout indices) that
%               belong to each pData file, as cell array of vectors
%     FROM INPUT: delay, minYawThresh, turnDur, spikerateParams, 
%       postStimExclDur
%
% CREATED: 7/17/23 - HHY
%
% UPDATED:
%   7/17/23 - HHY
%   8/2/23 - HHY - update comments
%
function saveSpikerate_bouts(delay, minYawThresh, turnDur, ...
    spikerateParams, postStimExclDur, pDataPath, saveFilePath, ...
    saveFileName)

    % prompt user to select pData files
    [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
        'Select pData files for 1 fly', pDataPath, 'MultiSelect', 'on');
    
    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end

    % preallocate outputs
    pDataFiles.names = pDataFNames;
    pDataFiles.inds = [];

    allSpikerate = [];
    allYaw = [];
    allFwd = [];


    rmvInd = []; % indices of pData files to remove
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

        % get variables saved in pData file
        pDatVars = whos('-file', pDataFullPath);
    
        pDatVarsNames = cell(size(pDatVars));
        
        % convert pDatVars into cell array of just names
        for j = 1:length(pDatVars)
            pDatVarsNames{j} = pDatVars(j).name;
        end

        % check if this pData file has legSteps, fictracSmo, moveNotMove
        %  structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'ephysSpikes')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracSmo')) || ...
                ~any(strcmpi(pDatVarsNames, 'moveNotMove')))
            rmvInd = [rmvInd; i];
            continue;
        end

        % load variables from pData
        % check if opto or iInj trial, load those as appropriate
        % also, save opto or iInj times as stim times; no stim, empty
        if (any(strcmpi(pDatVarsNames, 'iInj'))) % current inj
            load(pDataFullPath, 'ephysSpikes', 'moveNotMove', 'fictracSmo', ...
                 'fictracProc', 'iInj');

            stimStartTimes = iInj.startTimes;
            stimEndTimes = iInj.endTimes;
        else % no stimulation
            load(pDataFullPath, 'ephysSpikes', 'moveNotMove', ...
                'fictracProc', 'fictracSmo');

            stimStartTimes = [];
            stimEndTimes = [];
        end

        % get yaw velocity peaks for right turns
        [rightPeakTimes, rightStartTimes, rightEndTimes, ...
            rightPeakInd, rightStartInd, rightEndInd] = ...
            findYawVelPeaksFT(fictracSmo, minYawThresh, turnDur, ...
            moveNotMove, true);

        % get yaw velocity peaks for left turns
        [leftPeakTimes, leftStartTimes, leftEndTimes, ...
            leftPeakInd, leftStartInd, leftEndInd] = ...
            findYawVelPeaksFT(fictracSmo, minYawThresh, turnDur, ...
            moveNotMove, false);

        % remove turning bouts that overlap with stimulation bouts
        if (~isempty(stimStartTimes)) % only run on those with stimulation
            % right turns
            % get bouts that overlap
            boutInd2RmvRight = findBoutsDurStim2Rmv(rightStartTimes, ...
                rightEndTimes, stimStartTimes, stimEndTimes, ...
                postStimExclDur);
            % remove overlapping bouts from all bout times/inds
            rightPeakTimes(boutInd2RmvRight) = [];
            rightStartTimes(boutInd2RmvRight) = [];
            rightEndTimes(boutInd2RmvRight) = [];
            rightPeakInd(boutInd2RmvRight) = [];
            rightStartInd(boutInd2RmvRight) = [];
            rightEndInd(boutInd2RmvRight) = [];

            % left turns
            % get bouts that overlap
            boutInd2RmvLeft = findBoutsDurStim2Rmv(leftStartTimes, ...
                leftEndTimes, stimStartTimes, stimEndTimes, ...
                postStimExclDur);
            % remove overlapping bouts from all bout times/inds
            leftPeakTimes(boutInd2RmvLeft) = [];
            leftStartTimes(boutInd2RmvLeft) = [];
            leftEndTimes(boutInd2RmvLeft) = [];
            leftPeakInd(boutInd2RmvLeft) = [];
            leftStartInd(boutInd2RmvLeft) = [];
            leftEndInd(boutInd2RmvLeft) = [];
        end
        

        % check if this pData file contributes any turns
        % if not, skip and move to next pData file
        if (isempty(rightPeakInd) && isempty(leftPeakInd))
            rmvInd = [rmvInd;i];
            continue;
        end

        % number of bouts for this trial
        thisNumBouts = length(rightPeakInd) + length(leftPeakInd);

        % get indices for bouts for this trial
        thisTrialStartInd = 1 + countNumBouts;
        thisTrialEndInd = thisTrialStartInd + thisNumBouts - 1;

        thisTrialInds = [thisTrialStartInd thisTrialEndInd];
    
        pDataFiles.inds = [pDataFiles.inds; thisTrialInds];

        % update counter
        countNumBouts = countNumBouts + thisNumBouts;


        % get spike rate, yaw vel, fwd velocity for each bout
        % use fictracProc

        % spike rate, right
        [spikerateRight, tBins] = getSpikerateFromBouts(ephysSpikes, ...
            spikerateParams, rightPeakTimes, rightStartTimes, ...
            rightEndTimes, delay);
        % spike rate, left
        [spikerateLeft, ~] = getSpikerateFromBouts(ephysSpikes, ...
            spikerateParams, leftPeakTimes, leftStartTimes, ...
            leftEndTimes, delay);

        % yaw vel, right
        [yawVelRight, ~] = getFictracFromBouts(fictracProc, 'yawAngVel', ...
            spikerateParams, rightPeakInd, rightStartInd, rightEndInd);
        % yaw vel, left
        [yawVelLeft, ~] = getFictracFromBouts(fictracProc, 'yawAngVel', ...
            spikerateParams, leftPeakInd, leftStartInd, leftEndInd);

        % fwd vel, right
        [fwdVelRight, ~] = getFictracFromBouts(fictracProc, 'fwdVel', ...
            spikerateParams, rightPeakInd, rightStartInd, rightEndInd);
        % fwd vel, left
        [fwdVelLeft, ~] = getFictracFromBouts(fictracProc, 'fwdVel', ...
            spikerateParams, leftPeakInd, leftStartInd, leftEndInd);

        % concatenate all
        allSpikerate = cat(2, allSpikerate, spikerateRight);
        allSpikerate = cat(2, allSpikerate, spikerateLeft);

        allYaw = cat(2, allYaw, yawVelRight);
        allYaw = cat(2, allYaw, yawVelLeft);

        allFwd = cat(2, allFwd, fwdVelRight);
        allFwd = cat(2, allFwd, fwdVelLeft);
    end

    % remove unused pData files
    pDataFiles.names(rmvInd) = [];

    numBouts = countNumBouts;

    % save output file
    fullSavePath = [saveFilePath filesep saveFileName '.mat'];

    save(fullSavePath, 'allSpikerate', 'allYaw', 'allFwd',...
         'numBouts', 'pDataFiles', 'delay', 'minYawThresh', 'turnDur', ...
         'spikerateParams', 'postStimExclDur','-v7.3');
end