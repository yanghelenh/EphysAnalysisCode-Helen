% saveSpikerateCond_bouts.m
%
% Function that saves spike rate over time for turning bouts, defined by
%  yaw velocity peaks. User can specify FicTrac conditions that yaw
%  velocity peaks must meet as well as the initial forward velocity range 
%  and the change in forward velocity between the start and peak of the
%  turn
% User selects one or more pData files through GUI, for 1 fly
%
% Adaptation of saveBallLegStepParamCond_bouts(), which analyzes the leg
%  step parameters for ball-walking behavior
% 
% INPUTS:
%   cond - struct of conditions for yaw velocity peak, if multiple 
%     conditions, treats it as AND
%       whichParam - cell array (even if 1 element) on which fictracSmo
%           field to condition on, one for each condition
%       cond - cell array of strings to condition on, for eval(); same size
%           as whichParam
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
%   spikerateParams - struct of parameters on spike rate
%       maxDuration - time in seconds to consider on each side 
%       interpFrameRate - frame rate to interpolate to, in Hz
%   flipLR - boolean for whether to flip left/right (inverts yaw and
%       lateral velocity)
%   postStimExclDur - additional time, in sec, after stimulation (iInj) 
%       to exclude turning bouts from overlapping with
%   pDataFNames - cell array of pData file names or [] if select through
%       GUI
%   pDataPath - full path to pData directory
%   saveFilePath - directory in which to save output file
%   saveFileName - name of output file, without .mat part
%
% OUTPUTS:
%   none, but saves output file with name saveFileName in saveFilePath
%       allSpikerate - matrix of numTimePts x numBouts for aligned spike
%         rate
%       allYaw - matrix of numTimePts x numBouts for aligned yaw velocity
%       allFwd - matrix of numTimePts x numBouts for aligned fwd velocity
%       initFwd - numBouts length vector for initial forward velocity
%       changeFwd - numBouts length vector for change in forward velocity
%         from start to peak
%       peakYaw - numBouts length vector for peak forward velocity
%       t - time points for spike rate, yaw, fwd. 0 is yaw velocity peak
%       numBouts - total number of bouts (max n for stepParam (if no NaNs))
%       pDataFiles - struct of info on pData files
%           names - name of each pData file with at least 1 valid step, as
%               cell array
%           inds - indices (corresponding to bout indices) that
%               belong to each pData file, as cell array of vectors
%     FROM INPUT: cond, fwdVelCond, spikerateParams, postStimExclDur
%
% CREATED: 8/23/23 - HHY
%
% UPDATED:
%   8/23/23 - HHY
%
function saveSpikerateCond_bouts(cond, fwdVelCond, spikerateParams, ...
    flipLR, postStimExclDur, pDataFNames, pDataPath, saveFilePath, ...
    saveFileName)

    % check if pData files already specified
    if isempty(pDataFNames)
        % prompt user to select pData files
        [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
            'Select pData files for 1 fly', pDataPath, 'MultiSelect', 'on');
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
    pDataFiles.names = pDataFNames;
    pDataFiles.inds = [];

    allSpikerate = [];
    allYaw = [];
    allFwd = [];
    allInitFwd = [];
    allChangeFwd = [];
    allPeakYaw = [];


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
        % check if iInj trial, load those as appropriate
        % also, save iInj times as stim times; no stim, empty
        if (any(strcmpi(pDatVarsNames, 'iInj'))) % current inj
            load(pDataFullPath, 'ephysSpikes', 'moveNotMove', ...
                'fictracSmo', 'fictracProc', 'iInj');

            stimStartTimes = iInj.startTimes;
            stimEndTimes = iInj.endTimes;
        else % no stimulation
            load(pDataFullPath, 'ephysSpikes', 'moveNotMove', ...
                'fictracProc', 'fictracSmo');

            stimStartTimes = [];
            stimEndTimes = [];
        end

        % get yaw velocity peaks, with conditioning
        % specifies whether left or right turns
        [yawVelPeakTimes, boutStartTimes, boutEndTimes, ...
            yawVelPeakInd, boutStartInd, boutEndInd] = ...
            findCondYawVelPeaksFT(fictracSmo, cond, fwdVelCond, ...
            moveNotMove, ~flipLR);

        % remove turning bouts that overlap with stimulation bouts
        if (~isempty(stimStartTimes)) % only run on those with stimulation
            % get bouts that overlap
            boutInd2Rmv = findBoutsDurStim2Rmv(boutStartTimes, ...
                boutEndTimes, stimStartTimes, stimEndTimes, ...
                postStimExclDur);
            % remove overlapping bouts from all bout times/inds
            yawVelPeakTimes(boutInd2Rmv) = [];
            boutStartTimes(boutInd2Rmv) = [];
            boutEndTimes(boutInd2Rmv) = [];
            yawVelPeakInd(boutInd2Rmv) = [];
            boutStartInd(boutInd2Rmv) = [];
            boutEndInd(boutInd2Rmv) = [];
        end
        

        % check if this pData file contributes any turns
        % if not, skip and move to next pData file
        if (isempty(yawVelPeakInd))
            rmvInd = [rmvInd;i];
            continue;
        end

        % number of bouts for this trial
        thisNumBouts = length(yawVelPeakInd);

        % get indices for bouts for this trial
        thisTrialStartInd = 1 + countNumBouts;
        thisTrialEndInd = thisTrialStartInd + thisNumBouts - 1;

        thisTrialInds = [thisTrialStartInd thisTrialEndInd];
    
        pDataFiles.inds = [pDataFiles.inds; thisTrialInds];

        % update counter
        countNumBouts = countNumBouts + thisNumBouts;

        % get initial forward vel, change in forward vel, peak yaw vel
        % operates on fictracSmo
        initFwd = fictracSmo.fwdVel(boutStartInd);
        changeFwd = fictracSmo.fwdVel(yawVelPeakInd) - ...
            fictracSmo.fwdVel(boutStartInd);
        peakYaw = fictracSmo.yawAngVel(yawVelPeakInd);

        % concatenate
        allInitFwd = [allInitFwd; initFwd];
        allChangeFwd = [allChangeFwd; changeFwd];
        allPeakYaw = [allPeakYaw; peakYaw];



        % get spike rate, yaw vel, fwd velocity over time for each bout
        % use fictracProc
        % spike rate
        [spikerate, t] = getSpikerateFromBouts(ephysSpikes, ...
            spikerateParams, yawVelPeakTimes, boutStartTimes, ...
            boutEndTimes, 0);

        % yaw vel
        [yawVel, ~] = getFictracFromBouts(fictracProc, 'yawAngVel', ...
            spikerateParams, yawVelPeakInd, boutStartInd, boutEndInd);

        % fwd vel
        [fwdVelRight, ~] = getFictracFromBouts(fictracProc, 'fwdVel', ...
            spikerateParams, yawVelPeakInd, boutStartInd, boutEndInd);

        % concatenate all
        allSpikerate = cat(2, allSpikerate, spikerate);
        allYaw = cat(2, allYaw, yawVel);
        allFwd = cat(2, allFwd, fwdVelRight);
    end

    % remove unused pData files
    pDataFiles.names(rmvInd) = [];

    numBouts = countNumBouts;

    % save output file
    fullSavePath = [saveFilePath filesep saveFileName '.mat'];

    save(fullSavePath, 'allSpikerate', 'allYaw', 'allFwd',...
        'allInitFwd', 'allChangeFwd', 'allPeakYaw', 't', ...
         'numBouts', 'pDataFiles', 'cond', 'fwdVelCond', ...
         'spikerateParams', 'postStimExclDur','-v7.3');
end