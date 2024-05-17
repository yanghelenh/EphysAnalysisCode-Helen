% getCorrelationEphysFwdAccYawVel_cond_fly.m
%
% Function to extract correlation between spike rate or median filtered
%  membrane potential and forward acceleration and yaw velocity.
% Adaptation of getCorrelationEphysContParam_cond_fly()
%
% User-specified filtering parameters to compute forward acceleration
% Normalization options on spike rate: fraction of max, change from
%  baseline (not walking), change relative to baseline as fraction of max
% Allows conditioning one or more specified parameters (FicTrac,
%  legStepsCont, or binary for whether step is moving forward or
%  backwards). 
% Removes not-moving times and window around transitions (user specified).
% Removes current injection and window after injection (user specified).
% Returns correlation between variables with user specified delay between
%  ephys and behavior
% Operates user specified sampling rate
% Runs on multiple pData files from one fly (specified either
%  through input or through GUI
%
% INPUTS:
%   ephysParam - string for which ephys parameter ('spikeRate' or
%     'medFiltV')
%   cond - struct of conditions that time points have to meet to be
%     included. Multiple conditions are AND. [] if no conditions
%       whichParam - cell array (even if 1 element) on which fictracProc or
%           legStepsCont fields to condition on. Use 'stepFwdBool' for
%           boolean of whether step is moving backwards or forwards. 
%       legs - cell array for which leg for legStepsCont and 'stepFwdBool'.
%         [] for FicTrac elements
%       cond - cell array of strings to condition on, for eval(); same size
%           as whichParam
%   normSpikeRate - struct of parameters for normalizing spike rate. Empty
%     if no normalization across flies
%       method - 'max' for fraction of max, 'baseline' for change relative
%         to not walking, 'maxBase' as change relative to walking as
%         fraction of max
%       val - for 'max', use val percentile as max. For 'baseline', time in
%         sec at beginning and end of not moving bout to exclude from
%         consideration. As 2 element vector. For 'maxBase', 3 element
%         vector with percentile first, times second
%   fwdAccFilt - struct of filtering parameters for computing forward
%     acceleration (Gaussian filtering)
%       padLen - pad length, in s
%       sigmaVel - sigma for forward velocity, applied before
%       sigmaAcc - sigma for forward accleration, applied after
%   tDelay - time offset b/w behavior and ephys, in seconds. Negative for 
%     ephys before behavior
%   sampRate - sampling rate, in Hz
%   notMoveExclDur - additional time, in sec, before and after not
%     moving bout to exclude from consideration. Length 2 vector for before
%     and after, respectively
%   postStimExclDur - additional time, in sec, after iInj stimulation 
%     to exclude from consideration
%   pDataPath - full path to pData files
%   pDataFNames - string or cell array of strings for pData files to operate
%     on. If [], then select pData files through GUI. One fly only 
%   saveFileName - name of output file
%   saveFileDir - full path to directory to save output file
%
% OUTPUTS:
%   fwdAccVals
%   yawVelVals
%   ephysVals
%   physValsNorm
%   saves output file
%
% CREATED: 5/15/24 - HHY
%
% UPDATED:
%   5/15/24 - HHY
%
function [fwdAccVals, yawVelVals, ephysVals, ephysValsNorm] = ...
    getCorrelationEphysFwdAccYawVel_cond_fly(ephysParam, cond, ...
    normSpikeRate, fwdAccFilt, tDelay, sampRate, notMoveExclDur, ...
    postStimExclDur, pDataPath, pDataFNames, saveFileName, saveFileDir)
    
    % prompt user to select pData files
    if isempty(pDataFNames)
        [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
            'Select pData files', pDataPath, 'MultiSelect', 'on');
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


    % preallocate
    ephysVals = []; % running tracker of all ephys time points
    fwdAccVals = []; % running tracker of all fwd acceleration time points
    yawVelVals = []; % running tracker of all yaw velocity time points
    % running tracker of all ephys time points during not moving
    notMoveEphysVals = []; 


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

        % check if this pData file has legSteps, bodytraj, moveNotMove
        %  structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'fictracProc')) || ...
                ~any(strcmpi(pDatVarsNames, 'moveNotMove')) || ...
                ~any(strcmpi(pDatVarsNames, 'ephysSpikes')))
            continue;
        end

        % load variables from pData
        if (any(strcmpi(pDatVarsNames, 'iInj')))
            load(pDataFullPath, 'moveNotMove', 'fictracProc', ...
                'ephysSpikes', 'legSteps', 'iInj');
        else
            load(pDataFullPath, 'moveNotMove', 'fictracProc', ...
                'ephysSpikes', 'legSteps');
        end


        % get forward acceleration - gradient of forward velocity, with
        %  gaussian filtering applied on velocity and on acceleration
        % get filtering parameters, in samples
        ftIFI = median(diff(fictracProc.t));
        padLen = round(fwdAccFilt.padLen * 1/ftIFI);
        sigmaVel = round(fwdAccFilt.sigmaVel * 1/ftIFI);
        sigmaAcc = round(fwdAccFilt.sigmaAcc * 1/ftIFI);

        filtFwdVel = gaussSmooth(fictracProc.fwdVel,padLen,sigmaVel);
        unFiltFwdAcc = gradient(filtFwdVel) .* (1/ftIFI);
        filtFwdAcc = gaussSmooth(unFiltFwdAcc,padLen,sigmaAcc);

        % get sample times
        sampT = fictracProc.t(1):(1/sampRate):fictracProc.t(end);
        sampT = sampT'; % column vector


        % get forward acceleration and yaw velocity, interpolated to sample
        %  rate
        thisFwdAcc = interp1(fictracProc.t, filtFwdAcc, sampT, 'spline');
        thisYawVel = interp1(fictracProc.t, fictracProc.yawAngVel, ...
            sampT, 'spline');

        % get ephys values interpolated to sample rate
        % also, incorporate time offset b/w ephys and behavior
        thisEphysVal = interp1(ephysSpikes.t - tDelay, ...
            ephysSpikes.(ephysParam), sampT, 'spline', nan);
        % initalize ephys vals for not moving
        if contains(normSpikeRate.method,'base', 'IgnoreCase',true)
            thisEphysValNotMove = thisEphysVal;
        end

        % get not moving times, including additional exclusion duration
        %  set those to NaN
        notMoveStartTimes = ...
            moveNotMove.legT(moveNotMove.legNotMoveBout(:,1)) - ...
            notMoveExclDur(1);
        notMoveEndTimes = ...
            moveNotMove.legT(moveNotMove.legNotMoveBout(:,2)) + ...
            notMoveExclDur(2);


        % loop through all bouts, set appropriate inds to NaN
        for j = 1:length(notMoveStartTimes)
            thisNotMoveLog = (sampT > notMoveStartTimes(j)) & ...
                (sampT < notMoveEndTimes(j));
            % set to NaN
            thisFwdAcc(thisNotMoveLog) = nan;
            thisYawVel(thisNotMoveLog) = nan;
            thisEphysVal(thisNotMoveLog) = nan;
        end

        % when spike rate normalized to baseline during not walking
        % get ephys val when not moving
        if contains(normSpikeRate.method,'base', 'IgnoreCase',true)          
            % keep track of all not moving, as logical
            allNotMoveLog = false(size(thisEphysValNotMove));
    
            % not move start and end times, account for buffer
            if(strcmpi(normSpikeRate.method, 'baseline'))
                notMoveStartTimes = ...
                    moveNotMove.legT(moveNotMove.legNotMoveBout(:,1)) + ...
                    normSpikeRate.val(1);
                notMoveEndTimes = ...
                    moveNotMove.legT(moveNotMove.legNotMoveBout(:,2)) - ...
                    normSpikeRate.val(2);
            elseif(strcmpi(normSpikeRate.method, 'maxBase'))
                notMoveStartTimes = ...
                    moveNotMove.legT(moveNotMove.legNotMoveBout(:,1)) + ...
                    normSpikeRate.val(2);
                notMoveEndTimes = ...
                    moveNotMove.legT(moveNotMove.legNotMoveBout(:,2)) - ...
                    normSpikeRate.val(3);
            end

            % loop through all not moving bouts
            for j = 1:length(notMoveStartTimes)
                thisNotMoveLog = (sampT > notMoveStartTimes(j)) & ...
                    (sampT < notMoveEndTimes(j));

                % update logical for all not moving
                allNotMoveLog = allNotMoveLog | thisNotMoveLog;
            end
    
            % update ephys vals for not moving (NaN all ind not in not 
            %  moving)
            thisEphysValNotMove(~allNotMoveLog) = nan;
        end

        % set current injection times to NaN
        if (any(strcmpi(pDatVarsNames, 'iInj')))
            % current injection times; only add to end
            iInjStartTimes = iInj.startTimes;
            iInjEndTimes = iInj.endTimes + postStimExclDur;

            % loop through all iInj bouts, set appropriate inds to NaN
            for j = 1:length(iInjStartTimes)
                thisIInjLog = (sampT > iInjStartTimes(j)) & ...
                    (sampT < iInjEndTimes(j));

                % set to NaN
                thisFwdAcc(thisIInjLog) = nan;
                thisYawVel(thisIInjLog) = nan;
                thisEphysVal(thisIInjLog) = nan;
                thisEphysValNotMove(thisIInjLog) = nan;               
            end
        end

        % set points that don't meet conditions to NaN
        if ~isempty(cond)
            % initialize logical to track which points meet condition
            % loop through all conditions
            condLog = true(size(sampT));
            for j = 1:length(cond.whichParam)
                % check which data structures conditions belong to
                % check first for stepFwdBool, as it's not a field of any
                %  struct
                if (strcmpi(cond.whichParam{j}, 'stepFwdBool'))
                    thisFwdLog = getStepFwdLogical(legSteps, sampT,...
                        cond.legs{j});
                    if ~(eval(cond.cond{j})) % if target is false, invert
                        thisLog = ~thisFwdLog;
                    else
                        thisLog = thisFwdLog;
                    end
                % if FicTrac    
                elseif isfield(fictracProc,cond.whichParam{j})
                    thisFT = interp1(fictracProc.t, ...
                        fictracProc.(cond.whichParam{j}), sampT, ...
                        'spline');
                    thisLog = eval(['thisFT' cond.cond{j}]); 
                % if step parameter    
                elseif isfield(legStepsCont,cond.whichParam{j})
                    thisLegInd = legSteps.legIDs.ind(strcmpi(legSteps.legIDs.names, ...
                        cond.legs{j}));
                    thisStep = legStepsCont.(cond.whichParam{j})(:,thisLegInd);
                    thisLog = eval(['thisStep' cond.cond{j}]);
                end
                % update logical for all conditions
                condLog = condLog & thisLog;
            end
            % invert condition logical
            condFalseLog = ~condLog;
            % set every time when condition logcial false to nan
            thisFwdAcc(condFalseLog,:) = nan;
            thisYawVel(condFalseLog,:) = nan;
            thisEphysVal(condFalseLog) = nan;

        end



        % remove NaNs
        rmInd = [];
        for j = 1:size(thisFwdAcc,1)
            % get all sample indices to remove
            if (isnan(thisFwdAcc(j)) || isnan(thisYawVel(j))...
                    || isnan(thisEphysVal(j)))
                rmInd = [rmInd; j];
            end
        end

        thisFwdAcc(rmInd) = [];
        thisYawVel(rmInd) = [];
        thisEphysVal(rmInd) = [];


        % add to variable tracker across pData files
        fwdAccVals = cat(1,fwdAccVals, thisFwdAcc);
        yawVelVals = cat(1,yawVelVals, thisYawVel);
        ephysVals = cat(1,ephysVals, thisEphysVal);

        % if relative to baseline
        if contains(normSpikeRate.method,'base', 'IgnoreCase',true)
            % remove NaNs from not move
            thisEphysValNotMove(isnan(thisEphysValNotMove)) = [];
            notMoveEphysVals = cat(1,notMoveEphysVals, thisEphysValNotMove);
        end
    end


    % remove outliers
    [~,fwdAccOutliers] = rmoutliers(fwdAccVals,"percentiles",[0.1 99.9]);
    [~,yawVelOutliers] = rmoutliers(yawVelVals,"percentiles",[0.1 99.9]);
    [~,ephysOutliers] = rmoutliers(ephysVals, "percentiles",[0 99.5]);
    % eliminate any time points where ephys spike rate < 0
    if (strcmpi(ephysParam, 'spikeRate'))
        ephysL0 = ephysVals < 0;
        if contains(normSpikeRate.method,'base', 'IgnoreCase',true)
            notMoveEphysVals(notMoveEphysVals < 0) = [];
        end
    else
        ephysL0 = [];
    end

    outRmLog = fwdAccOutliers | yawVelOutliers | ephysOutliers | ephysL0;

    fwdAccVals(outRmLog) = [];
    yawVelVals(outRmLog) = [];
    ephysVals(outRmLog) = [];


    % normalize spike rate
    % check that we want to first
    if (strcmpi(ephysParam, 'spikeRate') && ~isempty(normSpikeRate))
        % normalization as fraction of max value (percentile)
        if (strcmpi(normSpikeRate.method, 'max'))
            maxSpikeRate = prctile(ephysVals, normSpikeRate.val);

            % norm rate as fraction of max rate
            ephysValsNorm = ephysVals / maxSpikeRate;

        % normalization relative to baseline spike rate, i.e. avg spike
        %  rate during not moving period
        elseif (strcmpi(normSpikeRate.method, 'baseline'))
            baseSpikeRate = mean(notMoveEphysVals);
    
            % norm rate as baseline subtracted 
            ephysValsNorm = ephysVals - baseSpikeRate;

        % normalization as fraction of max (percentile), baseline subtracted    
        elseif (strcmpi(normSpikeRate.method, 'maxBase'))
            baseSpikeRate = mean(notMoveEphysVals);
            maxSpikeRate = prctile(ephysVals, normSpikeRate.val(1));


            % norm rate as fraction of max, baseline subtracted
            ephysValsNorm = (ephysVals - baseSpikeRate) / ...
                (maxSpikeRate - baseSpikeRate);
        end
    end


    % save output file
    fullSavePath = [saveFileDir filesep saveFileName '.mat'];

    save(fullSavePath, 'fwdAccVals', 'yawVelVals', 'ephysVals', ...
        'ephysValsNorm', 'ephysParam', 'tDelay', 'fwdAccFilt', 'sampRate', ...
        'notMoveExclDur', 'postStimExclDur', 'cond', 'normSpikeRate', ...
        '-v7.3');
end
