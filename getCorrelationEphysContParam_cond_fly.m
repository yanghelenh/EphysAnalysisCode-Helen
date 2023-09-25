% getCorrelationEphysContParam_cond_fly.m
%
% Variant of getCorrelationEphysContParam_cond() that operates on one fly
%  at a time and allows normalization of spike rate (as
%  fraction of max spike rate (max as percentile) or as change from
%  baseline (when fly isn't walking). Note that this also means that mean
%  subtracted projections of step parameters will differ across flies
% 
% Allows conditioning one or more specified parameters (FicTrac,
%  legStepsCont, or binary for whether step is moving forward or
%  backwards).
% Function to extract correlation between spike rate or median filtered
%  membrane potential and one or more continous parameters (from continuous
%  estimate of step parameters or FicTrac variables. If more than one 
%  parameter is chosen (i.e. more than 1 parameter and 1 leg), correlate 
%  specified ephys parameter with the equal-weighted, mean-subtracted 
%  combination of them. Assumes that the chosen set of parameters uses the
%  same units.
% Removes not-moving times and window around transitions (user specified).
% Returns correlation between variables with user specified delay between
%  ephys and steps
% Operates on timescale of legTrack frame rate (250 Hz)
% Runs on multiple pData files from one fly (specified either
%  through input or through GUI
%
% INPUTS:
%   ephysParam - string for which ephys parameter ('spikeRate' or
%     'medFiltV')
%   behParams - single string or cell array of strings for which behavioral
%     parameters to correlate. Must be members of FicTrac or continuous
%     estimate of step parameters
%   legs - single string or cell array of strings for which legs, when step
%     parameters. [] when FicTrac param
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
%   tDelay - time offset b/w behavior and ephys. Negative for ephys before
%     behavior
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
%   none, but saves output file with following variables:
%
% CREATED: 8/11/23 - HHY
%
% UPDATED:
%   8/11/23 - HHY
%   9/25/23 - HHY - add handling of circular parameters
%
function [behVals, ephysVals, ephysValsNorm] = ...
    getCorrelationEphysContParam_cond_fly(ephysParam, behParams, ...
    legs, cond, normSpikeRate, tDelay, notMoveExclDur, postStimExclDur, ...
    pDataPath, pDataFNames, saveFileName, saveFileDir)
    
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
    behVals = []; % running tracker of all behavior time points
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
        if (~any(strcmpi(pDatVarsNames, 'legStepsCont')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracProc')) || ...
                ~any(strcmpi(pDatVarsNames, 'moveNotMove')) || ...
                ~any(strcmpi(pDatVarsNames, 'ephysSpikes')))
            continue;
        end

        % load variables from pData
        if (any(strcmpi(pDatVarsNames, 'iInj')))
            load(pDataFullPath, 'moveNotMove', 'fictracProc', ...
                'legStepsCont', 'ephysSpikes', 'legSteps', 'iInj');
        else
            load(pDataFullPath, 'moveNotMove', 'fictracProc', ...
                'legStepsCont', 'ephysSpikes', 'legSteps');
        end


        % get behavioral variable(s)
        % if multiple behavioral parameters to correlate ephys with
        if iscell(behParams)
            % preallocate tracker for this trial's variables
            thisBehVars = [];

            % loop through all behavioral parameters
            for j = 1:length(behParams)
                % get this behavioral variable's values
                % leg or FicTrac parameter - no leg val for FicTrac
                if isempty(legs{j})
                    % interpolate FicTrac values to leg timing
                    thisVal = interp1(fictracProc.t, ...
                        fictracProc.(behParams{j}), legStepsCont.t, ...
                        'spline');
                else
                    % get appropriate leg
                    thisLegInd = legStepsCont.legIDs.ind(...
                        strcmpi(legStepsCont.legIDs.names, legs{j}));

                    thisVal = legStepsCont.(behParams{j})(:,thisLegInd);
                end

                % add to running tracker for this trial
                thisBehVars = cat(2,thisBehVars,thisVal);
            end

        % if single behavioral parameter to correlate with
        else
            % FicTrac param
            if isempty(legs)
                % interpolate FicTrac values to leg timing
                thisBehVars = interp1(fictracProc.t, ...
                    fictracProc.(behParams), legStepsCont.t, ...
                    'spline');
            else
                % get appropriate leg
                thisLegInd = legStepsCont.legIDs.ind(...
                    strcmpi(legStepsCont.legIDs.names, legs));

                thisBehVars = legStepsCont.(behParams)(:,thisLegInd);
            end
        end

        % get ephys values interpolated to leg timing
        thisEphysVal = interp1(ephysSpikes.t - tDelay, ephysSpikes.(ephysParam), ...
            legStepsCont.t, 'spline', nan);
        % initalize ephys vals for not moving
        if contains(normSpikeRate.method,'base', 'IgnoreCase',true)
            thisEphysValNotMove = thisEphysVal;
        end

        % get not moving times, set those to NaN
        notMoveStartTimes = moveNotMove.legT(moveNotMove.legNotMoveBout(:,1)) - ...
            notMoveExclDur(1);
        notMoveEndTimes = moveNotMove.legT(moveNotMove.legNotMoveBout(:,2)) + ...
            notMoveExclDur(2);


        % loop through all bouts, set appropriate inds to NaN
        for j = 1:length(notMoveStartTimes)
            thisNotMoveLog = (legStepsCont.t > notMoveStartTimes(j)) & ...
                (legStepsCont.t < notMoveEndTimes(j));
            % set to NaN
            if (size(thisNotMoveLog,2) > 1)
                thisBehVars(thisNotMoveLog,:) = nan;
            else
                thisBehVars(thisNotMoveLog) = nan;
            end
            thisEphysVal(thisNotMoveLog) = nan;
        end

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
                thisNotMoveLog = (legStepsCont.t > notMoveStartTimes(j)) & ...
                    (legStepsCont.t < notMoveEndTimes(j));

                % update logical for all not moving
                allNotMoveLog = allNotMoveLog | thisNotMoveLog;
            end
    
            % update ephys vals for not moving (NaN all ind not in not moving)
            thisEphysValNotMove(~allNotMoveLog) = nan;
        end

        % set current injection times to NaN
        if (any(strcmpi(pDatVarsNames, 'iInj')))
            % current injection times; only add to end
            iInjStartTimes = iInj.startTimes;
            iInjEndTimes = iInj.endTimes + postStimExclDur;

            % loop through all iInj bouts, set appropriate inds to NaN
            for j = 1:length(iInjStartTimes)
                thisIInjLog = (legStepsCont.t > iInjStartTimes(j)) & ...
                    (legStepsCont.t < iInjEndTimes(j));

                % set to NaN
                if (size(thisIInjLog,2) > 1)
                    thisBehVars(thisIInjLog,:) = nan;
                else
                    thisBehVars(thisIInjLog) = nan;
                end
                thisEphysVal(thisIInjLog) = nan;
                thisEphysValNotMove(thisIInjLog) = nan;               
            end
        end

        % set points that don't meet conditions to NaN
        if ~isempty(cond)
            % initialize logical to track which points meet condition
            % loop through all conditions
            condLog = true(size(legStepsCont.t));
            for j = 1:length(cond.whichParam)
                % check which data structures conditions belong to
                % check first for stepFwdBool, as it's not a field of any
                %  struct
                if (strcmpi(cond.whichParam{j}, 'stepFwdBool'))
                    thisFwdLog = getStepFwdLogical(legSteps, legStepsCont.t,...
                        cond.legs{j});
                    if ~(eval(cond.cond{j})) % if target is false, invert
                        thisLog = ~thisFwdLog;
                    else
                        thisLog = thisFwdLog;
                    end
                % if FicTrac    
                elseif isfield(fictracProc,cond.whichParam{j})
                    thisFT = interp1(fictracProc.t, ...
                        fictracProc.(cond.whichParam{j}), legStepsCont.t, ...
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
            thisBehVars(condFalseLog,:) = nan;
            thisEphysVal(condFalseLog) = nan;

        end



        % remove NaNs
        rmInd = [];
        for j = 1:size(thisBehVars,1)
            % get all sample indices to remove
            if (any(isnan(thisBehVars(j,:))) || isnan(thisEphysVal(j)))
                rmInd = [rmInd; j];
            end
        end

        thisBehVars(rmInd,:) = [];
        thisEphysVal(rmInd) = [];


        % add to variable tracker across pData files
        behVals = cat(1,behVals, thisBehVars);
        ephysVals = cat(1,ephysVals, thisEphysVal);
        if contains(normSpikeRate.method,'base', 'IgnoreCase',true)
            % remove NaNs from not move
            thisEphysValNotMove(isnan(thisEphysValNotMove)) = [];
            notMoveEphysVals = cat(1,notMoveEphysVals, thisEphysValNotMove);
        end
    end


    % remove outliers
    [~,behOutliers] = rmoutliers(behVals,"percentiles",[0.1 99.9], 1);
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

    outRmLog = behOutliers | ephysOutliers | ephysL0;

    behVals(outRmLog,:) = [];
    ephysVals(outRmLog) = [];



    % if multiple behavioral variables, get projection to equal weighting
    if iscell(behParams)
        % handle circular parameter correctly
        % assume all behParams are circular
        if any(strcmpi(behParams, 'stepDirection'))
            numParams = length(behParams);

            meanSubBehVals = behVals;
            for i = 1:numParams
                thisBehVals = behVals(:,i);
                meanSubBehVals(:,i) = thisBehVals - ...
                    rad2deg(circ_mean(deg2rad(thisBehVals)));
            end

            behVals1D = wrapTo180(rad2deg(circ_mean(deg2rad(meanSubBehVals), [], 2)));

        % not circular parameter, project    
        else
            % get coefficients
            coeffs = ones(length(behParams), 1) * (1/length(behParams));
            behVals1D = getLinProj(behVals, coeffs);
        end
    % otherwise,    
    else
        % just original behavioral value 
%         behVals1D = behVals;

        % deal with circular parameters appropriately
        if strcmpi(behParams,'stepDirection')
            behVals1D = behVals - rad2deg(circ_mean(deg2rad(behVals)));
        else      
            % mean subtracted behavioral value
            behVals1D = behVals - mean(behVals);
        end
    end

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

    save(fullSavePath, 'behVals', 'ephysVals', 'ephysValsNorm', ...
        'behVals1D', 'ephysParam', 'behParams', 'legs', 'tDelay', ...
        'notMoveExclDur', 'postStimExclDur', 'cond', 'normSpikeRate', ...
        '-v7.3');
end
