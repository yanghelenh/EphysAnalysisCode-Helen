% getCorrelationEphysContParam_cond.m
%
% Variant of getCorrelationEphysContParam() (description below) but it also
%  allows conditioning one or more specified parameters (FicTrac,
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
% Runs on multiple pData files from multiple flies (specified either
%  through input or through GUI
%
% INPUTS:
%   ephysParam - string for which ephys parameter ('spikeRate' or
%     'medFiltV')
%   behParam - single string or cell array of strings for which behavioral
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
%   tDelay - time offset b/w behavior and ephys. Negative for ephys before
%     behavior
%   notMoveExclDur - additional time, in sec, before and after not
%     moving bout to exclude from consideration. Length 2 vector for before
%     and after, respectively
%   postStimExclDur - additional time, in sec, after iInj stimulation 
%     to exclude from consideration
%   pDataPath - full path to pData files
%   pDataFNames - string or cell array of strings for pData files to operate
%     on. If [], then select pData files through GUI. 
%   saveFileName - name of output file
%   saveFileDir - full path to directory to save output file
%
% OUTPUTS:
%   none, but saves output file with following variables:
%
% CREATED: 8/7/23 - HHY
%
% UPDATED:
%   8/7/23 - HHY
%
function getCorrelationEphysContParam_cond(ephysParam, behParams, legs, ...
    cond, tDelay, notMoveExclDur, postStimExclDur, pDataPath, ...
    pDataFNames, saveFileName, saveFileDir)
    
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
    end

    % remove outliers
    [~,behOutliers] = rmoutliers(behVals,"percentiles",[0.1 99.9], 1);
    [~,ephysOutliers] = rmoutliers(ephysVals, "percentiles",[0.1 99.5]);
    if (strcmpi(ephysParam, 'spikeRate'))
        ephysL0 = ephysVals < 0;
    else
        ephysL0 = [];
    end

    outRmLog = behOutliers | ephysOutliers | ephysL0;

    behVals(outRmLog,:) = [];
    ephysVals(outRmLog) = [];


    % if multiple behavioral variables, get projection to equal weighting
    if iscell(behParams)
        % get coefficients
        coeffs = ones(length(behParams), 1) * (1/length(behParams));
        behVals1D = getLinProj(behVals, coeffs);
    % otherwise, just original behavioral value    
    else
        behVals1D = behVals;
    end

    % get linear model
    linMdl = fitlm(behVals1D,ephysVals);

    % save output file
    fullSavePath = [saveFileDir filesep saveFileName '.mat'];

    save(fullSavePath, 'behVals', 'ephysVals', 'behVals1D', 'linMdl', ...
        'ephysParam', 'behParams', 'legs', 'tDelay', 'notMoveExclDur', ...
        'postStimExclDur', 'cond', '-v7.3');
end
