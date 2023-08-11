% multiLinRegEphysFTContParam.m
%
% Function to perform multiple linear regression where an ephys parameter
%  is the independent variable and FicTrac and continuous step parameters 
%  are the dependent variables. Specifically, dependent variables are one
%  FicTrac parameter and one or more continuous step parameters. If there's
%  more than one step parameter, use the equal-weighted, mean-subtracted
%  combination of them. 
% This analyses asks if adding step parameter info can explain more of the
%  ephys response than the FicTrac parameter alone.
% Assumes that the chosen set of step parameters uses the same units.
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
%   fictracParam - string for which FicTrac parameter. Must be member of
%     fictracProc
%   stepParams - single string or cell array of strings for which step
%     parameters to correlate. Must be members of continuous
%     estimate of step parameters
%   legs - single string or cell array of strings for which legs for step
%     parameters.
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
% CREATED: 8/9/23 - HHY
%
% UPDATED:
%   8/9/23 - HHY
%
function multiLinRegEphysFTContParam(ephysParam, fictracParam, ...
    stepParams, legs, tDelay, notMoveExclDur, postStimExclDur, ...
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
    stepVals = []; % running tracker of all step time points
    fictracVals = []; % running tracker of all FicTrac time points


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
                'legStepsCont', 'ephysSpikes', 'iInj');
        else
            load(pDataFullPath, 'moveNotMove', 'fictracProc', ...
                'legStepsCont', 'ephysSpikes');
        end


        % get step parameter(s)
        % if multiple step parameters to input into regression
        if iscell(stepParams)
            % preallocate tracker for this trial's variables
            thisStepVars = [];

            % loop through all behavioral parameters
            for j = 1:length(stepParams)
                % get appropriate leg
                thisLegInd = legStepsCont.legIDs.ind(...
                    strcmpi(legStepsCont.legIDs.names, legs{j}));
    
                thisVal = legStepsCont.(stepParams{j})(:,thisLegInd);

                % add to running tracker for this trial
                thisStepVars = cat(2,thisStepVars,thisVal);
            end

        % if single step parameter to correlate with
        else
            % get appropriate leg
            thisLegInd = legStepsCont.legIDs.ind(...
                strcmpi(legStepsCont.legIDs.names, legs));

            thisStepVars = legStepsCont.(stepParams)(:,thisLegInd);
        end

        % get ephys values interpolated to leg timing
        thisEphysVal = interp1(ephysSpikes.t - tDelay, ephysSpikes.(ephysParam), ...
            legStepsCont.t, 'spline', nan);

        % get FicTrac values interpolated to leg timing
        thisFicTracVal = interp1(fictracProc.t, ...
            fictracProc.(fictracParam), legStepsCont.t, 'spline', nan);


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
                thisStepVars(thisNotMoveLog,:) = nan;
            else
                thisStepVars(thisNotMoveLog) = nan;
            end
            thisEphysVal(thisNotMoveLog) = nan;
            thisFicTracVal(thisNotMoveLog) = nan;
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
                    thisStepVars(thisIInjLog,:) = nan;
                else
                    thisStepVars(thisIInjLog) = nan;
                end
                thisEphysVal(thisIInjLog) = nan;
                thisFicTracVal(thisIInjLog) = nan;
            end
        end


        % remove NaNs
        rmInd = [];
        for j = 1:size(thisStepVars,1)
            % get all sample indices to remove
            if (any(isnan(thisStepVars(j,:))) || isnan(thisEphysVal(j)))
                rmInd = [rmInd; j];
            end
        end

        thisStepVars(rmInd,:) = [];
        thisEphysVal(rmInd) = [];
        thisFicTracVal(rmInd) = [];


        % add to variable tracker across pData files
        stepVals = cat(1,stepVals, thisStepVars);
        ephysVals = cat(1,ephysVals, thisEphysVal);
        fictracVals = cat(1,fictracVals, thisFicTracVal);
    end

    % remove outliers
    [~,stepOutliers] = rmoutliers(stepVals,"percentiles",[0.1 99.9], 1);
    [~,fictracOutliers] = rmoutliers(fictracVals,"percentiles",[0.1 99.9], 1);
    [~,ephysOutliers] = rmoutliers(ephysVals, "percentiles",[0.1 99.5]);
    if (strcmpi(ephysParam, 'spikeRate'))
        ephysL0 = ephysVals < 0;
    else
        ephysL0 = [];
    end

    outRmLog = stepOutliers | fictracOutliers | ephysOutliers | ephysL0;

    stepVals(outRmLog,:) = [];
    ephysVals(outRmLog) = [];
    fictracVals(outRmLog) = [];


    % if multiple step parameters, get projection to equal weighting
    if iscell(stepParams)
        % get coefficients
        coeffs = ones(length(stepParams), 1) * (1/length(stepParams));
        stepVals1D = getLinProj(stepVals, coeffs);
    % otherwise, just original behavioral value    
    else
        stepVals1D = stepVals;
    end

    % get linear model for equal projection of step parameters
    linMdl1D = fitlm(cat(2,fictracVals,stepVals1D),ephysVals);

    % get linear model, considering each step parameter individually
    linMdlAll = fitlm(cat(2,fictracVals,stepVals), ephysVals);

    % get robust linear fit, equal projections of step parameters
    [robustFit1D.coeff, robustFit1D.stats] = robustfit(...
        cat(2,fictracVals,stepVals1D), ephysVals);

    % get robust linear fit, each step parameter individually
    [robustFit.coeff, robustFit.stats] = robustfit(...
        cat(2,fictracVals,stepVals), ephysVals);

    % get names of step parameters (step parameters + legs)
    stepParamNames = cell(size(stepParams));
    for i = 1:length(stepParams)
        stepParamNames{i} = sprintf('%s %s', stepParams{i}, legs{i});
    end

    % save output file
    fullSavePath = [saveFileDir filesep saveFileName '.mat'];

    save(fullSavePath, 'stepVals', 'ephysVals', 'fictracVals', ...
        'stepVals1D', 'linMdl1D', 'linMdlAll','robustFit1D', ...
        'robustFit', 'ephysParam', 'stepParams', 'legs', ...
        'fictracParam', 'stepParamNames', ...
        'tDelay', 'notMoveExclDur', 'postStimExclDur','-v7.3');
end
