% getContParamAutocorr_cell.m
%
% Function that computes the autocorrelation of a continuous leg step 
%  parameter or FicTrac parameterfor a single cell, only during walking 
%  and only under additional conditions on behavioral parameters, 
%  specified by cond.
% Across trials within the same fly, weight by number of valid points
%  during that trial
% Select input pData files for single cell through GUI or as input
%
% INPUTS:
%   behParams - single string or cell array of strings for which behavioral
%     parameters to get autocorrelation for. Must be members of FicTrac or 
%     continuous estimate of step parameters
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
%   maxLag - max lag, in sec, to compute autocorrelation for
%   notMoveExclDur - additional time, in sec, before and after not
%     moving bout to exclude from consideration. Length 2 vector for before
%     and after, respectively
%   postStimExclDur - additional time, in sec, after iInj stimulation 
%     to exclude from consideration
%   pDataFNames - cell array of pData file names or [] if select through
%       GUI
%   pDataPath - path to folder containing pData files
%   saveFileDir - full path to folder in which to save output file
%
% OUTPUTS:
%   none, but saves file
%
% CREATED: 9/1/23 - HHY
%
% UPDATED:
%   9/1/23 - HHY
%
function getContParamAutocorr_cell(behParams, legs, cond, maxLag, ...
    notMoveExclDur, postStimExclDur, pDataFNames, pDataPath, saveFileDir)

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

    % preallocate
    aCorrAll = [];
    trialCount = zeros(1, numPDataFiles);


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

        % save fly name as first pDataName's date, fly, cell (19 characters)
        if (i == 1)
            flyName = pDataName(1:19);
        end

        % check if this pData file has legSteps, fictracProc, moveNotMove
        %  structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legStepsCont')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracProc')) || ...
                ~any(strcmpi(pDatVarsNames, 'moveNotMove')) || ...
                ~any(strcmpi(pDatVarsNames, 'legSteps')))
            continue;
        end

        % load variables from pData
        if (any(strcmpi(pDatVarsNames, 'iInj')))
            load(pDataFullPath, 'moveNotMove', 'fictracProc', ...
                'legStepsCont', 'legSteps', 'iInj');
        else
            load(pDataFullPath, 'moveNotMove', 'fictracProc', ...
                'legStepsCont', 'legSteps');
        end

        % get behavioral variable(s)
        % if multiple behavioral parameters to autocorrelate
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

            % get equal weighting of behavioral params
            coeffs = ones(length(behParams), 1) * (1/length(behParams));
            thisParamVal = getLinProj(thisBehVars, coeffs);
            thisT = legStepsCont.t;

        % if single behavioral parameter to correlate with
        else
            % FicTrac param
            if isempty(legs)
                thisParamVal = fictracProc.(behParams);
                thisT = fictracProc.t;
            else
                % get appropriate leg
                thisLegInd = legStepsCont.legIDs.ind(...
                    strcmpi(legStepsCont.legIDs.names, legs));

                    thisParamVal = legStepsCond.(behParams)(:,thisLegInd);
                    thisT = legStepsCond.t;
            end
        end

        % get not moving times, set those to NaN
        notMoveStartTimes = moveNotMove.ftT(moveNotMove.ftNotMoveBout(:,1)) - ...
            notMoveExclDur(1);
        notMoveEndTimes = moveNotMove.ftT(moveNotMove.ftNotMoveBout(:,2)) + ...
            notMoveExclDur(2);


        % loop through all bouts, set appropriate inds to NaN
        for j = 1:length(notMoveStartTimes)
            thisNotMoveLog = (thisT > notMoveStartTimes(j)) & ...
                (thisT < notMoveEndTimes(j));
            % set to NaN
            thisParamVal(thisNotMoveLog) = nan;
        end


        % set current injection times to NaN
        if (any(strcmpi(pDatVarsNames, 'iInj')))
            % current injection times; only add to end
            iInjStartTimes = iInj.startTimes;
            iInjEndTimes = iInj.endTimes + postStimExclDur;

            % loop through all iInj bouts, set appropriate inds to NaN
            for j = 1:length(iInjStartTimes)
                thisIInjLog = (thisT > iInjStartTimes(j)) & ...
                    (thisT < iInjEndTimes(j));

                % set to NaN
                thisParamVal(thisIInjLog) = nan;         
            end
        end

        % set points that don't meet conditions to NaN
        if ~isempty(cond)
            % initialize logical to track which points meet condition
            % loop through all conditions
            condLog = true(size(thisT));
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
                        fictracProc.(cond.whichParam{j}), thisT, ...
                        'spline');
                    thisLog = eval(['thisFT' cond.cond{j}]); 
                % if step parameter    
                elseif isfield(legStepsCont,cond.whichParam{j})
                    thisLegInd = legSteps.legIDs.ind(strcmpi(legSteps.legIDs.names, ...
                        cond.legs{j}));
                    thisStep = interp1(legStepsCont.t, ...
                        legStepsCont.(cond.whichParam{j})(:,thisLegInd), ...
                        thisT, 'spline');
                    thisLog = eval(['thisStep' cond.cond{j}]);
                end
                % update logical for all conditions
                condLog = condLog & thisLog;
            end
            % invert condition logical
            condFalseLog = ~condLog;
            % set every time when condition logcial false to nan
            thisParamVal(condFalseLog) = nan;
        end


        % get autocorrelation

        % get number of lags
        ifi = median(diff(thisT));
        % number of lags to one side
        numLags1Side = ceil(maxLag / ifi);
        % number of lags, total
        numLags = numLags1Side * 2 + 1;

        % autocorrelation, with NaNs allowed
        [allCorr, lags] = xcorrWGaps(thisParamVal, thisParamVal, numLags);

        % number of points this trial contributes
        trialCount(i) = length(thisParamVal(~isnan(thisParamVal)));

        % add autocorr to running tracker across trials
%         size(aCorrAll)
%         size(allCorr)
        aCorrAll = [aCorrAll, allCorr];
    end

    % get weighting for each trial
    totPts = sum(trialCount);
    trialWeight = trialCount / totPts;

    % matrix of weighting, same value in each column
    trialWeightMatrix = repmat(trialWeight, numLags, 1);

    % weight autocorrelations
    weightAutoCorrAll = aCorrAll .* trialWeightMatrix;

    % sum across rows to get weighted autocorrelation
    autoCorr = zeros(size(weightAutoCorrAll,1),1);
    for i = 1:size(weightAutoCorrAll,1)
        thisRow = weightAutoCorrAll(i,:);
        autoCorr(i) = sum(thisRow(~isnan(thisRow)));
    end
%     autoCorr = sum(weightAutoCorrAll, 2);

    % convert lags in samples to t
    lagsT = lags * ifi;

    % save data to output file
    saveFileFullName = [saveFileDir filesep flyName '_autocorr.mat'];
    save(saveFileFullName, 'behParams', 'legs', 'autoCorr', 'lags', 'lagsT', ...
        'cond', 'maxLag', 'notMoveExclDur', 'postStimExclDur', 'aCorrAll',...
        '-v7.3');
end