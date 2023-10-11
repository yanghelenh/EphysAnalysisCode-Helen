% getCorrCoeffBehContParam_cond_fly.m
% 
% Function to get Pierson's correlation coefficient between two specified
%  behavioral parameters, while conditioning on others, for a single fly.
% Modification of getCorrelationEphysContParam_cond_fly()
% 
% INPUTS:
%   behParams1 - struct of parameters for first behavioral parameter to
%     correlate
%       whichParams - cell array of string(s) for which behavioral
%         parameters to correlate. Must be members of FicTrac or continuous
%         estimate of step parameters
%       legs - single string or cell array of strings for which legs, when step
%         parameters. [] when FicTrac param
%   behParams2 - struct of parameters for second behavioral parameter to
%     correlate. Same fields as behParams1
%   cond - struct of conditions that time points have to meet to be
%     included. Multiple conditions are AND. [] if no conditions
%       whichParam - cell array (even if 1 element) on which fictracProc or
%           legStepsCont fields to condition on. Use 'stepFwdBool' for
%           boolean of whether step is moving backwards or forwards. 
%       legs - cell array for which leg for legStepsCont and 'stepFwdBool'.
%         [] for FicTrac elements
%       cond - cell array of strings to condition on, for eval(); same size
%           as whichParam
%   notMoveExclDur - additional time, in sec, before and after not
%     moving bout to exclude from consideration. Length 2 vector for before
%     and after, respectively
%   postStimExclDur - additional time, in sec, after iInj stimulation 
%     to exclude from consideration
%   pDataPath - full path to pData files
%   pDataFNames - string or cell array of strings for pData files to operate
%     on. 
%
% OUTPUTS:
%   r - Pierson's correlation coefficient between two behavioral parameters
%
% CREATED: 10/7/23 - HHY
%
% UPDATED:
%   10/7/23 - HHY
%
function r = getCorrCoeffBehContParam_cond_fly(behParams1, behParams2, ...
    cond, notMoveExclDur, postStimExclDur, pDataPath, pDataFNames)

    % if only 1 pData file, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end


    % preallocate
    behVals1 = []; % running tracker of all behavior time points
    behVals2 = [];


    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end

        pDataFullPath = [pDataPath filesep pDataName];

        % get variables saved in pData file
        pDatVars = whos('-file', pDataFullPath);
    
        pDatVarsNames = cell(size(pDatVars));
        
        % convert pDatVars into cell array of just names
        for j = 1:length(pDatVars)
            pDatVarsNames{j} = pDatVars(j).name;
        end

        % check if this pData file has legStepsCont, fictracProc, moveNotMove
        %  moveNotMove structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'legStepsCont')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracProc')) || ...
                ~any(strcmpi(pDatVarsNames, 'moveNotMove')))
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


        % get behavioral variable 1
        % if multiple behavioral parameters to correlate ephys with
        if iscell(behParams1.whichParams)
            % preallocate tracker for this trial's variables
            thisBehVars1 = [];

            % loop through all behavioral parameters
            for j = 1:length(behParams1.whichParams)
                % get this behavioral variable's values
                % leg or FicTrac parameter - no leg val for FicTrac
                if isempty(behParams1.legs{j})
                    % interpolate FicTrac values to leg timing
                    thisVal = interp1(fictracProc.t, ...
                        fictracProc.(behParams1.whichParams{j}), ...
                        legStepsCont.t, 'spline');
                else
                    % get appropriate leg
                    thisLegInd = legStepsCont.legIDs.ind(...
                        strcmpi(legStepsCont.legIDs.names, ...
                        behParams1.legs{j}));

                    thisVal = legStepsCont.(behParams1.whichParams{j})(:,thisLegInd);
                end

                % add to running tracker for this trial
                thisBehVars1 = cat(2,thisBehVars1,thisVal);
            end

        % if single behavioral parameter to correlate with
        else
            % FicTrac param
            if isempty(behParams1.legs)
                % interpolate FicTrac values to leg timing
                thisBehVars1 = interp1(fictracProc.t, ...
                    fictracProc.(behParams1.whichParams), legStepsCont.t, ...
                    'spline');
            else
                % get appropriate leg
                thisLegInd = legStepsCont.legIDs.ind(...
                    strcmpi(legStepsCont.legIDs.names, legs));

                thisBehVars1 = legStepsCont.(behParams1.whichParams)(:,thisLegInd);
            end
        end

        % get behavioral variable 2
        % if multiple behavioral parameters to correlate ephys with
        if iscell(behParams2.whichParams)
            % preallocate tracker for this trial's variables
            thisBehVars2 = [];

            % loop through all behavioral parameters
            for j = 1:length(behParams2.whichParams)
                % get this behavioral variable's values
                % leg or FicTrac parameter - no leg val for FicTrac
                if isempty(behParams2.legs{j})
                    % interpolate FicTrac values to leg timing
                    thisVal = interp1(fictracProc.t, ...
                        fictracProc.(behParams2.whichParams{j}), ...
                        legStepsCont.t, 'spline');
                else
                    % get appropriate leg
                    thisLegInd = legStepsCont.legIDs.ind(...
                        strcmpi(legStepsCont.legIDs.names, ...
                        behParams2.legs{j}));

                    thisVal = legStepsCont.(behParams2.whichParams{j})(:,thisLegInd);
                end

                % add to running tracker for this trial
                thisBehVars2 = cat(2,thisBehVars2,thisVal);
            end

        % if single behavioral parameter to correlate with
        else
            % FicTrac param
            if isempty(behParams2.legs)
                % interpolate FicTrac values to leg timing
                thisBehVars2 = interp1(fictracProc.t, ...
                    fictracProc.(behParams2.whichParams), legStepsCont.t, ...
                    'spline');
            else
                % get appropriate leg
                thisLegInd = legStepsCont.legIDs.ind(...
                    strcmpi(legStepsCont.legIDs.names, legs));

                thisBehVars2 = legStepsCont.(behParams2.whichParams)(:,thisLegInd);
            end
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
                thisBehVars1(thisNotMoveLog,:) = nan;
                thisBehVars2(thisNotMoveLog,:) = nan;
            else
                thisBehVars1(thisNotMoveLog) = nan;
                thisBehVars2(thisNotMoveLog) = nan;
            end
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
                    thisBehVars1(thisIInjLog,:) = nan;
                    thisBehVars2(thisIInjLog,:) = nan;
                else
                    thisBehVars1(thisIInjLog) = nan;
                    thisBehVars2(thisIInjLog) = nan;
                end         
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
            thisBehVars1(condFalseLog,:) = nan;
            thisBehVars2(condFalseLog,:) = nan;
        end



        % remove NaNs
        rmInd = [];
        for j = 1:size(thisBehVars1,1)
            % get all sample indices to remove
            if (any(isnan(thisBehVars1(j,:))) || ...
                    any(isnan(thisBehVars2(j,:))))
                rmInd = [rmInd; j];
            end
        end

        thisBehVars1(rmInd,:) = [];
        thisBehVars2(rmInd,:) = [];


        % add to variable tracker across pData files
        behVals1 = cat(1,behVals1, thisBehVars1);
        behVals2 = cat(1,behVals2, thisBehVars2);
    end


    % remove outliers
    [~,beh1Outliers] = rmoutliers(behVals1,"percentiles",[0.1 99.9], 1);
    [~,beh2Outliers] = rmoutliers(behVals2,"percentiles",[0.1 99.9], 1);

    outRmLog = beh1Outliers | beh2Outliers;

    behVals1(outRmLog,:) = [];
    behVals2(outRmLog,:) = [];



    % if multiple behavioral variables, get projection to equal weighting
    % for behParams1
    if iscell(behParams1.whichParams)
        % handle circular parameter correctly
        % assume all behParams are circular
        if any(strcmpi(behParams1.whichParams, 'stepDirection'))
            numParams = length(behParams1.whichParams);

            meanSubBehVals = behVals1;
            for i = 1:numParams
                thisBehVals = behVals1(:,i);
                meanSubBehVals(:,i) = thisBehVals - ...
                    rad2deg(circ_mean(deg2rad(thisBehVals)));
            end

            behVals11D = wrapTo180(rad2deg(circ_mean(deg2rad(meanSubBehVals), [], 2)));

        % not circular parameter, project    
        else
            % get coefficients
            coeffs = ones(length(behParams1.whichParams), 1) * ...
                (1/length(behParams1.whichParams));
            behVals11D = getLinProj(behVals1, coeffs);
        end
    % otherwise,    
    else
        % just original behavioral value 
%         behVals1D = behVals;

        % deal with circular parameters appropriately
        if strcmpi(behParams1.whichParams,'stepDirection')
            behVals11D = behVals1 - rad2deg(circ_mean(deg2rad(behVals1)));
        else      
            % mean subtracted behavioral value
            behVals11D = behVals1 - mean(behVals1);
        end
    end

    % if multiple behavioral variables, get projection to equal weighting
    % for behParams2
    if iscell(behParams2.whichParams)
        % handle circular parameter correctly
        % assume all behParams are circular
        if any(strcmpi(behParams2.whichParams, 'stepDirection'))
            numParams = length(behParams2.whichParams);

            meanSubBehVals = behVals2;
            for i = 1:numParams
                thisBehVals = behVals2(:,i);
                meanSubBehVals(:,i) = thisBehVals - ...
                    rad2deg(circ_mean(deg2rad(thisBehVals)));
            end

            behVals21D = wrapTo180(rad2deg(circ_mean(deg2rad(meanSubBehVals), [], 2)));

        % not circular parameter, project    
        else
            % get coefficients
            coeffs = ones(length(behParams2.whichParams), 1) * ...
                (1/length(behParams2.whichParams));
            behVals21D = getLinProj(behVals2, coeffs);
        end
    % otherwise,    
    else
        % just original behavioral value 
%         behVals1D = behVals;

        % deal with circular parameters appropriately
        if strcmpi(behParams2.whichParams,'stepDirection')
            behVals21D = behVals2 - rad2deg(circ_mean(deg2rad(behVals2)));
        else      
            % mean subtracted behavioral value
            behVals21D = behVals2 - mean(behVals2);
        end
    end

    % get correlation between two 1D behavioral parameters
    r = corr(behVals11D, behVals21D);
end