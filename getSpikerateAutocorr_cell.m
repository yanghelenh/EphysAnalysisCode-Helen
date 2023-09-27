% getSpikerateAutocorr_cell.m
%
% Function that computes the autocorrelation of the spike rate for a single
%  cell, only during walking and only under additional conditions on
%  behavioral parameters, specified by cond.
% Across trials within the same fly, weight by number of valid points
%  during that trial
% Select input pData files for single cell through GUI or as input
%
% INPUTS:
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
% CREATED: 8/31/23 - HHY
%
% UPDATED:
%   8/31/23 - HHY
%
function getSpikerateAutocorr_cell(cond, maxLag, notMoveExclDur, ...
    postStimExclDur, pDataFNames, pDataPath, saveFileDir)

    DSF = 20; % downsampling factor

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

        % downsample spike rate
        newT = downsample(ephysSpikes.t, DSF);
        
        thisEphysVal = interp1(ephysSpikes.t, ephysSpikes.spikeRate, ...
            newT, 'linear');

%         get spike rate, initialize with everything, NaN later
%         thisEphysVal = ephysSpikes.spikeRate;

        % get not moving times, set those to NaN
        notMoveStartTimes = moveNotMove.ftT(moveNotMove.ftNotMoveBout(:,1)) - ...
            notMoveExclDur(1);
        notMoveEndTimes = moveNotMove.ftT(moveNotMove.ftNotMoveBout(:,2)) + ...
            notMoveExclDur(2);


        % loop through all bouts, set appropriate inds to NaN
        for j = 1:length(notMoveStartTimes)
            thisNotMoveLog = (newT > notMoveStartTimes(j)) & ...
                (newT < notMoveEndTimes(j));
            % set to NaN
            thisEphysVal(thisNotMoveLog) = nan;
        end


        % set current injection times to NaN
        if (any(strcmpi(pDatVarsNames, 'iInj')))
            % current injection times; only add to end
            iInjStartTimes = iInj.startTimes;
            iInjEndTimes = iInj.endTimes + postStimExclDur;

            % loop through all iInj bouts, set appropriate inds to NaN
            for j = 1:length(iInjStartTimes)
                thisIInjLog = (newT > iInjStartTimes(j)) & ...
                    (newT < iInjEndTimes(j));

                % set to NaN
                thisEphysVal(thisIInjLog) = nan;         
            end
        end

        % set points that don't meet conditions to NaN
        if ~isempty(cond)
            % initialize logical to track which points meet condition
            % loop through all conditions
            condLog = true(size(newT));
            for j = 1:length(cond.whichParam)
                % check which data structures conditions belong to
                % check first for stepFwdBool, as it's not a field of any
                %  struct
                if (strcmpi(cond.whichParam{j}, 'stepFwdBool'))
                    thisFwdLog = getStepFwdLogical(legSteps, newT,...
                        cond.legs{j});
                    if ~(eval(cond.cond{j})) % if target is false, invert
                        thisLog = ~thisFwdLog;
                    else
                        thisLog = thisFwdLog;
                    end
                % if FicTrac    
                elseif isfield(fictracProc,cond.whichParam{j})
                    thisFT = interp1(fictracProc.t, ...
                        fictracProc.(cond.whichParam{j}), newT, ...
                        'spline');
                    thisLog = eval(['thisFT' cond.cond{j}]); 
                % if step parameter    
                elseif isfield(legStepsCont,cond.whichParam{j})
                    thisLegInd = legSteps.legIDs.ind(strcmpi(legSteps.legIDs.names, ...
                        cond.legs{j}));
                    thisStep = interp1(legStepsCont.t, ...
                        legStepsCont.(cond.whichParam{j})(:,thisLegInd), ...
                        newT, 'spline');
                    thisLog = eval(['thisStep' cond.cond{j}]);
                end
                % update logical for all conditions
                condLog = condLog & thisLog;
            end
            % invert condition logical
            condFalseLog = ~condLog;
            % set every time when condition logcial false to nan
            thisEphysVal(condFalseLog) = nan;
        end


        % get autocorrelation

        % get number of lags
        ifi = median(diff(newT));
        % number of lags to one side
        numLags1Side = ceil(maxLag / ifi);
        % number of lags, total
        numLags = numLags1Side * 2 + 1;

        % autocorrelation, with NaNs allowed
        [allCorr, lags] = xcorrWGaps(thisEphysVal, thisEphysVal, numLags);

        % number of points this trial contributes
        trialCount(i) = length(thisEphysVal(~isnan(thisEphysVal)));

        % add autocorr to running tracker across trials
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
    autoCorr = sum(weightAutoCorrAll, 2);

    % convert lags in samples to t
    lagsT = lags * ifi;

    % save data to output file
    saveFileFullName = [saveFileDir filesep flyName '_autocorr.mat'];
    save(saveFileFullName, 'autoCorr', 'lags', 'lagsT', 'cond', 'maxLag', ...
        'notMoveExclDur', 'postStimExclDur', '-v7.3');
end