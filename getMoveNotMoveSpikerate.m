% getMoveNotMoveSpikerate.m
%
% Function to get mean spike rate during moving and not moving periods,
%  with transitions excluded.
% Operates on pData files from multiple flies (merges data from pData files
%  from same fly but separates data from different flies)
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
%   transExclDur - additional time, in sec, before and after not
%     transition to exclude from consideration. Length 2 vector for before
%     and after, respectively
%   postStimExclDur - additional time, in sec, after iInj stimulation 
%     to exclude from consideration
%   pDataFNames - cell array of pData file names or [] if select through
%       GUI
%   pDataPath - path to folder containing pData files
%   saveFileName - name to save output file, without .mat
%   saveFileDir - full path to folder in which to save output file
%
% OUTPUTS:
%   none, but saves file
%
% CREATED: 9/1/23 - HHY
%
% UPDATED:
%   9/1/23 - HHY
%   9/29/23 - HHY - add ability to condition on step and fictrac parameters
%
function getMoveNotMoveSpikerate(cond, transExclDur, postStimExclDur, ...
    pDataFNames, pDataPath, saveFileName, saveFileDir)

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
    % to track mean spike rate, 1 for each fly
    % move as column 1, not move as column 2
    meanSpikerate = []; 
    % to track fly names
    allFlyNames = {};

    % initialize
    prevFlyName = ''; % previous fly's name

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

        % save fly name as first pDataName's date, fly (12 characters)
        thisFlyName = pDataName(1:12);

        % check if this pData file has legSteps, fictracProc, moveNotMove
        %  structs, if not, skip
        if (~any(strcmpi(pDatVarsNames, 'moveNotMove')) || ...
                ~any(strcmpi(pDatVarsNames, 'ephysSpikes')))
            continue;
        end

        % new fly
        if ~(strcmp(thisFlyName, prevFlyName))
            % if we're starting new fly and this isn't the first pData file
            % save info for previous fly
            if (i~=1)
                % remove NaNs from tracking vectors
                moveSpikerate = moveSpikerate(~isnan(moveSpikerate));
                notMoveSpikerate = notMoveSpikerate(~isnan(notMoveSpikerate));

                % get mean and add to tracking 
                thisMoveMean = mean(moveSpikerate);
                thisNotMoveMean = mean(notMoveSpikerate);

                meanSpikerate = [meanSpikerate; ...
                    [thisMoveMean, thisNotMoveMean]];

                % save fly name
                allFlyNames = [allFlyNames; prevFlyName];
            end

            % initialize/re-initialize vectors for keeping track of spike
            %  rate across pData files
            moveSpikerate = [];
            notMoveSpikerate = [];
        end

        % load variables from pData
        if (any(strcmpi(pDatVarsNames, 'iInj')))
            load(pDataFullPath, 'moveNotMove', 'ephysSpikes', 'iInj',...
                'legSteps', 'legStepsCont', 'fictracProc');
        else
            load(pDataFullPath, 'moveNotMove', 'ephysSpikes', ...
                'legSteps', 'legStepsCont', 'fictracProc');
        end

        % downsample spike rate
        thisT = downsample(ephysSpikes.t, DSF);
        
        thisMoveEphysVal = interp1(ephysSpikes.t, ephysSpikes.spikeRate, ...
            thisT, 'linear');
        thisNotMoveEphysVal = thisMoveEphysVal;

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
                thisMoveEphysVal(thisIInjLog) = nan;
                thisNotMoveEphysVal(thisIInjLog) = nan;
            end
        end

        % set points for moving ephys that don't meet conditions to NaN
        if ~isempty(cond)
            % initialize logical to track which points meet condition
            % loop through all conditions
            condLog = true(size(thisT));
            for j = 1:length(cond.whichParam)
                % check which data structures conditions belong to
                % check first for stepFwdBool, as it's not a field of any
                %  struct
                if (strcmpi(cond.whichParam{j}, 'stepFwdBool'))
                    thisFwdLog = getStepFwdLogical(legSteps, thisT,...
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
            thisMoveEphysVal(condFalseLog) = nan;
        end

        
        % loop through all moving bouts, add spike rate to tracking vector
        for j = 1:size(moveNotMove.ftMoveBout,1)
            % start and end time of moving bout, with adjustment to exclude
            %  extra time around transitions
            thisStartTime = moveNotMove.ftT(moveNotMove.ftMoveBout(j,1)) + ...
                transExclDur(2);
            thisEndTime = moveNotMove.ftT(moveNotMove.ftMoveBout(j,2)) + ...
                transExclDur(1);

            % get logical into spike rate vector, portion during this 
            %  moving bout
            thisBoutLog = (thisT >= thisStartTime) & (thisT <= thisEndTime);

            % if there are any valid points, grab and add to tracking
            %  vector
            if any(thisBoutLog)
                thisSpikerate = thisMoveEphysVal(thisBoutLog);
                moveSpikerate = [moveSpikerate; thisSpikerate];
            end
        end

        % loop through all not moving bouts, add spike rate to tracking vector
        for j = 1:size(moveNotMove.ftNotMoveBout,1)
            % start and end time of not moving bout, with adjustment to 
            %  exclude extra time around transitions
            thisStartTime = moveNotMove.ftT(moveNotMove.ftNotMoveBout(j,1)) + ...
                transExclDur(2);
            thisEndTime = moveNotMove.ftT(moveNotMove.ftNotMoveBout(j,2)) + ...
                transExclDur(1);

            % get logical into spike rate vector, portion during this 
            %  moving bout
            thisBoutLog = (thisT >= thisStartTime) & (thisT <= thisEndTime);

            % if there are any valid points, grab and add to tracking
            %  vector
            if any(thisBoutLog)
                thisSpikerate = thisNotMoveEphysVal(thisBoutLog);
                notMoveSpikerate = [notMoveSpikerate; thisSpikerate];
            end
        end

        % update fly name
        prevFlyName = thisFlyName;
    end

    % add last fly
    % remove NaNs from tracking vectors
    moveSpikerate = moveSpikerate(~isnan(moveSpikerate));
    notMoveSpikerate = notMoveSpikerate(~isnan(notMoveSpikerate));
    
    % get mean and add to tracking 
    thisMoveMean = mean(moveSpikerate);
    thisNotMoveMean = mean(notMoveSpikerate);
    
    meanSpikerate = [meanSpikerate; ...
        [thisMoveMean, thisNotMoveMean]];

    % save fly name
    allFlyNames = [allFlyNames; prevFlyName];

    % save output file
    saveFileFullName = [saveFileDir filesep saveFileName '.mat'];
    save(saveFileFullName, 'meanSpikerate', 'allFlyNames',...
        'transExclDur', 'postStimExclDur', 'cond', '-v7.3');
end