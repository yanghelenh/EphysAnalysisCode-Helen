% getMoveNotMoveSpikerate.m
%
% Function to get mean spike rate during moving and not moving periods,
%  with transitions excluded.
% Operates on pData files from multiple flies (merges data from pData files
%  from same fly but separates data from different flies)
%
% INPUTS:
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
%
function getMoveNotMoveSpikerate(transExclDur, postStimExclDur, ...
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
            load(pDataFullPath, 'moveNotMove', 'ephysSpikes', 'iInj');
        else
            load(pDataFullPath, 'moveNotMove', 'ephysSpikes');
        end

        % downsample spike rate
        thisT = downsample(ephysSpikes.t, DSF);
        
        thisEphysVal = interp1(ephysSpikes.t, ephysSpikes.spikeRate, ...
            thisT, 'linear');

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
                thisEphysVal(thisIInjLog) = nan;
            end
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
                thisSpikerate = thisEphysVal(thisBoutLog);
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
                thisSpikerate = thisEphysVal(thisBoutLog);
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
        'transExclDur', 'postStimExclDur', '-v7.3');
end