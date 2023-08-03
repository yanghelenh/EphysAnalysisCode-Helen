% savePCAinterpStepParamEphys.m
%
% Function that takes in one or two sets of step parameters, performs PCA 
%  on each, and compares the two against each other.
% Unlike pcaCompStepParams(), this function does not operate on turn bouts.
%  Instead, it works on step parameters interpolated to the specified frame
%  rate.
% Select pData files through GUI
% Saves output file with value of PC1 for each set for all time points.
%  Also saves values of all other step parameters, FicTrac velocity, and 
%  ephys data values at these same time points. Ephys data is at specified
%  delay.
% Removes not moving, FicTrac dropInd, stimulation bouts from consideration
%
% INPUTS:
%   set1 - struct of first set of step params
%     name - string for name of this set
%     params - cell array of step param names, 1 per variable
%     legs - cell array of leg names (R1-3, L1-3), 1 per variable
%     whichPhase - which phase ('swing' or 'stance')
%   set2 - struct of second set of step params or [] if no 2nd set
%     name - string for name of this set
%     params - cell array of step param names, 1 per variable
%     legs - cell array of leg names (R1-3, L1-3), 1 per variable
%     whichPhase - which phase ('swing' or 'stance')
%   interpFrameRate - sampling rate to interpolate step parameters to, in
%       Hz
%   ephysDelay - time in seconds to shift matched ephys data, negative for
%       ephys b/f behavior
%   postStimExclDur - additional time, in sec, after stimulation (opto or 
%       iInj to exclude turning bouts from overlapping with
%   pDataPath - full path to pData directory
%   saveFileName - name of output file
%   saveFileDir - full path to directory to save output file
%
% OUTPUTS:
%   none, but saves output file
%
% CREATED: 7/31/23 - HHY
%
% UPDATED:
%   7/31/23 - HHY
%
function savePCAinterpStepParamEphys(set1, set2, interpFrameRate, ...
    ephysDelay, postStimExclDur, pDataPath, pDataFNames, saveFileName, ...
    saveFileDir)

    % names of all step parameters to save
    stepParamNames = {'stepLengths', 'stepXLengths',...
        'stepYLengths', 'stepDirections', 'stepDurations', 'stepSpeeds',...
        'stepVelX', 'stepVelY', 'stepAEPX', 'stepAEPY', 'stepPEPX', ...
        'stepPEPY'};

    % fictrac parameters to save
    fictracParamNames = {'fwdVel', 'yawAngVel', 'slideVel'};

    % ephys parameters to save
    ephysParamNames = {'spikeRate', 'medFiltV'};

    % number of step parameter variables
    set1NumVars = length(set1.params);
    if ~isempty(set2)
        set2NumVars = length(set2.params);
    end

    % interpolation ifi
    ifi = 1/interpFrameRate;

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

    % preallocate outputs
    pDataFiles.names = pDataFNames;
    pDataFiles.inds = [];
    rmvPDataInd = [];
    currNumPts = 0;

    % matrices on which to run PCA; each column is variable, each row is a
    %  time point; pooled across pData files
    set1PCVars = [];
    set2PCVars = [];

    % for step parameters at matched data points to PC scores
    for i = 1:length(stepParamNames)
        matchStanceStepParams.(stepParamNames{i}) = [];
        matchSwingStepParams.(stepParamNames{i}) = [];
    end

    % for FicTrac parameters at matched data points
    for i = 1:length(fictracParamNames)
        matchFictracParams.(fictracParamNames{i}) = [];
    end

    % for ephys parameters at matched data points
    for i = 1:length(ephysParamNames)
        matchEphysParams.(ephysParamNames{i}) = [];
    end

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
        if (~any(strcmpi(pDatVarsNames, 'legSteps')) || ...
                ~any(strcmpi(pDatVarsNames, 'fictracProc')) || ...
                ~any(strcmpi(pDatVarsNames, 'moveNotMove')) || ...
                ~any(strcmpi(pDatVarsNames, 'ephysSpikes')))
            rmvPDataInd = [rmvPDataInd; i];
            continue;
        end

        % load variables from pData
        if (any(strcmpi(pDatVarsNames, 'iInj')))
            load(pDataFullPath, 'legTrack', 'moveNotMove', 'fictracProc', ...
                'legSteps', 'ephysSpikes', 'iInj');
        else
            load(pDataFullPath, 'legTrack', 'moveNotMove', 'fictracProc', ...
                'legSteps', 'ephysSpikes');
        end

        % legIDs
        legIDs.ind = legSteps.legIDs.ind;
        legIDs.name = legSteps.legIDs.names;

        % get sample times to interpolate to
        startTime = max([legTrack.t(1), fictracProc.t(1)]);
        endTime = min([legTrack.t(end), fictracProc.t(end)]);

        sampTime = startTime:ifi:endTime;
        sampTime = sampTime';

        % check that there are at least 2 steps for each leg for this pData
        %  file
        tooFewSteps = false;
        for j = 1:length(legIDs.ind)
            if (sum(legSteps.stepWhichLeg == legIDs.ind(j)) < 2)
                tooFewSteps = true;
            end
        end

        % if not, skip to next pData file
        if (tooFewSteps)
            rmvPDataInd = [rmvPDataInd; i];
            continue;
        end

        % preallocate matrix of var vals for PCA
        if ~isempty(set2)
            pcaVarsMat = zeros(length(sampTime), set1NumVars + set2NumVars);
        else
            pcaVarsMat = zeros(length(sampTime), set1NumVars);
        end
        
        % loop through all set 1 variables
        for j = 1:set1NumVars
            pcaVarsMat(:,j) = stepParam2VectorSpline(legSteps, ...
                legTrack.t, sampTime, ...
                moveNotMove, set1.legs{j}, set1.params{j}, ...
                set1.whichPhase{j});
        end

        % loop through all set 2 variables, if needed
        if ~isempty(set2)
            for j = 1:set2NumVars
                pcaVarsMat(:,j+set1NumVars) = stepParam2VectorSpline(...
                    legSteps, legTrack.t, sampTime, moveNotMove, ...
                    set2.legs{j}, set2.params{j}, set2.whichPhase{j});
            end
        end

        % get interpolated values for all step parameters
        for j = 1:length(stepParamNames)
            % reinialize for each step parameter
            thisStanceStepVal = zeros(length(sampTime), length(legIDs.ind));
            thisSwingStepVal = zeros(length(sampTime), length(legIDs.ind));

            % loop through all legs
            for k = 1:length(legIDs.ind)
                thisStanceStepVal (:,k) = ...
                stepParam2VectorSpline(legSteps, legTrack.t, sampTime, ...
                moveNotMove, legIDs.name{k}, stepParamNames{j}, ...
                'stance');

                thisSwingStepVal (:,k) = ...
                stepParam2VectorSpline(legSteps, legTrack.t, sampTime, ...
                moveNotMove, legIDs.name{k}, stepParamNames{j}, ...
                'swing');
            end

            % structs for this pData
            thisStanceMatchStepParams.(stepParamNames{j}) = thisStanceStepVal;
            thisSwingMatchStepParams.(stepParamNames{j}) = thisSwingStepVal;
        end

        % get interpolated values for all FicTrac parameters
        for j = 1:length(fictracParamNames)
            thisMatchFictracParams.(fictracParamNames{j}) = ...
                interp1(fictracProc.t, ...
                fictracProc.(fictracParamNames{j}), sampTime, 'spline');
        end

        % get interpolated values for all ephys parameters, with delay
        for j = 1:length(ephysParamNames)
            thisMatchEphysParams.(ephysParamNames{j}) = ...
                interp1(ephysSpikes.t - ephysDelay, ...
                ephysSpikes.(ephysParamNames{j}), sampTime, 'spline', ...
                'extrap');
        end

        % find stimulation times, flag for removal
        allStimRmvInd = [];

        if (any(strcmpi(pDatVarsNames, 'iInj')))
            for j = 1:length(iInj.startTimes)
                thisStartTime = iInj.startTimes(j);
                thisEndTime = iInj.endTimes(j) + postStimExclDur;
    
                thisRmvInd = find((sampTime >= thisStartTime) & ...
                    (sampTime <= thisEndTime));
    
                allStimRmvInd = [allStimRmvInd; thisRmvInd];
            end
        end


        % find all NaNs in PCA variables, note which rows those are found
        % in
        nanRows = [];

        for j = 1:size(pcaVarsMat,1)
            if (any(isnan(pcaVarsMat(j,:))))
                nanRows = [nanRows; j];
            end
        end

%         % remove outliers
%         [~, rmvOutLog] = rmoutliers(pcaVarsMat);

        % merge stim and NaN removal indices
%         allRmvInd = unique([nanRows; allStimRmvInd; find(rmvOutLog)]);
        allRmvInd = unique([nanRows; allStimRmvInd]);

        % remove rows containing NaNs from PCA matrix as well as all
        %  interpolated variables
        pcaVarsMat(allRmvInd,:) = [];

        for j = 1:length(stepParamNames)
            thisStanceMatchStepParams.(stepParamNames{j})(allRmvInd,:) = [];
            thisSwingMatchStepParams.(stepParamNames{j})(allRmvInd,:) = [];
        end

        for j = 1:length(fictracParamNames)
            thisMatchFictracParams.(fictracParamNames{j})(allRmvInd,:) = [];
        end

        for j = 1:length(ephysParamNames)
            thisMatchEphysParams.(ephysParamNames{j})(allRmvInd,:) = [];
        end

        
        % add these values to running tracker across pData files
        set1PCVars = [set1PCVars; pcaVarsMat(:,1:set1NumVars)];
        if ~isempty(set2)
            set2PCVars = [set2PCVars; pcaVarsMat(:, ...
                (set1NumVars+1):(set1NumVars + set2NumVars))];
        end

        for j = 1:length(stepParamNames)
            matchStanceStepParams.(stepParamNames{j}) = ...
                [matchStanceStepParams.(stepParamNames{j}); ...
                thisStanceMatchStepParams.(stepParamNames{j})];

            matchSwingStepParams.(stepParamNames{j}) = ...
                [matchSwingStepParams.(stepParamNames{j}); ...
                thisSwingMatchStepParams.(stepParamNames{j})];
        end

        for j = 1:length(fictracParamNames)
            matchFictracParams.(fictracParamNames{j}) = ...
                [matchFictracParams.(fictracParamNames{j}); ...
                thisMatchFictracParams.(fictracParamNames{j})];
        end

        for j = 1:length(ephysParamNames)
            matchEphysParams.(ephysParamNames{j}) = ...
                [matchEphysParams.(ephysParamNames{j}); ...
                thisMatchEphysParams.(ephysParamNames{j})];
        end

        % update pData tracker with start and end indices for this pData
        %  file
        pDataFiles.inds = [pDataFiles.inds; ...
            currNumPts + 1, currNumPts + size(pcaVarsMat,1)];
        currNumPts = currNumPts + size(pcaVarsMat,1);

    end

    % remove pData files that didn't contribute any pts
    pDataFiles.names(rmvPDataInd) = [];

    % remove outliers
    % NOTE % this current screws up pDataFiles tracking. Temp fix is just
    %  to record which ones are removed
    [~, rmvOutLog1] = rmoutliers(set1PCVars);
    if ~isempty(set2)
        [~, rmvOutLog2] = rmoutliers(set2PCVars);
    end
    rmvOutLog = rmvOutLog1 | rmvOutLog2;
    
    % remove outliers from PC vars
    set1PCVars(rmvOutLog,:) = [];
    if ~isempty(set2)
        set2PCVars(rmvOutLog,:) = [];
    end
    outlierRmvInd = find(rmvOutLog);

    % remove outliers from matched vars
    for j = 1:length(stepParamNames)
        matchStanceStepParams.(stepParamNames{j})(outlierRmvInd,:) = [];
        matchSwingStepParams.(stepParamNames{j})(outlierRmvInd,:) = [];
    end

    for j = 1:length(fictracParamNames)
        matchFictracParams.(fictracParamNames{j})(outlierRmvInd,:) = [];
    end

    for j = 1:length(ephysParamNames)
        matchEphysParams.(ephysParamNames{j})(outlierRmvInd,:) = [];
    end

    % perform PCA - across all pData files
    [set1.coeff, set1Score, set1.latent, set1.tsquared, ...
        set1.explained, set1.mu] = pca(set1PCVars);

    if ~isempty(set2)
        [set2.coeff, set2Score, set2.latent, set2.tsquared, ...
            set2.explained, set2.mu] = pca(set2PCVars);
    end

    % save output
    saveFileFullPath = [saveFileDir filesep saveFileName '.mat'];

    if ~isempty(set2)
        save(saveFileFullPath, 'set1', 'set2', 'set1PCVars', 'set2PCVars',...
            'set1Score', 'set2Score', 'matchStanceStepParams', ...
            'matchSwingStepParams', 'matchFictracParams', ...
            'matchEphysParams', 'outlierRmvInd', 'ephysDelay', ...
            'postStimExclDur', ...
            'interpFrameRate', 'legIDs', 'pDataFiles', '-v7.3');
    else
        save(saveFileFullPath, 'set1', 'set1PCVars',...
            'set1Score', 'set2Score', 'matchStanceStepParams', ...
            'matchSwingStepParams', 'matchFictracParams', ...
            'matchEphysParams', 'outlierRmvInd', 'ephysDelay', ...
            'postStimExclDur', ...
            'interpFrameRate', 'legIDs', 'pDataFiles', '-v7.3');
    end

    % print some outputs to screen
    fprintf('Set 1, PC1, variance explained = %.2f%%\n', set1.explained(1));
    if ~isempty(set2)
        fprintf('Set 2, PC1, variance explained = %.2f%%\n', ...
            set2.explained(1));
    end
end