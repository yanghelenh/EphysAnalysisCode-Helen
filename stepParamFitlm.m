% stepParamFitlm.m
%
% Function that takes in one or more pData files from same fly (specified
%  through GUI), a list of step parameters to use as predictors, a time
%  delay, and a sample rate and uses the MATLAB fitlm function to produce a
%  multiple linear regression model that predicts spike rate or membrane
%  voltage from those step parameters. Saves model output plus info about
%  input parameters into date_fly_cell_stepFitlm.mat file
%
% INPUTS:
%   paramList - cell array of names of step parameters to use as
%       predictors, matches the names of fields of legSteps. Or single
%       string
%   legList - cell array of names of legs to consider for each predictor.
%       Or single string
%   phaseList - cell array of phases (swing, stance) to consider. Or single
%       string
%   ephysVarName - name of ephys variable to predict: spikeRate or medFiltV
%   timeDelay - offset between times of step parameters and ephys, in sec;
%       negative values for ephys before behavior, positive for ephys after
%       behavior
%   sampRate - sample rate at which to compute model, in Hz
%   savePath - full path to folder in which to save output, extract name of
%       fly from first pData file
%
% OUTPUTS:
%   none, but saves to output file
%   
% CREATED: 7/11/22 - HHY
%
% UPDATED:
%   7/11/22 - HHY
%   7/13/22 - HHY - add variable names to fitlm call
%   7/27/22 - HHY - step parameter vector uses spline interpolation instead
%       of constant value for whole step
%
function stepParamFitlm(paramList, legList, phaseList, ephysVarName, ...
    timeDelay, sampRate, savePath)

    % prompt user to select pData files
    disp('Select pData files');
    [pDataFNames, pDataPath] = uigetfile('*.mat', 'Select pData files', ...
        pDataDir(), 'MultiSelect', 'on');
    
    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end

    % initialize variables for predictors and ephysVar
    % to allow concatenation across pData files
    predMatrix = [];
    ephysVar = [];
    
    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end
        
        pDataFullPath = [pDataPath pDataName];

        if (i==1) % get fly name from first pData file
            flyName = pDataName(1:(end-18)); % remove trial, pData, and mat
        end

        % check that pData file has ephysSpikes and legSteps structs
        varNames = who('-file', pDataFullPath);
        if (ismember('ephysSpikes',varNames) && ...
                ismember('legSteps',varNames))
            % load relevant data from pData
            load(pDataFullPath, 'ephysSpikes','legSteps','moveNotMove');
        % if pData file doesn't contain both structs, skip to next file
        else 
            fprintf('%s does not contain both ephysSpikes and legSteps. Skipping to next pData file.\n',...
                pDataName);
            continue;
        end

        % get time vector: 0 to final time point of ephys var
        startTime = 0;
        endTime = ephysSpikes.t(end);
        ifi = 1/sampRate;
        t = startTime:ifi:endTime;
        t = t';

        % adjusted time vector for ephys var, incorporate time delay
        adjTime = ephysSpikes.t - timeDelay;

        % get ephys var values at time t points, using interpolation
        %  (method: spline)
        thisEphysVar = interp1(adjTime, ephysSpikes.(ephysVarName),...
            t,'spline',nan);
        if ~iscolumn(thisEphysVar) % make sure it's a column vector
            thisEphysVar = thisEphysVar';
        end

        % append this ephys var values to vector across pData
        ephysVar = [ephysVar; thisEphysVar];

        % get number of step parameters
        if iscell(paramList) % cell array
            numStepParams = length(paramList);
        else % only 1 string
            numStepParams = 1;
        end

        % get number of legs
        if iscell(legList) % cell array
            numLegs = length(legList);
        else % only one string
            numLegs = 1;
        end

        % get number of phases
        if iscell(phaseList) % cell array
            numPhases = length(phaseList);
        else
            numPhases = 1;
        end

        % initialize matrix for this pData file
        % total number of predictors, when considering all parameters,
        %  legs, phases separately
        totalNumPred = numStepParams * numLegs * numPhases;

        % initialize cell array for all variable names (+1 for response
        %  var: ephysVar)
        allVarNames = cell(1,totalNumPred + 1);

        % each column of matrix is predictor, each row is time point of t
        thisPredMatrix = nan(length(t),totalNumPred);

        % counter for which predictor we're on
        counter = 1;

        % loop through all step parameters
        for j = 1:numStepParams
            % handle list of parameters vs. 1 parameter
            if (numStepParams == 1) % only 1 string
                thisStepParam = paramList;
            else
                thisStepParam = paramList{j};
            end

            % loop through all legs
            for k = 1:numLegs
                % handle list of legs vs. 1 leg only
                if (numLegs == 1) % only 1 leg
                    thisLeg = legList;
                else
                    thisLeg = legList{k};
                end

                for l = 1:numPhases
                    % handle list of phases vs. 1  phase only
                    if (numPhases == 1) % only one phase
                        thisPhase = phaseList;
                    else
                        thisPhase = phaseList{l};
                    end

                    thisPredName = sprintf('%s_%s_%s',thisStepParam,...
                        thisLeg, thisPhase);

                    % get param value for this param, leg, phase
                    thisParamVal = stepParam2Vector(legSteps, t, ...
                        thisLeg, thisStepParam, thisPhase);
%                     thisParamVal = stepParam2VectorSpline(legSteps, t, ...
%                         moveNotMove, thisLeg, thisStepParam, thisPhase);
                    
                    % add this value to predictor matrix
                    thisPredMatrix(:,counter) = thisParamVal;

                    % add this name to predictor name list
                    allVarNames{counter} = thisPredName;

                    counter = counter + 1; % increment
                end
            end
        end

        % append this predictor matrix to one for all pData
        predMatrix = [predMatrix; thisPredMatrix];
    end

    allVarNames{end} = ephysVarName; % add ephys variable name

    % perform multiple linear regression, returns LinearModel object
    mdl = fitlm(predMatrix, ephysVar, 'VarNames',allVarNames);

    % put all inputs into struct
    mdlInputs.paramList = paramList;
    mdlInputs.legList = legList;
    mdlInputs.phaseList = phaseList;
    mdlInputs.ephysVarName = ephysVarName;
    mdlInputs.timeDelay = timeDelay;
    mdlInputs.sampRate = sampRate;

    % save file full path
    saveFileFullPath = [savePath filesep flyName '_stepFitlm.mat'];

    % save to output file
    save(saveFileFullPath, 'mdl', 'mdlInputs', '-v7.3');
end