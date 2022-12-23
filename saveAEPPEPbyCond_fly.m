% saveAEPPEPbyCond_fly.m
%
% Function that takes in user selected pData files (for 1 fly) on which
%  sortLegStepsByIInj or sortLegStepsByOpto has been run and saves mean, 
%  and std dev for AEPs and PEPs for each condition for that fly
%
% INPUTS:
%   stim - struct of parameters about which stimulation conditions to plot,
%     if any, with fields:
%       whichStim - which type of stimulation to condition on, specified as
%           string 'iInj' or 'opto'. [] for no conditioning
%       durs - stim durations to include, for iInj or opto
%       NDs or amps - specify stimulus amplitude, NDs for opto, amps for
%           iInj
%   whichPhase - which phase of leg movement, specified as string 'swing'
%     or 'stance'
%   flipLegsLR - logical for whether to flip leg assignments over midline
%   saveDir - folder in which to save file, will be named date_fly_cell
%       after 1st pData file
%
% OUTPUTS:
%   none, but generates saved file
%
% CREATED: 12/2/22 - HHY
%
% UPDATED:
%   12/2/22 - HHY
%   12/5/22 - HHY - add SEM
%
function saveAEPPEPbyCond_fly(stim, whichPhase, flipLegsLR, saveDir)

    % if swap right and left sides
    if (flipLegsLR)
        legInd = [4 5 6 1 2 3]; % reassign indices
    else
        legInd = 1:6; % indicies into raw position matricies for legs
    end


    NUM_LEGS = 6; % number of legs

    % prompt user to select pData files
    [pDataFNames, pDataPath] = uigetfile('*.mat', 'Select pData files', ...
        pDataDir(), 'MultiSelect', 'on');
    
    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end

    % initialize cell array for keeping track of steps in each category
    % generate key for mapping b/w cell array indices and specific durs and
    %  amps/NDs
    if (isempty(stim.whichStim))
        numCats = 1;
    else
        switch lower(stim.whichStim)
            case 'iinj'
                numCats = length(stim.durs) * length(stim.amps) + 1;
    
                % generate key for mapping b/w indices and amps and durs
                iInjCatAmps = zeros(1,numCats);
                iInjCatDurs = zeros(1,numCats);
                % counter index into vectors, skip 1 for 0,0, no iInj
                counter = 2; 
        
                % assign amps to indices
                for i = 1:length(stim.amps)
                    for j = 1:length(stim.durs)
                        iInjCatAmps(counter) = stim.amps(i);
                        iInjCatDurs(counter) = stim.durs(j);
        
                        counter = counter + 1;
                    end
                end
            case 'opto'
                numCats = length(stim.durs) * length(stim.NDs) + 1;
    
                % generate key for mapping b/w indices and NDs and durs
                optoCatNDs = ones(1,numCats) * -1;
                optoCatDurs = ones(1,numCats) * -1;
                % counter index into vectors, skip 1 for -1, -1, no opto
                counter = 2; 
        
                % assign NDs, durs to indices
                for i = 1:length(stim.NDs)
                    for j = 1:length(stim.durs)
                        optoCatNDs(counter) = stim.NDs(i);
                        optoCatDurs(counter) = stim.durs(j);
        
                        counter = counter + 1;
                    end
                end
        end
    end
    allAEPXvals = cell(numCats,1); % values for X
    allAEPYvals = cell(numCats,1); % values for Y
    allAEPLegs = cell(numCats,1); % which leg the steps belong to

    allPEPXvals = cell(numCats,1); % values for X
    allPEPYvals = cell(numCats,1); % values for Y
    allPEPLegs = cell(numCats,1); % which leg the steps belong to

    % convert whichPhase to index
    switch lower(whichPhase)
        case 'swing'
            phaseInd = -1;
        case 'stance'
            phaseInd = 1;
    end

    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end
        
        pDataFullPath = [pDataPath pDataName];

        % save fly name as first pDataName's date, fly, cell (19 characters)
        if (i == 1)
            flyName = pDataName(1:19);
        end

        % get variables in pData file
        pDataMatObj = matfile(pDataFullPath);
        pDataVarsStrct = whos(pDataMatObj);
        pDataVars = struct2cell(pDataVarsStrct);
        pDataVarsNames = pDataVars(1,:); % cell array of names     

        % for conditioned on current injection
        if (strcmpi(stim.whichStim, 'iInj'))
            % check that pData has appropriate structs; otherwise, skip
            if (any(contains(pDataVarsNames, 'legSteps')) && ...
                    any(contains(pDataVarsNames, 'legStepsByIinj')))

                % load pData
                load(pDataFullPath, 'legSteps', 'legStepsByIinj');

                % AEP and PEP, also selects phase
                stepAEPX = groupStepParamsBySwingStance(...
                    legSteps.stepAEPX, ...
                    legSteps.stepSwingStance, phaseInd);
                stepAEPY = groupStepParamsBySwingStance(...
                    legSteps.stepAEPY, ...
                    legSteps.stepSwingStance, phaseInd);
                stepPEPX = groupStepParamsBySwingStance(...
                    legSteps.stepPEPX, ...
                    legSteps.stepSwingStance, phaseInd);
                stepPEPY = groupStepParamsBySwingStance(...
                    legSteps.stepPEPY, ...
                    legSteps.stepSwingStance, phaseInd);

                % get step categories for this phase
                thisLegStepCat = groupStepParamsBySwingStance(...
                    legStepsByIinj.stepIinjCat, legSteps.stepSwingStance,...
                    phaseInd);

                % loop through all conditions
                for j = 1:length(iInjCatAmps)
                    % AEP
                    % get values for just this condition
                    [stepValX, stepValY, theseWhichLeg] = ...
                        getCatStepValsXY(stepAEPX, stepAEPY, iInjCatAmps(j),...
                        iInjCatDurs(j), legSteps.stepWhichLeg, ...
                        thisLegStepCat, legStepsByIinj.iInjCatAmps, ...
                        legStepsByIinj.iInjCatDurs);
                    % append to values from across pData files
                    allAEPXvals{j} = [allAEPXvals{j}; stepValX];
                    allAEPYvals{j} = [allAEPYvals{j}; stepValY];
                    allAEPLegs{j} = [allAEPLegs{j}; theseWhichLeg];

                    % PEP
                    % get values for just this condition
                    [stepValX, stepValY, theseWhichLeg] = ...
                        getCatStepValsXY(stepPEPX, stepPEPY, iInjCatAmps(j),...
                        iInjCatDurs(j), legSteps.stepWhichLeg, ...
                        thisLegStepCat, legStepsByIinj.iInjCatAmps, ...
                        legStepsByIinj.iInjCatDurs);
                    % append to values from across pData files
                    allPEPXvals{j} = [allPEPXvals{j}; stepValX];
                    allPEPYvals{j} = [allPEPYvals{j}; stepValY];
                    allPEPLegs{j} = [allPEPLegs{j}; theseWhichLeg];
                end
            end
        % for conditioned on opto stim
        elseif (strcmpi(stim.whichStim, 'opto'))
            % check that pData has appropriate structs; otherwise, skip
            if (any(contains(pDataVarsNames, 'legSteps')) && ...
                    any(contains(pDataVarsNames, 'legStepsByOpto')) && ...
                    any(contains(pDataVarsNames, 'opto')))

                % load pData
                load(pDataFullPath, 'legSteps', 'legStepsByOpto', 'opto');

                % AEP PEP, also selects phase, returns single vector
                stepAEPX = groupStepParamsBySwingStance(...
                    legSteps.stepAEPX, ...
                    legSteps.stepSwingStance, phaseInd);
                stepAEPY = groupStepParamsBySwingStance(...
                    legSteps.stepAEPY, ...
                    legSteps.stepSwingStance, phaseInd);
                stepPEPX = groupStepParamsBySwingStance(...
                    legSteps.stepPEPX, ...
                    legSteps.stepSwingStance, phaseInd);
                stepPEPY = groupStepParamsBySwingStance(...
                    legSteps.stepPEPY, ...
                    legSteps.stepSwingStance, phaseInd);

                % get step categories for this phase
                thisLegStepCat = groupStepParamsBySwingStance(...
                    legStepsByOpto.stepOptoCat, legSteps.stepSwingStance,...
                    phaseInd);

                % create ND key for this pData (only 1 ND per trial, match
                %  length of durs key)
                ndKey = repmat(opto.stimParams.ndFilter,...
                    size(legStepsByOpto.optoCat));
                ndKey(1) = -1; % first index, for no stim

                % loop through all conditions
                for j = 1:length(optoCatNDs)
                    % AEP
                    % get values for just this condition
                    [stepValX, stepValY, theseWhichLeg] = ...
                        getCatStepValsXY(stepAEPX, stepAEPY, optoCatNDs(j),...
                        optoCatDurs(j), legSteps.stepWhichLeg, ...
                        thisLegStepCat, ndKey, legStepsByOpto.optoCat);
                    % append to values from across pData files
                    allAEPXvals{j} = [allAEPXvals{j}; stepValX];
                    allAEPYvals{j} = [allAEPYvals{j}; stepValY];
                    allAEPLegs{j} = [allAEPLegs{j}; theseWhichLeg];

                    % PEP
                    % get values for just this condition
                    [stepValX, stepValY, theseWhichLeg] = ...
                        getCatStepValsXY(stepPEPX, stepPEPY, optoCatNDs(j),...
                        optoCatDurs(j), legSteps.stepWhichLeg, ...
                        thisLegStepCat, ndKey, legStepsByOpto.optoCat);
                    % append to values from across pData files
                    allPEPXvals{j} = [allPEPXvals{j}; stepValX];
                    allPEPYvals{j} = [allPEPYvals{j}; stepValY];
                    allPEPLegs{j} = [allPEPLegs{j}; theseWhichLeg];
                end
            end
        % invalid stim specification
        else
            disp('Invalid input for whichStim. Ending');
            return;
        end
    end

    % AEP
    % preallocate, matrices for means, stddev
    AEPxMeans = zeros(length(allAEPXvals), NUM_LEGS);
    AEPyMeans = zeros(length(allAEPYvals), NUM_LEGS);
    AEPxStdDev = zeros(length(allAEPXvals), NUM_LEGS);
    AEPyStdDev = zeros(length(allAEPYvals), NUM_LEGS);
    AEPxSEM = zeros(length(allAEPXvals), NUM_LEGS);
    AEPySEM = zeros(length(allAEPYvals), NUM_LEGS);

    % get means and std dev, by leg
    for i = 1:length(allAEPXvals)
        for j = 1:NUM_LEGS
            thisLegX = allAEPXvals{i}(allAEPLegs{i} == legInd(j));
            AEPxMeans(i,j) = mean(thisLegX);
            AEPxStdDev(i,j) = std(thisLegX);
            AEPxSEM(i,j) = std(thisLegX) / sqrt(length(thisLegX));

            thisLegY = allAEPYvals{i}(allAEPLegs{i} == legInd(j));
            AEPyMeans(i,j) = mean(thisLegY);
            AEPyStdDev(i,j) = std(thisLegY);
            AEPySEM(i,j) = std(thisLegY) / sqrt(length(thisLegY));
        end
    end

    % PEP
    % preallocate, matrices for means, stddev
    PEPxMeans = zeros(length(allPEPXvals), NUM_LEGS);
    PEPyMeans = zeros(length(allPEPXvals), NUM_LEGS);
    PEPxStdDev = zeros(length(allPEPXvals), NUM_LEGS);
    PEPyStdDev = zeros(length(allPEPXvals), NUM_LEGS);
    PEPxSEM = zeros(length(allPEPXvals), NUM_LEGS);
    PEPySEM = zeros(length(allPEPXvals), NUM_LEGS);

    % get means and std dev, by leg
    for i = 1:length(allPEPXvals)
        for j = 1:NUM_LEGS
            thisLegX = allPEPXvals{i}(allPEPLegs{i} == legInd(j));
            PEPxMeans(i,j) = mean(thisLegX);
            PEPxStdDev(i,j) = std(thisLegX);
            PEPxSEM(i,j) = std(thisLegX) / sqrt(length(thisLegX));

            thisLegY = allPEPYvals{i}(allPEPLegs{i} == legInd(j));
            PEPyMeans(i,j) = mean(thisLegY);
            PEPyStdDev(i,j) = std(thisLegY);
            PEPySEM(i,j) = std(thisLegY) / sqrt(length(thisLegY));
        end
    end

    % save into output file
    outFilePath = [saveDir filesep flyName '_condAEPPEP.mat'];

    save(outFilePath, 'allAEPXvals', 'allAEPYvals', 'allPEPXvals', ...
        'allPEPYvals', 'allAEPLegs', 'allPEPLegs', 'AEPxMeans', ...
        'AEPyMeans', 'AEPxStdDev', 'AEPyStdDev', 'AEPxSEM', 'AEPySEM', ...
        'PEPxMeans', 'PEPyMeans', 'PEPxStdDev', 'PEPyStdDev', 'PEPxSEM', ...
        'PEPySEM', 'stim', 'whichPhase', 'pDataFNames', ...
        'flipLegsLR', '-v7.3');
end