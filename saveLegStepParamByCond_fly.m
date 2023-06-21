% saveLegStepParamByCond_fly.m
%
% Function that takes in user selected pData files (for 1 fly) on which
%  sortLegStepsByIInj or sortLegStepsByOpto has been run and saves mean, 
%  and std dev for specified legStep parameter for each condition for that 
%  fly
% Can condition on a different value of legStep params
%
% INPUTS:
%   legStepParam - string for name of legStep parameter to extract data
%   stim - struct of parameters about which stimulation conditions to plot,
%     if any, with fields:
%       whichStim - which type of stimulation to condition on, specified as
%           string 'iInj' or 'opto'. [] for no conditioning
%       durs - stim durations to include, for iInj or opto
%       NDs or amps - specify stimulus amplitude, NDs for opto, amps for
%           iInj
%   cond - struct of parameters about selecting specific legStep parameter
%     to condition on. [] for no conditioning
%       whichParam - string of legStep parameter to condition on
%       cond - cell array of strings to condition on, for eval(), of size
%           numDurs * numNDs/amps + 1; matches order in stim struct (all  
%           durs for each ND/amp)
%   whichPhase - which phase of leg movement, specified as string 'swing'
%     or 'stance'
%   flipLegsLR - logical for whether to flip leg assignments over midline
%   saveDir - folder in which to save file, will be named date_fly_cell
%       after 1st pData file
%
% OUTPUTS:
%   none, but generates saved file
%
% CREATED: 2/14/23 - HHY
%
% UPDATED:
%   2/14/23 - HHY
%   6/14/23 - HHY - fix bug for computing mean and std dev on circular stat
%       (stepDirections)
%
function saveLegStepParamByCond_fly(legStepParam, stim, cond, ...
    whichPhase, flipLegsLR, saveDir)

    % all step parameters that are circular variables - need to use
    %  circular stats - 6/14/23 - HHY
    circStepParams = {'stepDirections'};

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

    allStepVals = cell(numCats,1); % all legStep values
    allStepLegs = cell(numCats,1); % which leg the steps belong to

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

                % get legStep param value, for given phase
                stepValAll = groupStepParamsBySwingStance(...
                    legSteps.(legStepParam),legSteps.stepSwingStance,...
                    phaseInd);

                % cond legStep vector, if present
                if ~isempty(cond)
                    stepCond = groupStepParamsBySwingStance(...
                        legSteps.(cond.whichParam),...
                        legSteps.stepSwingStance, phaseInd);
                end

                % get step categories for this phase
                thisLegStepCat = groupStepParamsBySwingStance(...
                    legStepsByIinj.stepIinjCat, legSteps.stepSwingStance,...
                    phaseInd);

                % loop through all conditions
                for j = 1:length(iInjCatAmps)
                    % get just for this amp and dur combo
                    % with conditioning
                    if ~isempty(cond)
                        [stepVal, theseWhichLeg] = getCatStepValsCond(...
                            stepValAll, iInjCatAmps(j), iInjCatDurs(j), ...
                            stepCond, cond.cond{j}, legSteps.stepWhichLeg, ...
                            thisLegStepCat, legStepsByIinj.iInjCatAmps, ...
                            legStepsByIinj.iInjCatDurs);
                    else % without conditioning
                        [stepVal, theseWhichLeg] = getCatStepVals(...
                            stepValAll, iInjCatAmps(j),iInjCatDurs(j),...
                            legSteps.stepWhichLeg, thisLegStepCat, ...
                            legStepsByIinj.iInjCatAmps, ...
                            legStepsByIinj.iInjCatDurs);
                    end
                    % append to values from across pData files
                    allStepVals{j} = [allStepVals{j}; stepVal];
                    allStepLegs{j} = [allStepLegs{j}; theseWhichLeg];
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

                % get legStep param value, for given phase
                stepValAll = groupStepParamsBySwingStance(...
                    legSteps.(legStepParam),legSteps.stepSwingStance,...
                    phaseInd);

                % cond legStep vector, if present
                if ~isempty(cond)
                    stepCond = groupStepParamsBySwingStance(...
                        legSteps.(cond.whichParam),...
                        legSteps.stepSwingStance, phaseInd);
                end

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
                    % with just this opto and dur combo
                    % with conditioning
                    if ~isempty(cond)
                        [stepVal, theseWhichLeg] = getCatStepValsCond(...
                            stepValAll, optoCatNDs(j), optoCatDurs(j),...
                            stepCond, cond.cond{j}, legSteps.stepWhichLeg,...
                            thisLegStepCat, ndKey, legStepsByOpto.optoCat);
                    else % no conditioning

                        [stepVal, theseWhichLeg] = getCatStepVals(...
                            stepValAll, optoCatNDs(j), optoCatDurs(j),...
                            legSteps.stepWhichLeg, thisLegStepCat, ...
                            ndKey, legStepsByOpto.optoCat);
                    end
                    % append to values from across pData files
                    allStepVals{j} = [allStepVals{j}; stepVal];
                    allStepLegs{j} = [allStepLegs{j}; theseWhichLeg];
                end
            end
        % invalid stim specification
        else
            disp('Invalid input for whichStim. Ending');
            return;
        end
    end

    % preallocate, matrices for means, stddev
    stepValMeans = nan(length(allStepVals), NUM_LEGS);
    stepValStdDev = nan(length(allStepVals), NUM_LEGS);
    stepValSEM = nan(length(allStepVals), NUM_LEGS);

    % get means and std dev, by leg
    for i = 1:length(allStepVals)
        for j = 1:NUM_LEGS
            thisLeg = allStepVals{i}(allStepLegs{i} == legInd(j));

            if (~isempty(allStepVals{i}))
            
                % if this needs circular stats
                if(any(strcmpi(legStepParam, circStepParams)))
                    % convert this parameter to radians
                    thisLegRad = deg2rad(thisLeg);
                    % get circular mean; convert back to degrees
                    thisMean = rad2deg(circ_mean(thisLegRad));
                    % get circular std; convert back to degrees
                    thisStd = rad2deg(circ_std(thisLegRad));
    
                    stepValMeans(i,j) = thisMean;
                    stepValStdDev(i,j) = thisStd;
                    stepValSEM(i,j) = thisStd / sqrt(length(thisLeg));
    
                % not circular stats
                else    
                    stepValMeans(i,j) = mean(thisLeg);
                    stepValStdDev(i,j) = std(thisLeg);
                    stepValSEM(i,j) = std(thisLeg) / sqrt(length(thisLeg));
                end
            else
                stepValMeans(i,j) = nan;
                stepValStdDev(i,j) = nan;
                stepValSEM(i,j) = nan;
            end
        end
    end

    % save into output file
    outFilePath = [saveDir filesep flyName '_cond_' legStepParam '.mat'];

    save(outFilePath, 'allStepVals', 'allStepLegs', 'stepValMeans', ...
        'stepValStdDev', 'stepValSEM', 'legStepParam', 'stim', ...
        'whichPhase', 'flipLegsLR', 'pDataFNames', 'cond', '-v7.3');
end