% plotLegStepsHist.m
%
% Function that plots histograms of leg step parameters. Multiple
%  histograms per figure, by specified conditions (iInj or opto).
% Works on multiple pData files. Select through GUI. Should be for 1 fly.
%
% INPUTS:
%   whichParam - which parameter to plot, as string; must match field of
%     legSteps
%   stim - struct of parameters about which stimulation conditions to plot,
%     if any, with fields:
%       whichStim - which type of stimulation to condition on, specified as
%           string 'iInj' or 'opto'. [] for no conditioning
%       durs - stim durations to include, for iInj or opto
%       NDs or amps - specify stimulus amplitude, NDs for opto, amps for
%           iInj
%   whichPhase - which phase of leg movement, specified as string 'swing'
%       or 'stance'
%   
% OUTPUTS:
%   none - but produces plot
%
% CREATED: 10/4/22 - HHY
%
% UPDATED:
%   10/4/22 - HHY
%
function plotLegStepsHist(whichParam, stim, whichPhase, binWidth)

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
        otherwise
            numCats = 1;
    end
    allStepVals = cell(numCats,1); % values for X
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

                % leg step parameter as specified by input
                % returns single vector
                stepVal = groupStepParamsBySwingStance(...
                    legSteps.(whichParam),...
                    legSteps.stepSwingStance, phaseInd);

                % get step categories for this phase
                thisLegStepCat = groupStepParamsBySwingStance(...
                    legStepsByIinj.stepIinjCat, legSteps.stepSwingStance,...
                    phaseInd);

                % loop through all conditions
                for j = 1:length(iInjCatAmps)
                    % get values for just this condition
                    [theseStepVal, theseWhichLeg] = getCatStepVals(...
                        stepVal, iInjCatAmps(j),...
                        iInjCatDurs(j), legSteps.stepWhichLeg, ...
                        thisLegStepCat, legStepsByIinj.iInjCatAmps, ...
                        legStepsByIinj.iInjCatDurs);
                    % append to values from across pData files
                    allStepVals{j} = [allStepVals{j}; theseStepVal];
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

                % leg step parameter as specified by input
                % returns single vector
                stepVal = groupStepParamsBySwingStance(...
                    legSteps.(whichParam),...
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
                    % get values for just this condition
                    [theseStepVal, theseWhichLeg] = ...
                        getCatStepVals(stepVal, optoCatNDs(j),...
                        optoCatDurs(j), legSteps.stepWhichLeg, ...
                        thisLegStepCat, ndKey, legStepsByOpto.optoCat);
                    % append to values from across pData files
                    allStepVals{j} = [allStepVals{j}; theseStepVal];
                    allStepLegs{j} = [allStepLegs{j}; theseWhichLeg];
                end
            end
        % no conditioning
        elseif (isempty(stim.whichStim))
            % check that pData has appropriate structs; otherwise, skip
            if (any(contains(pDataVarsNames, 'legSteps')))

                % load pData
                load(pDataFullPath, 'legSteps');

                % leg step parameter as specified by input
                % returns single vector
                stepVal = groupStepParamsBySwingStance(...
                    legSteps.(whichParam),...
                    legSteps.stepSwingStance, phaseInd);

                % append all values, cell array size 1
                allStepVals{1} = [allStepVals{1}; stepVal];
                allStepLegs{1} = [allStepLegs{1}; legSteps.stepWhichLeg];

            end
        % invalid stim specification
        else
            disp('Invalid input for whichStim. Ending');
            return;
        end
    end

    % PLOTTING

    % initialize cell array for legend
    legendStr = cell(1,length(allStepVals));

    figure;

%     % use colormap lines
%     c = colormap('lines');
    for l = 1:6
        for i = 1:length(allStepVals)
            thisLegStepVals = allStepVals{i}(allStepLegs{i}==l);
            
            subplot(2,3,l);
            histogram(thisLegStepVals, 'Normalization', 'probability',...
                'BinWidth', binWidth);
    
            hold on;
    
            if (i==1)
                legendStr{i} = 'No stim';
            else
                if (strcmpi(stim.whichStim, 'iInj'))
                    legendStr{i} = sprintf('Amp = %.1f pA, Dur = %.1f s', ...
                        iInjCatAmps(i), iInjCatDurs(i));
                elseif (strcmpi(stim.whichStim, 'opto'))
                    legendStr{i} = sprintf('ND = %.1f, Dur = %.1f s', ...
                        optoCatNDs(i), optoCatDurs(i));
                end
            end
        end

    end
    % title of plot
    ttlStr = sprintf('%s during %s', whichParam, whichPhase);
    sgtitle(ttlStr);
    
    legend(legendStr); % legend 
end

% helper function for returning steps for 1 pData file for 1 category
function [theseStepVal, theseWhichLeg] = getCatStepVals(stepVal, ...
    whichAmp, whichDur, stepWhichLeg, legStepCat, ampsKey, dursKey)

    % turn whichAmp and whichDur into index value for legStepCat
    thisInd = intersect(find(ampsKey == whichAmp), ...
        find(dursKey == whichDur));

    % get indices into steps that match this categories
    valInd = find(legStepCat == thisInd);

    % use indices to select values of extreme point, which leg that are
    %  valid
    theseStepVal = stepVal(valInd);
    theseWhichLeg = stepWhichLeg(valInd);   
end