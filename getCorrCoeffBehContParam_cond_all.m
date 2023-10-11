% getCorrCoeffBehContParam_cond_all.m
%
% Function that gets all pairwise correlation coefficients for all input
%  continuous behavioral parameters, across all flies (specified as list of
%  pData file names).
% Calls getCorrCoeffBehContParam_cond_fly()
%
% INPUTS:
%   allBehParams - array of structs for behavioral params to get
%     correlations for; each struct has following fields:
%       whichParams - cell array of string(s) for which behavioral
%         parameters to correlate. Must be members of FicTrac or continuous
%         estimate of step parameters
%       legs - single string or cell array of strings for which legs, when step
%         parameters. [] when FicTrac param
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
%     on. For all flies
%   saveFileName - suffix to include in name of output file, after cell
%     name
%   saveFileDir - full path to directory to save output file
%
% CREATED: 10/7/23 - HHY
%
% UPDATED:
%   10/7/23 - HHY
%
function getCorrCoeffBehContParam_cond_all(allBehParams, cond, ...
    notMoveExclDur, postStimExclDur, pDataPath, pDataFNames, ...
    saveFileName, saveFileDir)


    % prompt user to select pData files if none given
    if isempty(pDataFNames)
        [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
            'Select pData files', pDataPath, 'MultiSelect', 'on');
    else
        pDataDirPath = pDataPath;
    end

    % first 12 characters defines name of fly
    for i = 1:length(pDataFNames)
        allCellNames{i} = pDataFNames{i}(1:12);
    end

    % get unique cell names
    uniCellNames = unique(allCellNames);

    % loop through all unique cells
    % generate name to save
    for i = 1:length(uniCellNames)
        % get which pDataFiles belong to this cell
        thisCellPDataNames = pDataFNames(contains(pDataFNames, ...
            uniCellNames{i}));

        % preallocate matrix for all pairwise correlations
        allCorr = zeros(length(allBehParams));

        % loop through all behParams as 1st behParam
        for j = 1:length(allBehParams)
            % loop through all behParams as 2nd behParam
            for k = 1:length(allBehParams)
                behParams1 = allBehParams(j);
                behParams2 = allBehParams(k);

                allCorr(j,k) = getCorrCoeffBehContParam_cond_fly(...
                    behParams1, behParams2, cond, notMoveExclDur, ...
                    postStimExclDur, pDataDirPath, thisCellPDataNames);
            end
        end


        % save all pairwise correlations for this fly
        fullSavePath = ...
            [saveFileDir filesep uniCellNames{i} '_' saveFileName '.mat'];

        save(fullSavePath, 'allCorr', 'allBehParams', 'cond',...
        'notMoveExclDur', 'postStimExclDur', '-v7.3');

        fprintf('Saved for %s!\n', uniCellNames{i});
    end
end
