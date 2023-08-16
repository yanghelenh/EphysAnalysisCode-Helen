% getCorrelationEphysContParam_cond_all.m
%
% Wrapper function for getCorrelationEphysContParam_cond_fly() that takes
%  in pData files from multiple flies (either through GUI or input), parses
%  them into separate flies, generates output file name, and feeds them 
%  into the single cell function (called single fly, but really single
%  cell)
%
% INPUTS:
%   same as getCorrelationEphysContParam_cond_fly(), but saveFileName is
%   suffix after cell name and pDataFNames is for all cells, not just one
%
% OUTPUTS:
%   none, but generates one saved file per cell
%
% CREATED: 8/11/23 - HHY
%
% UPDATED:
%   8/11/23 - HHY
%
function getCorrelationEphysContParam_cond_all(ephysParam, behParams, ...
    legs, cond, normSpikeRate, tDelay, notMoveExclDur, postStimExclDur, ...
    pDataPath, pDataFNames, saveFileName, saveFileDir)

    
    % prompt user to select pData files if none given
    if isempty(pDataFNames)
        [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
            'Select pData files', pDataPath, 'MultiSelect', 'on');
    else
        pDataDirPath = pDataPath;
    end

    % first 19 characters defines name of fly
    for i = 1:length(pDataFNames)
        allCellNames{i} = pDataFNames{i}(1:12);
    end

    % get unique cell names
    uniCellNames = unique(allCellNames);

    % initialize
    ephysVals = []; % running tracker of all ephys time points
    behVals = []; % running tracker of all behavior time points
    ephysValsNorm = [];

    % loop through all unique cells
    % generate name to save
    for i = 1:length(uniCellNames)
        % get which pDataFiles belong to this cell
        thisCellPDataNames = pDataFNames(contains(pDataFNames, ...
            uniCellNames{i}));

        % generate name to save this file
        thisSaveName = [uniCellNames{i} '_' saveFileName];

        % get output for this cell
        [thisBehVals, thisEphysVals, thisEphysValsNorm] = ...
            getCorrelationEphysContParam_cond_fly(ephysParam, behParams, ...
            legs, cond, normSpikeRate, tDelay, notMoveExclDur, ...
            postStimExclDur, pDataDirPath, thisCellPDataNames, ...
            thisSaveName, saveFileDir);

        % add to variable tracker across pData files
        behVals = cat(1,behVals, thisBehVals);
        ephysVals = cat(1,ephysVals, thisEphysVals);
        ephysValsNorm = cat(1,ephysVals, thisEphysValsNorm);

        fprintf('Saved for %s!\n', uniCellNames{i});
    end

    % if multiple behavioral variables, get projection to equal weighting
    if iscell(behParams)
        % get coefficients
        coeffs = ones(length(behParams), 1) * (1/length(behParams));
        behVals1D = getLinProj(behVals, coeffs);
    % otherwise, just original behavioral value    
    else
        behVals1D = behVals;
    end
    
    % save combination across flies
    fullSavePath = [saveFileDir filesep saveFileName '_all.mat'];

    save(fullSavePath, 'behVals', 'ephysVals', 'ephysValsNorm', ...
        'behVals1D', 'ephysParam', 'behParams', 'legs', 'tDelay', ...
        'notMoveExclDur', 'postStimExclDur', 'cond', 'normSpikeRate', ...
        '-v7.3');

end