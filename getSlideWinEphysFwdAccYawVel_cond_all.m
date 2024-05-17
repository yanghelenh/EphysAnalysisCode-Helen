% getSlideWinEphysFwdAccYawVel_cond_all.m
%
% Wrapper function for getSlideWinEphysFwdAccYawVel_cond_fly() that 
%  takes in pData files from multiple flies (either through GUI or input), 
%  parses them into separate flies, generates output file name, and feeds 
%  them into the single cell function (called single fly, but really single
%  cell)
%
% INPUTS:
%   same as getCorrelationEphysFwdAccYawVel_cond_fly(), but saveFileName is
%   suffix after cell name and pDataFNames is for all cells, not just one
%
% OUTPUTS:
%   none, but generates one saved file per cell
%
% CREATED: 5/15/24 - HHY
%
% UPDATED:
%   5/15/24 - HHY
%
function getSlideWinEphysFwdAccYawVel_cond_all(ephysParam, cond, ...
    normSpikeRate, fwdAccParams, winParams, tDelay, ...
    notMoveExclDur, postStimExclDur, pDataPath, pDataFNames, ...
    saveFileName, saveFileDir)

    
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
    fwdAccVals = []; % running tracker of all fwd acceleration time points
    yawVelVals = []; % running tracker of all yaw velocity time points
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
        [thisFwdAccVals, thisYawVelVals, thisEphysVals, thisEphysValsNorm] = ...
            getSlideWinEphysFwdAccYawVel_cond_fly(ephysParam, cond, ...
            normSpikeRate, fwdAccParams, winParams, tDelay, ...
            notMoveExclDur, postStimExclDur, pDataDirPath, ...
            thisCellPDataNames, thisSaveName, saveFileDir);

        % add to variable tracker across pData files
        fwdAccVals = cat(1,fwdAccVals, thisFwdAccVals);
        yawVelVals = cat(1,yawVelVals, thisYawVelVals);
        ephysVals = cat(1,ephysVals, thisEphysVals);
        ephysValsNorm = cat(1,ephysVals, thisEphysValsNorm);

        fprintf('Saved for %s!\n', uniCellNames{i});
    end
    
    % save combination across flies
    fullSavePath = [saveFileDir filesep saveFileName '_all.mat'];

    save(fullSavePath, 'fwdAccVals', 'yawVelVals', 'ephysVals', ...
        'ephysValsNorm', 'ephysParam', 'tDelay', 'fwdAccParams', ...
        'winParams', 'notMoveExclDur', 'postStimExclDur', 'cond', ...
        'normSpikeRate', '-v7.3');

end