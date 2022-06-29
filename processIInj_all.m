% processIInj_all.m
%
% Function that processes current injection data from ephyData in
%  user-selected pData file(s)
% Calls preprocessIInj(). Meant to be run on pData processed before
%  preprocessIInj() was added to preprocess()
%
% INPUTS:
%   none, but prompts user to select pData files to process
%
% OUTPUTS:
%   none, but saves processed iInj data back into pData file, as struct
%       iInj
%
% CREATED: 6/28/22 - HHY
%
% UPDATED:
%   6/28/22 - HHY
%
function processIInj_all()

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
    
    % loop through all pData files
    for i = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{i};
        else
            pDataName = pDataFNames;
        end
        
        pDataFullPath = [pDataPath pDataName];
        
        % load pData
        load(pDataFullPath, 'exptCond');
        
        % check that pData has FicTrac data, otherwise, skip
        if (contains(exptCond, 'IInj', 'IgnoreCase', true))
            % load ephysData struct only if experiment has iInj data
            load(pDataFullPath, 'ephysData', 'inputParams');

            % call preprocessIInj to do processing
            iInj = preprocessIInj(ephysData, inputParams);

            % save iInj into pData file
            save(pDataFullPath, 'iInj', '-append');
            
            % display
            fprintf('Saved iInj for %s!\n', pDataName);
        else
            % display
            fprintf('%s does not have iInj data\n', pDataName);
        end
    end
end