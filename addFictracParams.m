% addFictracParams.m
%
% Quick function to save FicTrac parameters from filtFictrac_all() into the
%  pData file. Corrects for bug (fixed 9/26/20) where it wasn't saved.
%
% INPUTS:
%   none, but prompts user to select pData files to add FicTrac parameters
%    to
%
% OUTPUTS:
%   none, but updates pData file with FicTrac parameters
%
% CREATED: 9/26/20 - HHY
%
% UPDATED: 9/26/20 - HHY
%
function addFictracParams()
    % FicTrac parameters to add
    fictracParams.dsf = 20; % downsample to 1000 Hz;
    fictracParams.filtParams.padLen = int32(200);
    fictracParams.filtParams.sigmaPos = int32(100); % 100 ms
    fictracParams.filtParams.sigmaVel = int32(50); % 50 ms
    
    BALL_DIAM = 6.46;
    circum = BALL_DIAM * pi; % circumference of ball, in mm
    fictracParams.mmPerDeg = circum / 360; % mm per degree of ball
    fictracParams.degPerMM = 360 / circum; % deg per mm ball
    
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
        if (contains(exptCond, 'Fictrac', 'IgnoreCase', true))
            
            % save fictracParams into pData file
            save(pDataFullPath, 'fictracParams', '-append');
            
            % display
            fprintf('Saved fictracParams for %s!\n', pDataName);
        else
            % display
            fprintf('%s does not have FicTrac data\n', pDataName);
        end
    end
end