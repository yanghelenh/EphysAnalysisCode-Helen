% filtFictrac_all.m
%
% Function that downsamples and filters FicTrac data of all user-selected
%  pData files. Runs by calling dsFiltFictrac(), which operates on a single
%  trial.
%
% Must have run selectDroppedFicTrac() on selected pData first.
%
% INPUTS:
%   none, but prompts user to select pData files to process
%
% OUTPUTS:
%   none, but saves processed FicTrac data back into pData file, as struct
%       fictracProc
%
% CREATED: 9/11/20 - HHY
%
% UPDATED: 
%   9/11/20 - HHY
%   9/16/20 - HHY - correct bug with selecting single vs. multiple pData
%       files; single doesn't generate cell array
%   9/26/20 - HHY - save fictracParams
%
function filtFictrac_all()

    % some constants
    fictracParams.dsf = 20; % downsample to 1000 Hz;
    fictracParams.filtParams.padLen = 200;
    fictracParams.filtParams.sigmaPos = 200; % 200 ms
    fictracParams.filtParams.sigmaVel = 100; % 100 ms
    
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
            % load fictrac struct only if experiment has it
            load(pDataFullPath, 'fictrac');
            
            % call dsFiltFictrac() to perform processing
            fictracProc = dsFiltFictrac(fictracParams, fictrac) ;
            
            % save fictracProc into pData file
            save(pDataFullPath, 'fictracProc', 'fictracParams', '-append');
            
            % display
            fprintf('Saved fictracProc for %s!\n', pDataName);
        else
            % display
            fprintf('%s does not have FicTrac data\n', pDataName);
        end
    end
    

end