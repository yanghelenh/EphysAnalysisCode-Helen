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
%
function filtFictrac_all()

    % some constants
    fictracParams.dsf = 20; % downsample to 1000 Hz;
    fictracParams.filtParams.padLen = int32(100);
    fictracParams.filtParams.sigmaPos = int32(50); % 100 ms
    fictracParams.filtParams.sigmaVel = int32(25); % 50 ms
    
    BALL_DIAM = 6.46;
    circum = BALL_DIAM * pi; % circumference of ball, in mm
    fictracParams.mmPerDeg = circum / 360; % mm per degree of ball
    fictracParams.degPerMM = 360 / circum; % deg per mm ball

    % prompt user to select pData files
    [pDataFNames, pDataPath] = uigetfile('*.mat', 'Select pData files', ...
        pDataDir(), 'MultiSelect', 'on');
    
    % loop through all pData files
    for i = 1:length(pDataFNames)
    
        pDataFullPath = [pDataPath pDataFNames{i}];
        
        % load pData
        load(pDataFullPath, 'exptCond');
        
        % check that pData has FicTrac data, otherwise, skip
        if (contains(exptCond, 'Fictrac', 'IgnoreCase', true))
            % load fictrac struct only if experiment has it
            load(fullPath, 'fictrac');
            
            % call dsFiltFictrac() to perform processing
            fictracProc = dsFiltFictrac(fictracParams, fictrac) ;
            
            % save fictracProc into pData file
            save(pDataFullPath, 'fictracProc', '-append');
            
            % display
            fprintf('Saved fictracProc for %s!\n', pDataFNames{i});
        else
            % display
            fprintf('%s does not have FicTrac data\n', pDataFNames{i});
        end
    end
    

end