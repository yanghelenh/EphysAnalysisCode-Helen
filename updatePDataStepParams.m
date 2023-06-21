% updatePDataStepParams.m
%
% Quick script to update legSteps struct after update to
%  computeStepParameters()
%
% CREATED: 2/13/23 - HHY
%

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
    
    % get variables in pData file
    pDataMatObj = matfile(pDataFullPath);
    pDataVarsStrct = whos(pDataMatObj);
    pDataVars = struct2cell(pDataVarsStrct);
    pDataVarsNames = pDataVars(1,:); % cell array of names
    
    % check that pData has legTrack, FicTrac, ephys, and Iinj data, 
    %  otherwise, skip
    if (any(contains(pDataVarsNames, 'legTrack')) && ...
            any(contains(pDataVarsNames,'fictracProc')) && ...
            any(contains(pDataVarsNames,'ephysSpikes')) && ...
            any(contains(pDataVarsNames,'legSteps')))

        % load pData
        load(pDataFullPath, 'legTrack','fictracProc', ...
            'ephysSpikes','legSteps');

        % compute step parameters
        legSteps = computeStepParameters(legSteps, legTrack, fictracProc,...
            ephysSpikes);
        
        % get swing/stance, currently, use duration method
        legSteps = callSwingStanceSteps(legSteps, legTrack, ...
            'duration', fictracProc);
        
        % step parameters during swing/stance
        [stanceStepParams, swingStepParams] = getStepParamsSwingStance(...
            legSteps);
        
        % update pData file
        save(pDataFullPath, 'legSteps', 'stanceStepParams', ...
            'swingStepParams', '-append');

    end
end