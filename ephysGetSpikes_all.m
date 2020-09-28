% ephysGetSpikes_all.m
%
% Function that detects spikes, computes spike rate, and computes
%  median-filtered voltage trace on each preprocessed ephys trace. Runs on
%  selected pData files. Must have run preprocess() before running this.
%
% INPUT:
%   none, but prompts user to select pData file(s)
%
% OUTPUTS:
%   none, but saves spike info into same pData file
%
% CREATED: 9/14/20 - HHY
%
% UPDATED:
%   9/14/20 - HHY
%   9/16/20 - HHY - correct bug with selecting single vs. multiple pData
%       files; single doesn't generate cell array
%
function ephysGetSpikes_all()
    % some constants
    % threshold in 1st derivative to detect spikes
    ephysSpikes.params.dvdtThresh = 8000; 
    ephysSpikes.params.refractoryPeriod = 2/1000; % 2 ms min time b/w spikes
    ephysSpikes.params.medFiltOrder = 0.05; % 50 ms median filter to remove spikes

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
        if (contains(exptCond, 'Ephys', 'IgnoreCase', true))
            load(pDataFullPath, 'ephysData');
            
            % detect spikes
            ephysSpikes.startInd = detectSpikes(...
                ephysData.scaledVoltage, ephysData.t, ...
                ephysSpikes.params.dvdtThresh, ...
                ephysSpikes.params.refractoryPeriod);
            
            % compute spike rate
            ephysSpikes.spikeRate = computeSpikeRate(...
                ephysSpikes.startInd, ephysData.t);
            
            % get median-filtered voltage trace, no spikes
            ephysSpikes.medFiltV = filtEphysNoSpikes(...
                ephysData.scaledVoltage, ephysData.t, ...
                ephysSpikes.params.medFiltOrder);
            
            % copy over time vector
            ephysSpikes.t = ephysData.t;
            
            % save back into pData
            save(pDataFullPath, 'ephysSpikes', '-append');
            
            % display
            fprintf('Saved ephysSpikes for %s!\n', pDataName);
            
        else
            % display
            fprintf('%s does not have ephys data\n', pDataName);
        end
    end

end