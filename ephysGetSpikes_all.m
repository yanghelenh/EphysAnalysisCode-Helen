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
    
    for i = 1:length(pDataFNames)
        
        pDataFullPath = [pDataPath pDataFNames{i}];
        
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
            fprintf('Saved ephysSpikes for %s!\n', pDataFNames{i});
            
        else
            % display
            fprintf('%s does not have ephys data\n', pDataFNames{i});
        end
    end

end