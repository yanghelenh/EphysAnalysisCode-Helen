% preprocess.m
%
% First function to run on ephys/FicTrac/leg video/current injection data. 
%
% On ephys data:
%  Converts raw voltage traces collected on the DAQ into named 
%   voltage/current signals with the appropriate units. 
%  If option is selected, finds spikes, returns spike times and spike rate.
%   Also returns median filtered voltage trace, without spikes.
%  Associates pipette resistance, seal resistance, holding current, access 
%   resistance, input resistance, Vm rest with all trials for that cell. 
%  Cell attached trial gets its own pData.mat file.
% 
% On FicTrac data:
%  Converts raw voltage traces collected on experimental DAQ into position
%   and velocity traces for each of 3 axes (yaw/heading, forward, slide).
%   Also generates fictive path.
% 
% On leg tracking data:
%  Extracts frame time for each frame of leg tracking video. Saves this as
%   well as path, name, frame size info for legVid video, which is to be 
%   fed into APT. Copies leg video into user specified folder (for all leg
%   videos).
%
% All info saved into pData.mat file, which is saved into user specified
%  folder (for all pData.mat files).
%
% NOTE: assumes folder organization of date folder with one or more fly
% folders with cell folders.
%
% Adaptation of preprocess.m from 2PAnalysisCode-Helen
%
% CREATED: 9/3/20 - HHY
%
% UPDATED: 
%   9/10/20 - HHY
%   9/15/20 - HHY - fix bugs in handling cellAttachedTrial from prior to
%       change in preExperimentalRoutine (7/16/20)
%   10/4/21 - HHY - allow processing of FicTrac data through FicTrac vid
%
function preprocess()

    disp('Select metadata spreadsheet for experiment');
    [sprdshtName, sprdshtPath] = uigetfile('*.xlsx', ...
        'Select metadata spreadsheet');
    sprdshtFullPath = [sprdshtPath filesep sprdshtName];

    disp('Select a date folder to preprocess.');
    datePath = uigetdir;
    
    % get name of date folder
    fsLoc = strfind(datePath, filesep);
    dateFolder = datePath((fsLoc(end) + 1):end);
    
    curDir = pwd;
    cd(datePath);
    
    % get all fly folders in date directory    
    flyFolders = dir([datePath filesep 'fly*']);

    % loop through folders in date folder
    for i = 1:length(flyFolders)
        flyPath = [datePath filesep flyFolders(i).name];
        
        % get all cell folders in fly folder
        cellFolders = dir([flyPath filesep 'cell*']);
        
        % loop through cell folders in fly folder
        for j = 1:length(cellFolders)
            cellPath = [flyPath filesep cellFolders(j).name];
            
            % load metadata file for cell 
            metaDatFilePath = [cellPath filesep 'metaDat.mat'];
            load(metaDatFilePath, 'exptInfo', 'flyData', 'settings');
                
            % preprocess preExptTrials: load in ephys metadata about
            %  recording, preprocess cellAttached trial
            preExptPath = [cellPath filesep 'preExptTrials'];
            
            % make sure preExptTrials were actually collected before
            % analyzing
            if(isfolder(preExptPath))    
                % load pre-experimental data for cell
                preExptDataPath = [cellPath filesep 'preExptData.mat'];
                load(preExptDataPath, 'preExptData');
                
                cellAttMatPath = ...
                    [preExptPath filesep 'cellAttachedTrial.mat'];
                
                % if the file exists, preprocess cellAttached trial
                if(isfile(cellAttMatPath))
                    disp('Preprocessing cell attached trial');
                    
                    % to deal with data pre 7/16/20, when cellAttached
                    %  trial was changed to record behavior
                    % check variables in cellAttached.mat file (whether
                    %  processed or not) - check for rawData
                    cellAttVars = struct2cell(whos('-file', ...
                        cellAttMatPath));
                    cellAttVarsName = cellAttVars(1,:);
                    contRawData = sum(contains(cellAttVarsName, 'rawData'));
                    
                    % if cellAttached was not preprocessed (i.e. contains
                    %  rawData)
                    if (contRawData)
                        % load data
                        load(cellAttMatPath, 'inputParams', 'rawData', ...
                            'rawOutput');

                        % preprocess DAQ
                        [daqData, daqOutput, daqTime] = preprocessUserDaq(...
                            inputParams, rawData, rawOutput, settings);

                        % preprocess ephys
                        [ephysData, ephysMeta] = preprocessEphysData(...
                            daqData, daqOutput, daqTime, inputParams, settings);

                        % if there's behavioral data, preprocess that
                        % FicTrac
                        if(contains(inputParams.exptCond, 'Fictrac',...
                                'IgnoreCase', true))
                            fictrac = preprocessFicTrac(daqData, daqTime, ...
                                settings.bob.sampRate);
                        else % so writePData() has input
                            fictrac = [];
                        end

                        % leg video
                        if(contains(inputParams.exptCond, 'leg', ...
                                'IgnoreCase', true))
                            legVidPath = [preExptPath filesep ...
                                exptInfo.flyDir '_' exptInfo.cellDir ...
                                'cellAttachedTrial_legVid.mp4'];
                            leg = preprocessLegVid(legVidPath, daqData, ...
                                daqOutput, daqTime);
                        else % so writePData() has input
                            leg = [];
                        end

                        % update metadata spreadsheet
                        updateMetadataSprdsht(sprdshtFullPath, exptInfo, ...
                            flyData, inputParams, 'cellAttachedTrial', ...
                            preExptData);

                        % save pData
                        writePData(pDataDir(), settings, exptInfo, ...
                            preExptData, inputParams, ephysData, ephysMeta,...
                            fictrac, leg, 'cellAttachedTrial'); 

                        % clear variables specific to this trial
                        clearvars inputParams rawData rawOutput
                        clearvars ephysData ephysMeta fictrac leg
                    % otherwise, must have been ephysRecording, skip
                    %  preprocessing (only applies to data pre 7/16/20)
                    else
                        % load data
                        load(cellAttMatPath, 'ephysData', 'ephysMeta');
                        
                        % generate relevant fields of input params
                        inputParams.exptCond = 'ephysRecording'; 
                        inputParams.duration = ephysData.t(end);
                        inputParams.startTimeStamp = ...
                            ephysMeta.startTimeStamp;
                        inputParams.aInCh = {};
                        inputParams.aOutCh = {};
                        inputParams.dInCh = {};
                        inputParams.dOutCh = {};
                        
                        % empty vectors for fictrac and leg, since not
                        % present
                        fictrac = [];
                        leg = [];
                        
                        % update metadata spreadsheet
                        updateMetadataSprdsht(sprdshtFullPath, exptInfo, ...
                            flyData, inputParams, 'cellAttachedTrial', ...
                            preExptData);

                        % save pData
                        writePData(pDataDir(), settings, exptInfo, ...
                            preExptData, inputParams, ephysData, ephysMeta,...
                            fictrac, leg, 'cellAttachedTrial'); 
                        
                        % clear variables specific to this trial
                        clearvars inputParams ephysData ephysMeta
                        
                    end
                        
                end
            % if there's no preExpt data  (9/8/20 - not currently possible 
            %  if there are additional trials recorded properly, but handle
            %  in case recording done incorrectly or for future
            %  behavior-only experiments
            % updateMetadataSprdsht and writePData need preExptData input,
            %  but handle empty vector correctly
            else
                preExptData = [];
            end
            
            % get all trial files in cell folder
            trialFiles = dir([cellPath filesep 'trial*']);
            
            % loop through trials in cell folder
            for k = 1:length(trialFiles)
                trialPath = [cellPath filesep trialFiles(k).name];
                
                % name of trial, without .mat suffix
                trialName = trialFiles(k).name;
                periodInd = strfind(trialName,'.');
                trialName = trialName(1:(periodInd-1));
                
                % load trial 
                load(trialPath, 'inputParams', 'rawData', 'rawOutput');
                
                fprintf('\nPreprocessing %s\n', trialFiles(k).name);
                
                % preprocess DAQ
                [daqData, daqOutput, daqTime] = preprocessUserDaq(...
                    inputParams, rawData, rawOutput, settings);

                % if there's ephys data, preprocess that
                if(contains(inputParams.exptCond, 'ephys', 'IgnoreCase',...
                        true))
                    [ephysData, ephysMeta] = preprocessEphysData(...
                        daqData, daqOutput, daqTime, inputParams, ...
                        settings);
                else % so writePData() has input
                    ephysData = [];
                    ephysMeta = [];
                end

                % if there's behavioral data, preprocess that
                % FicTrac
                if(contains(inputParams.exptCond, 'Fictracvid', ...
                        'IgnoreCase', true))
                    
                    datFilePath = [cellPath filesep 'FicTrac' filesep ...
                        'Run3_newCalib' filesep trialName '_fictracVid.dat'];
                    fictrac = getFictracFromDat(datFilePath, daqData, ...
                        daqTime);
                    
                elseif(contains(inputParams.exptCond, 'Fictrac', ...
                        'IgnoreCase', true))
                    fictrac = preprocessFicTrac(daqData, daqTime, ...
                        settings.bob.sampRate);
                else % so writePData() has input
                    fictrac = [];
                end

                % leg video
                if(contains(inputParams.exptCond, 'leg', 'IgnoreCase',...
                        true))
                    legVidPath = [cellPath filesep exptInfo.flyDir ...
                        '_' exptInfo.cellDir '_' trialName ...
                        '_legVid.mp4'];
                    leg = preprocessLegVid(legVidPath, daqData, ...
                        daqOutput, daqTime);
                else % so writePData() has input
                    leg = [];
                end

                % update metadata spreadsheet
                updateMetadataSprdsht(sprdshtFullPath, exptInfo, ...
                    flyData, inputParams, trialName, ...
                    preExptData);

                % save pData
                writePData(pDataDir(), settings, exptInfo, ...
                    preExptData, inputParams, ephysData, ephysMeta,...
                    fictrac, leg, trialName);  
                
                % clear variables specific to this trial
                clearvars inputParams rawData rawOutput
                clearvars ephysData ephysMeta fictrac leg
                
            end
            
            % clear variabls specific to this cell
            clearvars exptInfo flyData settings preExptData
            
            % update display with what's happening; which cell done
            fprintf('Done preprocessing %s\n', cellFolders(j).name);
        end
        
        % update display with what's happening; which fly done
        fprintf('Done preprocessing %s\n', flyFolders(i).name);
    end
    
    % update display with what's happening; date folder done
    fprintf('Done preprocessing %s\n', dateFolder);
    
    cd(curDir);
end