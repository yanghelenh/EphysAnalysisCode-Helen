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
% UPDATED: 
%   9/3/20 - HHY
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
            
            % preprocess preExptTrials: load in ephys metadata about
            %  recording, preprocess cellAttached trial
            preExptPath = [cellPath filesep 'preExptTrials'];
            
            % make sure preExptTrials were actually collected before
            % analyzing
            if(isfolder(preExptPath)) 
                % load metadata file for cell 
                metaDatFilePath = [cellPath filesep 'metaDat.mat'];
                load(metaDatFilePath, 'exptInfo', 'flyData', 'settings');
                
                % load pre-experimental data for cell
                preExptDataPath = [cellPath filesep 'preExptData.mat'];
                load(preExptDataPath, 'preExptData');
                
                cellAttMatPath = ...
                    [preExptPath filesep 'cellAttachedTrial.mat'];
                
                % if the file exists, preprocess cellAttached trial
                if(isfile(cellAttMatPath))
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
                    if(contains(inputParams.exptCond, 'Fictrac'))
                        fictrac = preprocessFicTrac(daqData, daqTime, ...
                            settings.bob.sampRate);
                    end
                    % leg video
                    if(contains(inputParams.exptCond, 'leg'))
                        leg = preprocessLegVid(daqData, daqOutput, daqTime);
                    end
                    
                    % update metadata spreadsheet
                    
                end
                
                
                
                
            end
            
            
            % get all trial folders in cell folder
            trialFolders = dir([cellPath filesep 'trial*']);
            
            % loop through trial folders in FOV folder
            for k = 1:length(trialFolders)
                trialPath = [cellPath filesep trialFolders(k).name];
                exptName = [dateFolder '_' flyFolders(i).name '_' ...
                    cellFolders(j).name '_' trialFolders(k).name];
                cd (trialPath)
                
                % load data from experimental DAQ, metadata
                try
                    load([trialPath filesep 'userDaqDat.mat']);
                catch mExcep
                    % if there is no userDaqDat.mat file
                    if(strcmp(mExcep.identifier, ...
                            'MATLAB:load:couldNotReadFile'))
                        fprintf(['Error: no userDaqDat.mat file found '...
                            'for \n %s \n Skipping processing'], ... 
                            trialPath);
                        continue; % skips processing of this trial
                    % other errors, throw and stop preprocess completely    
                    else 
                        rethrow(ME); 
                    end
                end
                
                % display updates in command line
                fprintf('Preprocessing %s \n', trialPath);
                
                % process metadata from userDaqDat.mat
                [daqData, daqOutput, daqTime, settings] = ...
                    preprocessUserDaq(exptCond, flyData, inputParams, ...
                    rawData, rawOutput, settings, sprdshtFullPath,...
                    exptName);
                
                
                % if this experiment has imaging data
                if (contains(exptCond, 'Img'))
                    % name of ScanImage Tiff
                    tifFile = dir([trialPath filesep 'f*.tif']);

                    % only if imaging data exists
                    if ~isempty(tifFile)
                        disp('Preprocessing imaging data');
                        % align imaging data, save that and metadata in
                        %  trialPath
                        preprocessImaging(tifFile, daqData, daqTime);
%                         preprocessImaging2Ch(tifFile, daqData, daqTime);
                    else
                        fprintf(['Warning: Imaging data expected, but '
                            'no .tif file found for:\n %s \n'], trialPath);
                    end
                end
                
                % if this experiment has FicTrac data
                if (contains(exptCond, 'Fictrac'))
                    % Process FicTrac data
                    disp('Preprocessing FicTrac data');
                    preprocessFicTrac(daqData, daqTime, ...
                        settings.bob.sampRate);
                end
                
                % if this experiment has leg tracking data
                if (contains(exptCond, 'leg'))
                    % Process leg data
                    disp('Preprocessing leg video data');
                    preprocessLegVid(daqData, daqOutput, daqTime);  
                end
                
                % display updates in command line
                fprintf('Done preprocessing %s \n', trialPath);
                
            end
        end
    end
    
    cd(curDir);
end