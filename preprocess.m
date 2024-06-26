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
% On current injection data:
%  Run after preprocessing ephys data
%  Extracts start and end times, duration, and amplitude of each current
%   injection step. Saves parameters of current injection as well.
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
% On optogenetic stimulation data:
%  Extracts timing of optogenetic stimulation (start time, end time,
%   duration, logical for stim on/off). Also matches duration to commanded
%   duration to make averaging in later analyses easier. Saves this as well
%   as user-specified parameters about stimulation.
%
% On visual stimuli controlled in closed loop with voltage injection:
%  Extracts visual stimulus start and stop times and velocities. Saves
%   metadata about visual stimulus pattern and control
%  
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
%   3/25/22 - HHY - update handling of leg video files so that files with
%       and without the date in the file name are recognized (added date to
%       file name on 3/24/22)
%   3/28/22 - HHY - update to enable preprocessing of opto stim trials
%   6/28/22 - HHY - update to enable preprocessing of current injection
%       trials
%   1/16/24 - HHY - update to enable preprocessing of visual stimulus
%       trials
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

                        % preprocess current injection (unlike to exist
                        %  here, but include anyway)
                        if(contains(inputParams.exptCond, 'IInj',...
                                'IgnoreCase', true))
                            iInj = preprocessIInj(ephysData, inputParams);
                        else
                            iInj = [];
                        end

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
                            % changed leg video name to add date folder on
                            %  3/24/22
                            % accounts for name difference in leg video for
                            %  files made before and after that date
                            if (datetime(exptInfo.exptDate, ...
                                    'InputFormat', 'yyMMdd') < ...
                                datetime(2022,3,24))
    
                                legVidPath = [preExptPath filesep ...
                                    exptInfo.flyDir '_' exptInfo.cellDir ...
                                    'cellAttachedTrial_legVid.mp4'];
                            else
                                legVidPath = [preExptPath filesep ...
                                    dateFolder '_' exptInfo.flyDir '_' ...
                                    exptInfo.cellDir ...
                                    'cellAttachedTrial_legVid.mp4'];
                            end

                            leg = preprocessLegVid(legVidPath, daqData, ...
                                daqOutput, daqTime);
                        else % so writePData() has input
                            leg = [];
                        end

                        % opto (unlikely to actually exist here, but option
                        % anyway)
                        if(contains(inputParams.exptCond, 'opto',...
                                'IgnoreCase',true))
                            opto = preprocessOpto(daqData, daqTime, ...
                                inputParams);
                        else
                            opto = [];
                        end

                        % visstim - doesn't exist here
                        visstim = [];

                        % update metadata spreadsheet
                        updateMetadataSprdsht(sprdshtFullPath, exptInfo, ...
                            flyData, inputParams, 'cellAttachedTrial', ...
                            preExptData);

                        % save pData
                        writePData(pDataDir(), settings, exptInfo, ...
                            preExptData, inputParams, ephysData, ephysMeta,...
                            fictrac, leg, opto, iInj, visstim, ...
                            'cellAttachedTrial'); 

                        % clear variables specific to this trial
                        clearvars inputParams rawData rawOutput
                        clearvars ephysData ephysMeta fictrac leg opto visstim
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
                        
                        % empty vectors for not present variables
                        fictrac = [];
                        leg = [];
                        opto = [];
                        iInj = [];
                        visstim = [];
                        
                        % update metadata spreadsheet
                        updateMetadataSprdsht(sprdshtFullPath, exptInfo, ...
                            flyData, inputParams, 'cellAttachedTrial', ...
                            preExptData);

                        % save pData
                        writePData(pDataDir(), settings, exptInfo, ...
                            preExptData, inputParams, ephysData, ephysMeta,...
                            fictrac, leg, opto, iInj, visstim,...
                            'cellAttachedTrial'); 
                        
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

                % if there's current injection data, preprocess that
                if (contains(inputParams.exptCond, 'IInj', 'IgnoreCase',...
                        true))
                    iInj = preprocessIInj(ephysData, inputParams);
                else 
                    iInj = [];
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
                    % added date to leg video name on 3/24/22
                    % accounts for difference in file name before and after
                    %  that date
                    if (datetime(exptInfo.exptDate, ...
                        'InputFormat', 'yyMMdd') < datetime(2022,3,24))

                        legVidPath = [cellPath filesep exptInfo.flyDir ...
                            '_' exptInfo.cellDir '_' trialName ...
                            '_legVid.mp4'];
                    else
                        legVidPath = [cellPath filesep dateFolder '_' ...
                            exptInfo.flyDir '_' exptInfo.cellDir '_' ...
                            trialName '_legVid.mp4'];
                    end

                    leg = preprocessLegVid(legVidPath, daqData, ...
                        daqOutput, daqTime);
                else % so writePData() has input
                    leg = [];
                end

                % opto stimulation
                if(contains(inputParams.exptCond, 'opto', 'IgnoreCase',...
                        true))
                    opto = preprocessOpto(daqData, daqTime, inputParams);
                else % so writePData() has input
                    opto = [];
                end

                % visual stimulus controlled in closed loop with voltage
                %  injection from DAQ
                if (contains(inputParams.exptCond, 'visstim', ...
                        'IgnoreCase',true) && contains(...
                        inputParams.exptCond,'VInj','IgnoreCase',true))
                    visstim = preprocessVisstimVInj(daqData, daqTime, ...
                        inputParams, daqOutput);
                else % so writePData() has input
                    visstim = [];
                end


                % update metadata spreadsheet - different types for
                % behavior only vs. ephys
                % for ephys experiments
                if(contains(inputParams.exptCond, 'ephys','IgnoreCase',...
                        true))
                    updateMetadataSprdsht(sprdshtFullPath, exptInfo, ...
                        flyData, inputParams, trialName, ...
                        preExptData);
                % for behavior only experiments    
                else
                    updateBehMetadataSprdsht(sprdshtFullPath, exptInfo,...
                        flyData, inputParams, trialName);
                end

                % save pData
                writePData(pDataDir(), settings, exptInfo, ...
                    preExptData, inputParams, ephysData, ephysMeta,...
                    fictrac, leg, opto, iInj, visstim, trialName);  
                
                % clear variables specific to this trial
                clearvars inputParams rawData rawOutput
                clearvars ephysData ephysMeta fictrac leg opto
                
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