% logAnalysisStatus.m
%
% Function that reports on status of manual portions of analysis: selecting
%  dropped frames for FicTrac, getting filtered FicTrac signal, setting
%  threshold for ephys spike calls, loading in .trk data, leg step analysis
% Outputs as spreadsheet, each row is pData file.
% Select pData files to report on through GUI
%
% INPUTS:
%   sprdshtPath - full path for output spreadsheet
%   pDataPath - path to pData files, to start with
%
% OUTPUTS:
%   none, but outputs spreadsheet to sprdshtPath
%
% CREATED: 7/23/22 - HHY
%
% UPDATED:
%   7/23/22 - HHY
%   7/28/22 - HHY - updated to add opto
%
function logAnalysisStatus(sprdshtPath, pDataPath)
    
    % prompt user to select pData files
    [pDataFNames, pDataDirPath] = uigetfile('*.mat', ...
        'Select pData files', pDataPath, 'MultiSelect', 'on');
    
    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end

    % preallocate cell array that will become table written to spreadsheet
    tableArray = cell(numPDataFiles, 8);

    % loop through all pData files
    for i = 1:numPDataFiles
        % get full path of pData file(s)
        if (numPDataFiles == 1)
            thisPDataPath = [pDataDirPath filesep pDataFNames];
            thisPDataName = pDataFNames;
        else
            thisPDataPath = [pDataDirPath filesep pDataFNames{i}];
            thisPDataName = pDataFNames{i};
        end

        % all variables in pData file
        pDatVarNames = who('-file',thisPDataPath);

        % load exptCond, will inform which analyses need to be checked
        load(thisPDataPath,'exptCond');

        % FicTrac
        if contains(exptCond, 'fictrac','IgnoreCase',true)
            % need to load fictrac struct
            load(thisPDataPath,'fictrac');

            % check if fictrac struct as dropInd field
            if (isfield(fictrac,'dropInd')) % yes
                fictracDropInd = true;
            else
                fictracDropInd = false;
            end

            % check if filtFictrac_all has been run (fictracProc exists)

            % check if fictracProc is in pData variable names
            if (any(strcmp('fictracProc',pDatVarNames))) % yes
                fictracFilt = true;
            else
                fictracFilt = false;
            end
        else % NaN for indicator variables if this file has no fictrac 
            fictracDropInd = nan;
            fictracFilt = nan;
        end

        % ephysSpikes - generated when spikes called
        if contains(exptCond, 'ephys','IgnoreCase',true)
            % check if ephysSpikes is varaible in pData
            if (any(strcmp('ephysSpikes',pDatVarNames))) % yes
                ephysSpikeCalls = true;
            else
                ephysSpikeCalls = false;
            end
        else
            ephysSpikeCalls = nan;
        end

        % leg analysis

        if contains(exptCond, 'leg', 'IgnoreCase',true)
            % legTrack - generated when .trk file loaded in
            if (any(strcmp('legTrack',pDatVarNames))) % yes
                trkLoaded = true;
            else
                trkLoaded = false;
            end

            % moveNotMove - generated when moving/not moving called
            if (any(strcmp('moveNotMove',pDatVarNames)))
                notMoveCalled = true;
            else
                notMoveCalled = false;
            end

            % legSteps - generated when legSteps called
            if (any(strcmp('legSteps', pDatVarNames)))
                legStepsCalled = true;
            else
                legStepsCalled = false;
            end
        else
            trkLoaded = nan;
            notMoveCalled = nan;
            legStepsCalled = nan;
        end

        % opto

        if contains(exptCond, 'opto', 'IgnoreCase', true)
            if(any(strcmp('opto',pDatVarNames)))
                optoLoaded = true;
            else
                optoLoaded = false;
            end
        else
            optoLoaded = nan;
        end

        tableArray(i,:) = {thisPDataName, fictracDropInd, fictracFilt, ...
            ephysSpikeCalls, trkLoaded, notMoveCalled, legStepsCalled, ...
            optoLoaded};
    end

    % convert table array to table
    sprdshtTable = cell2table(tableArray);
    % add variable names to table
    sprdshtTable.Properties.VariableNames = {'pDataName',...
        'selectDroppedFictrac', 'filtFictrac_all', ...
        'interactEphysGetSpikes', 'loadTrk2PData','legMoveNotMove',...
        'legSteps', 'optoStim'};

    % write spreadsheet
    writetable(sprdshtTable, sprdshtPath, 'WriteVariableNames', true);
end