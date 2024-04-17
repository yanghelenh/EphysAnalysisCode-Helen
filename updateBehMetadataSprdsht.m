% preprocessOpto.m
%
% Function to update metadata spreadsheet for behavior-only experiments
%  run using the EphysCode-Helen codebase for each trial with all the
%  metadata that is saved during the experiment. Fields that require user
%  input posthoc are left blank.
% Assumes metadata spreadsheet has the arrangement of columns that matches
%  what this function writes.
% Adapted from updateMetadataSprdsht() in EphysAnalysisCode-Helen
%
% INPUTS:
%   sprdshtPath - full path to existing metadata spreadsheet
%   exptInfo - struct of experimental info (date, fly, cell) saved in 
%       metaDat.mat during runEphysExpt
%   flyData - struct of info about the fly, saved in metaDat.mat during
%       runEphysExpt
%   inputParams - struct of input parameters, specific to particular
%       experiment type, saved in trial.mat file when experiment run
%   trialName - string for name of trial (e.g. 'trial01' or 'cellAttached')
%
% OUTPUTS:
%   none, but writes to the metadata spreadsheet
%
% CREATED: 3/25/22 - HHY
%
% UPDATED:
%   3/28/22 - HHY
%   1/18/24 - HHY - update to account for different opto stims
%   
function updateBehMetadataSprdsht(sprdshtPath, exptInfo, flyData, ...
    inputParams, trialName)

    % get info about metadata spreadsheet
    opts = detectImportOptions(sprdshtPath);
    
    % read in first column of  metadata spreadsheet, to get number of rows
    opts.SelectedVariableNames = opts.VariableNames{1};
    sprdshtCol1 = readtable(sprdshtPath, opts);
    % number of rows in spreadsheet currently
    numSsRows = height(sprdshtCol1) + 1; % add 1 for heading row
    
    % read in fourth column of metadata spreadsheet, flyID
    opts.SelectedVariableNames = opts.VariableNames{4};
    sprdshtCol4 = readtable(sprdshtPath, opts);
    
    % figure out flyID for this row; either same as previous or own new one
    if(numSsRows == 1) % spreadsheet empty except for header
        thisFlyID = 1;
    else
        sprdshtCol4Arry = table2array(sprdshtCol4);
        % get flyID of last entry in spreadsheet
        lastFlyID = sprdshtCol4Arry(end);
        
        sprdshtCol1Arry = table2array(sprdshtCol1);
        % get experiment name for last entry in spreadsheet
        lastExptName = sprdshtCol1Arry{end};
        
        % this fly's date and fly number
        dateFlyNum = [flyData.dateDir '_' flyData.flyDir];
        
        % if the experiment name for the last entry contains the same date
        %  and fly number, it's the same fly (and has the same flyID)
        if(contains(lastExptName, dateFlyNum))
            thisFlyID = lastFlyID;
        else % if it doesn't contain it, then this is a new fly
            thisFlyID = lastFlyID + 1;
        end    
    end
    
    % generate experiment name
    exptName = [exptInfo.dateDir '_' exptInfo.flyDir '_' ...
        exptInfo.cellDir '_' trialName];
    
    % experiment type
    exptCond = inputParams.exptCond;

    % get info for opto experiments
    if (contains(inputParams.exptCond, 'opto','IgnoreCase',true))
        ndFilter = inputParams.optoStimParams.ndFilter;
        bpFilter = inputParams.optoStimParams.stimBPfilter;
        allStimDurs = num2str(inputParams.optoStimParams.allStimDurs);
        if isfield(inputParams.optoStimParams, 'durBwStims')
            durBwStims = inputParams.optoStimParams.durBwStims;
        elseif isfield(inputParams.optoStimParams, 'durBfStims')
            durBwStims = inputParams.optoStimParams.durBfStims + ...
                inputParams.optoStimParams.durAfStims;
        end
    else
        ndFilter = 'N/A';
        bpFilter = 'N/A';
        allStimDurs = 'N/A';
        durBwStims = 'N/A';
    end

    
    % starting cell to write new row to
    writeCell = sprintf('A%d', numSsRows + 1);
    
    % build cell array to write to excel spreadsheet
    rowArray = {exptName, exptCond, ...
        inputParams.startTimeStamp, thisFlyID, flyData.genotype, [], ...
        [], flyData.manipulation, flyData.prepType, ...
        flyData.age, flyData.ageUnits, flyData.dissectionNotes, ...
        ndFilter, bpFilter, allStimDurs, durBwStims, ...
        [], [], []};
    
    % convert cell array to table
    rowTable = cell2table(rowArray);
    
    
    % write to excel spreadsheet
%     xlswrite(sprdshtPath, rowArray, 1, writeCell); 
    %  xlswrite doesn't work on mac; use writetable instead
    writetable(rowTable, sprdshtPath, 'Range', writeCell, ...
        'WriteVariableNames', false);

end