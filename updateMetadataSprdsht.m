% updateMetadataSprdsht.m
%
% Function to update ephys metadata spreadsheet for each trial with all the
%  metadata that is saved during the experiment. Fields that require user
%  input posthoc are left blank.
% Assumes metadata spreadsheet has the arrangement of columns that matches
%  what this function writes.
% Adapted from the last section of preprocessUserDaq.m from
%  2PAnalysisCode-Helen
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
%   preExptData - struct of pre-experimental detail about ephys recording
%       (pipette resistance, etc.), saved in preExptData.mat during
%       runEphysExpt call to preExptRoutines; feed in [] if
%       experiment type doesn't generate preExptData file
%
% OUTPUTS:
%   none, but writes to the metadata spreadsheet
%
% CREATED: 9/6/20 - HHY
%
% UPDATED:
%   9/6/20 - HHY
%   9/8/20 - HHY - update to correct issue with exptCond for
%       legFictracEphysIInj prior to 9/8/20
%   
function updateMetadataSprdsht(sprdshtPath, exptInfo, flyData, ...
    inputParams, trialName, preExptData)

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
        % get flyID of last entry in spreadsheet
        lastFlyID = sprdshtCol4(end);
        
        % get experiment name for last entry in spreadsheet
        lastExptName = sprdshtCol1(end);
        
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
    
    % propagate which cell, pipette internal info; accounts for 9/6/20 
    %  update to runEphysExpt, preExptRoutine that asks for this data
    if (datetime(exptInfo.exptDate, 'InputFormat', 'yyMMdd') < ...
            datetime(2020,9,6))
        whichCell = []; % user has to enter manually
        intSln = [];
    else
        whichCell = exptInfo.cellInfo;
        % check if preExptData struct exists; if experiment type doesn't
        %  have, return N/A
        if(~isempty(preExptData))
            intSln = preExptData.internal;
        else
            intSln = 'N/A';
        end
    end
    
    % pre-experimental routine data; check if it exists, check if each
    %  field exists
    if(~isempty(preExptData))
        % pipette resistance
        if(isfield(preExptData, 'pipetteResistance'))
            pipetteResistance = preExptData.pipetteResistance;
        else
            pipetteResistance = 'N/A';
        end
        
        % seal resistance
        if(isfield(preExptData, 'sealResistance'))
            sealResistance = preExptData.sealResistance;
        else
            sealResistance = 'N/A';
        end
        
        % initial holding current
        if(isfield(preExptData, 'initialHoldingCurrent'))
            holdCurrent = preExptData.initialHoldingCurrent;
        else
            holdCurrent = 'N/A';
        end
        
        % access resistance
        if(isfield(preExptData, 'initialAccessResistance'))
            accessResistance = preExptData.initialAccessResistance;
        else
            accessResistance = 'N/A';
        end
        
        % input resistance
        if(isfield(preExptData, 'initialInputResistance'))
            inputResistance = preExptData.initialInputResistance;
        else
            inputResistance = 'N/A';
        end
        
        % Vm rest
        if(isfield(preExptData, 'initialRestingVoltage'))
            restVm = preExptData.initialRestingVoltage;
        else
            restVm = 'N/A';
        end
    % if experiment type doesn't have pre-experimental data, all fields are
    %  'N/A'
    else 
        pipetteResistance = 'N/A';
        sealResistance = 'N/A';
        holdCurrent = 'N/A';
        accessResistance = 'N/A';
        inputResistance = 'N/A';
        restVm = 'N/A';
    end
    
    % correct exptCond for legFictracEphysIInj prior to 9/8/20
    if ((datetime(exptInfo.exptDate, 'InputFormat', 'yyMMdd') < ...
        datetime(2020,9,8)) && ...
        (strcmpi(inputParams.exptCond,'legFictracEphys')) && ...
        (isfield(inputParams, 'iInjProtocol')))
        exptCond = 'legFictracEphysIInj';
    else
        exptCond = inputParams.exptCond;
    end
    
    % starting cell to write new row to
    writeCell = sprintf('A%d', numSsRows + 1);
    
    % build cell array to write to excel spreadsheet
    rowArray = {exptName, exptCond, ...
        inputParams.startTimeStamp, thisFlyID, flyData.genotype, [], ...
        whichCell, flyData.manipulation, flyData.prepType, ...
        flyData.age, flyData.ageUnits, flyData.dissectionNotes, ...
        intSln, pipetteResistance, sealResistance, holdCurrent, ...
        accessResistance, inputResistance, restVm, [], [], []};
    
    % convert cell array to table
    rowTable = cell2table(rowArray);
    
    
    % write to excel spreadsheet
%     xlswrite(sprdshtPath, rowArray, 1, writeCell); 
    %  xlswrite doesn't work on mac; use writetable instead
    writetable(rowTable, sprdshtPath, 'Range', writeCell, ...
        'WriteVariableNames', false);

end