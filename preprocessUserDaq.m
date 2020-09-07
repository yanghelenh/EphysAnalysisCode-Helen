% preprocessUserDaq.m
%
% Function to take raw data output from user DAQ for any experiment type
%  and convert it to a more useable form for later analyses
%
% Adapted from preprocessUserDaq.m from 2PAnalysisCode-Helen, but doesn't
%  incorporate information about fly or interact with a metadata
%  spreadsheet
%
% INPUTS:
%   inputParams - experiment-specific settings
%   rawData - data collected on DAQ during experiment
%   rawOutput - signal output on DAQ during experiment
%   settings - static settings about user DAQ
%
% OUTPUTS:
%   daqData - data collected on DAQ, with fields labeled
%   daqOutput - signal output on DAQ during experiment, with fields labeled
%   daqTime - timing vector for daqData and daqOutput
%   
%
% CREATED: 1/22/20
% UPDATED: 3/23/20 - HHY
%

function [daqData, daqOutput, daqTime] = preprocessUserDaq(...
    inputParams, rawData, rawOutput, settings)

    % initialize daqData, daqOutput, daqTime - prevents bug when there is
    % no input or output data
    daqData = [];
    daqOutput = [];
    daqTime = [];

    % convert rawData array into daqData struct, with named fields
    colNum = 1; % counter of columns in rawData array
    % analog input first
    for i = 1:length(inputParams.aInCh)
        % field names from inputParams
        daqData.(inputParams.aInCh{i}) = rawData(:, colNum);
        colNum = colNum + 1;
    end
    % digital input second
    for i = 1:length(inputParams.dInCh)
        daqData.(inputParams.dInCh{i}) = rawData(:, colNum);
        colNum = colNum + 1;
    end
    
    % convert rawOutput array into daqOutput struct, with name fields
    colNum = 1; % counter of columns in rawOutput array
    % analog input first
    for i = 1:length(inputParams.aOutCh)
        daqOutput.(inputParams.aOutCh{i}) = rawOutput(:, colNum);
        colNum = colNum + 1;
    end
    % digital output second
    for i = 1:length(inputParams.dOutCh)
        daqOutput.(inputParams.dOutCh{i}) = rawOutput(:, colNum);
        colNum = colNum + 1;
    end
    
    % extract timing vector for daqData, daqOutput
    numSamp = size(rawData,1);
    sampRate = settings.bob.sampRate;
    daqTime = (0:(numSamp-1))/sampRate;
       
end