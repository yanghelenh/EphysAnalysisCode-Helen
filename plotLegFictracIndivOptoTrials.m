% plotLegFictracIndivOptoTrials.m
%
% Function to plot leg and FicTrac data for individual opto stimulation
%  trials for a single fly
% User specifies which leg and FicTrac parameters, by name of fields in
%  legTrackOpto and fictracOpto structs
% Data expected: output from computeLegFictracOpto_1Fly(). Load in through
%  GUI
%
% INPUTS:
%   paramNames - cell array of names of leg and FicTrac parameters to plot,
%       must match names of fields of legTrackOpto and fictracOpto structs
%   dataDir - directory in which to look for data
%   whichTrials - vector of indices of which trials to plot; invalid
%       indices are ignored
%   whichND - scalar value, which ND to plot
%   whichDur - scalar value, which duration of stimulation to plot
%   legInd - vector of indices of leg points, in order [R1 R2 R3 L1 L2 L3]
%   yScale - y axis scale for all parameters, as n parameters by 2 matrix,
%       where column 1 is min and 2 is max
%
% OUTPUTS:
%   none, but produces plot
%
% CREATED: 5/16/22 - HHY
%
% UPDATED:
%   5/16/22 - HHY
%
function plotLegFictracIndivOptoTrials(paramNames, dataDir, whichTrials,...
    whichND, whichDur, legInd, yScale)

    % prompt user to select avgLegFictracOpto file
    [dataFileName, dataFilePath] = uigetfile('*_avgLegFictracOpto.mat', ...
        'Select avgLegFictracOpto file', dataDir, 'MultiSelect', 'off');

    % load avgLegFictracOpto file
    load([dataFilePath filesep dataFileName], 'legTrackOpto', ...
        'fictracOpto', 'durs', 'NDs');

    % convert whichND and whichDur to indices into legTrackOpto and
    %  fictracOpto cell arrays
    ndInd = find(whichND == NDs);
    durInd = find(whichDur == durs);

    % loop through all parameters to plot
end