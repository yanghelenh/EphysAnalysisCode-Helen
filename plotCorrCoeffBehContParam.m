% plotCorrCoeffBehContParam.m
%
% Function to plot output of getCorrCoeffBehContParam_cond_all()
% Heatmap (using MATLAB function) of mean correlation coefficient across
%  flies, where each one is its own file (assumes multiple files all have
%  same matrix of correlation coefficients)
% Select files through GUI
%
% INPUTS:
%   datDir - path to folder containing getCorrCoeffBehContParam_cond_all()
%     output files
%   cScale - [min max] correlation values for color scale
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 10/7/23 - HHY
%
% UPDATED:
%   10/7/23 - HHY
%
function plotCorrCoeffBehContParam(datDir, cScale)

    % select correlation files
    [corrFNames, corrPath] = uigetfile('*.mat', ...
        'Select output file', datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(corrFNames))
        numFlies = length(corrFNames);
    else
        numFlies = 1;
    end

    % initialize
    corrAllFlies = [];

    % loop through all files
    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(corrFNames))
            outName = corrFNames{i};
        else
            outName = corrFNames;
        end

        % load data from cond_bout file
        fullFilePath = [corrPath filesep outName];

        load(fullFilePath,'allCorr', 'allBehParams');

        % save correlation for this fly
        corrAllFlies = cat(3,corrAllFlies,allCorr);
    end

    % get mean correlation across flies
    meanCorr = mean(corrAllFlies,3);

    % turn allBehParams into labels
    labels = {};
    for i = 1:size(meanCorr,1)
        if (iscell(allBehParams(i).whichParams))
            if ~isempty(allBehParams(i).legs)
                labels{i} = [allBehParams(i).whichParams{1} ' ' ...
                    allBehParams(i).legs{1}];
            else
                labels{i} = [allBehParams(i).whichParams{1}];
            end
        else
            if ~isempty(allBehParams(i).legs)
                labels{i} = [allBehParams(i).whichParams ' ' ...
                    allBehParams(i).legs];
            else
                labels{i} = allBehParams(i).whichParams;
            end
        end
    end

    % generate heatmap figure
    figure;
    c = colormap('cool');
    heatmap(labels, labels, meanCorr, ...
        'Colormap',c,'ColorLimits',cScale);

end