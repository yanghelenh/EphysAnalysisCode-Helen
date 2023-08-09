% plotScatterEphysContParam.m
%
% Function to plot output of getCorrelationEphysContParam() as scatterplot
%  binned scatterplot, or linear model
%
% INPUTS:
%   datDir - path to folder containing savePCAinterpStepParam() output files
%   plotType - string for type of plot to make
%       'binScatter' - binned scatterplot, counts
%       'scatter' - regular scatterplot
%       'linMdl' - linear model
%   numBins - for 'binScatter' type plot, number of bins per dimension
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 8/8/23 - HHY
%
% UPDATED:
%   8/8/23 - HHY
%
function plotScatterEphysContParam(datDir, plotType, numBins)

    % prompt user to select getCorrelationEphysContParam() file
    [corrEphysFName, corrEphysPath] = uigetfile('*.mat', ...
        'Select output file', datDir, 'MultiSelect', 'off');

    % load data from cond_bout file
    fullFilePath = [corrEphysPath filesep corrEphysFName];

    % load variables
    load(fullFilePath, 'ephysVals', 'behVals1D', 'linMdl', ...
        'ephysParam', 'behParams', 'legs', 'tDelay');

    % generate figure
    figure;

    if (strcmpi(plotType, 'binScatter'))
        binscatter(behVals1D, ephysVals, numBins);
    elseif (strcmpi(plotType, 'scatter'))
        scatter(behVals1D, ephysVals, 'filled', ...
            'MarkerFaceAlpha', 0.03, 'MarkerEdgeAlpha', 0.03);
    elseif (strcmpi(plotType, 'linMdl'))
        plot(linMdl);
    end

    % ephys param label
    yLblStr = sprintf('%s, delay = %d ms', ephysParam, tDelay * 1000);
    ylabel(yLblStr);

    xLblStr = [];
    if iscell(behParams)
        for i = 1:length(behParams)
            if ~i==1
                xLblStr = [xLblStr ', ' behParams{i} ' ' legs{i}];
            else
                xLblStr = [behParams{i} ' ' legs{i}];
            end
        end
    else
        xLblStr = [behParams ' ' legs];
    end
    xlabel(xLblStr);
end