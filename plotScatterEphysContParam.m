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
%       'meanSEM' - mean +/- SEM, as line plus shading. Binned
%       'medMAD' - median +/- MAD, as line plus shadding. Binned
%   numBins - for 'binScatter' type plot, number of bins per dimension. For
%     'meanSEM' or 'medMAD', number of bins in x
%   xRange - range of x values to plot, or for 'meanSEM' and 'medMAD', to
%       bin across. As 2 element vector for start and end
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 8/8/23 - HHY
%
% UPDATED:
%   8/8/23 - HHY
%   8/11/23 - HHY - add 'meanSEM' and 'medMAD' plotting options
%
function plotScatterEphysContParam(datDir, plotType, numBins, xRange)

    % prompt user to select getCorrelationEphysContParam() file
    [corrEphysFName, corrEphysPath] = uigetfile('*.mat', ...
        'Select output file', datDir, 'MultiSelect', 'off');

    % load data from cond_bout file
    fullFilePath = [corrEphysPath filesep corrEphysFName];

    % load variables
    load(fullFilePath, 'ephysVals', 'behVals1D', 'linMdl', ...
        'ephysParam', 'behParams', 'legs', 'tDelay');

    % get boundaries of bins for 'meanSEM' and 'medMAD' type plots
    binSize = (xRange(2) - xRange(1)) / numBins;
    binEdges = xRange(1):binSize:xRange(2);
    binStarts = binEdges(1:(end-1));
    binEnds = binEdges(2:end);
    binMids = (binStarts + binEnds)/2;

    % generate figure
    figure;

    if (strcmpi(plotType, 'binScatter'))
        binscatter(behVals1D, ephysVals, numBins);
    elseif (strcmpi(plotType, 'scatter'))
        scatter(behVals1D, ephysVals, 'filled', ...
            'MarkerFaceAlpha', 0.03, 'MarkerEdgeAlpha', 0.03);
    elseif (strcmpi(plotType, 'linMdl'))
        plot(linMdl);
    elseif (strcmpi(plotType, 'meanSEM'))
        % preallocate vectors for mean and SEM
        allMeans = zeros(numBins,1);
        allSEMs = zeros(numBins,1);
        % loop through all bins, find mean and SEM
        for i = 1:numBins
            % get logical for which samples fall into this bin
            thisBinLog = (behVals1D >= binStarts(i)) & ...
                (behVals1D < binEnds(i));
            % get ephys values for these samples
            thisEphys = ephysVals(thisBinLog);
            % get mean and SEM for this bin
            allMeans(i) = mean(thisEphys);
            allSEMs(i) = std(thisEphys) / sqrt(length(thisEphys));
        end

        % plot
        plot_err_patch_v2(binMids, allMeans, allSEMs,...
            [0 0.4470 0.7410],[0.3010 0.7450 0.9330]);
    elseif (strcmpi(plotType, 'medMAD'))
        % preallocate vectors for median and MAD
        allMeds = zeros(numBins,1);
        allMADs = zeros(numBins,1);
        % loop through all bins, find mean and SEM
        for i = 1:numBins
            % get logical for which samples fall into this bin
            thisBinLog = (behVals1D >= binStarts(i)) & ...
                (behVals1D < binEnds(i));
            % get ephys values for these samples
            thisEphys = ephysVals(thisBinLog);
            % get mean and SEM for this bin
            allMeds(i) = median(thisEphys);
            allMADs(i) = mad(thisEphys,1);
        end

        % plot
        plot_err_patch_v2(binMids, allMeds, allMADs,...
            [0 0.4470 0.7410],[0.3010 0.7450 0.9330]);
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