% genCondBinnedSpikerateFictracPlot.m
%
% Function that generates plot of spike rate vs. FicTrac variable 1,
%  conditioned on the values of FicTrac variable 2. Plots points with error
%  bars. FicTrac variable 1 is binned finely, each point represents 1 bin.
%  FicTrac variable 2 is binned coarsely, each bin corresponding to another
%  line. Spike rate and FicTrac variables must be on same time scale.
% Does not introduce any temporal offset between spike rate and FicTrac
%  variables. If desired, do that in wrapper function for this function 
%  (plotCondBinnedSpikerateFicTrac()).
%
% INPUTS:
%   spikeRate - vector of spike rate values, length n
%   ftVar1 - vector of values of FicTrac variable 1, the one binned finely,
%       x-axis variable, length n
%   ftVar2 - vector of values of FicTrac variable 2, the one binned
%       coarsely, length n
%   ftVar1Bins - vector of length 3: [lower lim, upper lim, number bins]
%       for ftVar1
%   ftVar2Bins - matrix of size m x 2, where each row is [lower lim, upper
%       lim] for one bin; number of bins = m; will produce m separate lines
%   yScale - vector of length 2 for yScale of plot [lower upper]
%   xAxisName - string for x-axis label
%   yAxisName - string for y-axis label
%   ftVar1Name - string, name of ftVar1
%   ftVar2Name - string, name of ftVar2
%   ftVar2Units - string, units for ftVar2
%
% OUTPUTS:
%   f - handle to figure
%   meanSpikeRate - matrix of mean spike rate, size m x num bins ftVar1
%   stdDevSpikeRate - matrix of standard deviations of spike rate, size m x
%       num bins ftVar 1
%   binMids - vector of length number of bins of ftVar1, each value is
%       midpoint of bin
%
% CREATED: 10/26/22 - HHY
%
% UPDATED:
%   10/26/22 - HHY
%
function [f, meanSpikeRate, stdDevSpikeRate, semSpikeRate, binMids] = ...
    genCondBinnedSpikerateFictracPlot(spikeRate, ftVar1, ftVar2, ...
    ftVar1Bins, ftVar2Bins, yScale, xAxisName, yAxisName, ftVar1Name, ...
    ftVar2Name, ftVar2Units)

    % get start, mid, end points for each bin of ftVar1
    % width of bin
    binWidth = (ftVar1Bins(2) - ftVar1Bins(1)) / ftVar1Bins(3);

    binStarts = ftVar1Bins(1):binWidth:(ftVar1Bins(2)-binWidth);
    binEnds = (ftVar1Bins(1) + binWidth):binWidth:ftVar1Bins(2);

    % bin midpoints (can average bin starts and ends b/c linear)
    binMids = (binStarts + binEnds) / 2;

    % number of bins for ftVar2
    numBinsFT2 = size(ftVar2Bins,1);

    % preallocate cell array for categorizing each spikeRate point
    % each row as different condition of ftVar2
    % each column as different bin of ftVar1
    % will populate vector in each cell
    catSpikeRate = cell(numBinsFT2, ftVar1Bins(3));

    % loop through all points of spikeRate, put each one in appropriate bin
    for i = 1:length(spikeRate)

        var2BinInd = []; % rezero each time
        % determine which bin for ftVar2
        for j = 1:numBinsFT2
            if ((ftVar2(i) >= ftVar2Bins(j,1)) && (ftVar2(i) < ftVar2Bins(j,2)))
                var2BinInd = j;
            end
        end

        % determine which bin for ftVar1
        % find index of bin, going from both directions of start and end of
        %  bin
        whichBinStart = find(ftVar1(i) >= binStarts, 1, 'last');
        whichBinEnd = find(ftVar1(i) < binEnds, 1, 'first');

        % if either is empty, ftVar1 point exceeds x limits
        if (~isempty(whichBinStart) && ~isempty(whichBinEnd))
            var1BinInd = whichBinStart;
        else
            var1BinInd = [];
        end

        % if this point belongs to a bin for both var 1 and var 2
        if (~isempty(var2BinInd) && ~isempty(var1BinInd))
            % append this spike rate to appropriate vector of cell array
            catSpikeRate{var2BinInd, var1BinInd} = [...
                catSpikeRate{var2BinInd, var1BinInd}; spikeRate(i)];
        end
    end

    % convert vectors of of individual values for each bin into mean,
    %  standard deviation, and counts for each bin

    % preallocate
    meanSpikeRate = zeros(numBinsFT2, ftVar1Bins(3));
    stdDevSpikeRate = zeros(numBinsFT2, ftVar1Bins(3));
    countsSpikeRate = zeros(numBinsFT2, ftVar1Bins(3));
    semSpikeRate = zeros(numBinsFT2, ftVar1Bins(3));

    % loop over all bins, fill in mean, std dev, counts
    for i = 1:numBinsFT2
        for j = 1:ftVar1Bins(3)
            meanSpikeRate(i,j) = mean(catSpikeRate{i,j});
            stdDevSpikeRate(i,j) = std(catSpikeRate{i,j});
            countsSpikeRate(i,j) = length(catSpikeRate{i,j});
            semSpikeRate(i,j) = stdDevSpikeRate(i,j) / sqrt(countsSpikeRate(i,j));
        end
    end


    % PLOTTING %

    f = figure;

    c = colormap('lines');
    cs = (c + 0.2) ./ repmat(max(c+0.2,[],2),1,3);

    legendStr = cell(1,numBinsFT2);
    legendHdl = zeros(1,numBinsFT2);

    for i = 1:numBinsFT2

%         errorbar(binMids,meanSpikeRate(i,:),stdDevSpikeRate(i,:),...
%             'Marker','x', 'Color', c(i,:), 'LineWidth', 1.5);

%         errorbar(binMids,meanSpikeRate(i,:),semSpikeRate(i,:),...
%             'Marker','x', 'Color', c(i,:), 'LineWidth', 1.5);
%         plot(binMids, meanSpikeRate(i,:), 'Marker','x', 'Color', c(i,:), ...
%             'LineWidth', 1.5);

        valInd = ~isnan(meanSpikeRate(i,:));

        legendHdl(i) = plot_err_patch_v2(binMids(valInd), meanSpikeRate(i,valInd),semSpikeRate(i,valInd), ...
            c(i,:), cs(i,:));
%         legendHdl(i) = plot_err_patch_v2(binMids, meanSpikeRate(i,:),stdDevSpikeRate(i,:), ...
%             c(i,:), cs(i,:));

        hold on;
    
        legendStr{i} = sprintf('%.1f to %.1f %s', ftVar2Bins(i,1), ...
            ftVar2Bins(i,2), ftVar2Units);
    end

    xlim([binStarts(1)-binWidth, binEnds(end) + binWidth]);
    ylim(yScale);

    xlabel(xAxisName);
    ylabel(yAxisName);

    title(sprintf('Spike Rate vs. %s, conditioned on %s', ftVar1Name, ...
        ftVar2Name));

    legend(legendHdl, legendStr);
end