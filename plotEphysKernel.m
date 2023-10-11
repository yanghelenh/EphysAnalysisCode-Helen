% plotEphysKernel.m
%
% Function to plot kernel relating of spike rate to behavioral param 
%  (output of computeKernelFromAutoXcorr()). Plots the autocorrelation itself as well
%  as the FWHM, across conditions. Plots one line/pt per cell and mean +/-
%  SEM across cells
% Select output files through GUI
%
% INPUTS:
%   datDir - folder containing computeKernelFromAutoXcorr() output files
%   numCond - number of conditions
%   condNames - cell array of length numCond, string for name of each 
%       condition   
%   yScaleK - limits of y axis for plot of kernel
%   yScalePeakT - limits of y axis for plot of peak times
%
% OUTPUTS:
%   none, but generates plots
%
% CREATED: 10/5/23 - HHY
%
% UPDATED:
%   10/5/23 - HHY
%
function allCellsPeakT = plotEphysKernel(datDir, numCond, condNames, ...
    yScaleK, yScalePeakT)

    % preallocate 
    % cell array where each element will be numCells x numTPts matrix of
    %  kernels
    allCellsK = cell(numCond,1);
    % cell array where each element will be numCells length vector of peak
    %  times
    allCellsPeakT = cell(numCond,1);

    % cell array where each element will be 1 x numTPts vector of mean/SEM
    %  across flies
    allCondMeanK = cell(numCond,1);
    allCondSEMsK = cell(numCond,1);
    allCondMeanPeakT = zeros(numCond,1);
    allCondSEMsPeakT = zeros(numCond,1);
    allCondLagsT = cell(numCond,1);

    % loop across number of conditions, get data files for each cell
    % compute mean across flies
    for i = 1:numCond
        [outputFNames, outputPath] = uigetfile('*.mat', ...
            'Select computeKernelFromAutoXcorr() files', ...
            datDir, 'MultiSelect', 'on');

        % if only 1 file selected, not cell array; make sure loop still
        %  works 
        % num flies is number of cells
        if (iscell(outputFNames))
            numCells = length(outputFNames);
        else
            numCells = 1;
        end

        % preallocate
        thisCondK = [];
        thisCondPeakT = zeros(numCells,1);

        % loop through all cells
        for j = 1:numCells
            % handle whether it's a cell array or not
            if (iscell(outputFNames))
                outName = outputFNames{j};
            else
                outName = outputFNames;
            end

            outputFullPath = [outputPath outName];

            % load data file
            load(outputFullPath, 'kernel', 'lagsK', 'peakT');

            % save this peak time
            thisCondPeakT(j) = peakT;

            % save this kernel
            thisCondK = [thisCondK, kernel];
        end

        % crosscorr, peakT, lagsT across cond
        allCellsK{i} = thisCondK;
        allCellsPeakT{i} = thisCondPeakT;
        allCondLagsT{i} = lagsK;

        % get mean and SEM
        thisMean = zeros(size(thisCondK,1),1);
        thisSEM = zeros(size(thisCondK,1),1);
        for j = 1:size(thisCondK,1)
            thisRow = thisCondK(j,:);
            thisRow = thisRow(~isnan(thisRow));
            thisMean(j) = mean(thisRow);
            thisSEM(j) = std(thisRow) / sqrt(length(thisRow));
        end

        allCondMeanK{i} = thisMean;
        allCondSEMsK{i} = thisSEM;

%         allCondMeanXC{i} = mean(thisCondXC,2);
%         allCondSEMsXC{i} = std(thisCondXC,[],2) / sqrt(numCells);
        allCondMeanPeakT(i) = mean(thisCondPeakT);
        allCondSEMsPeakT(i) = std(thisCondPeakT) / sqrt(numCells);
    end

    % plot kernel
    figure;
    c = colormap('lines');
    hold on;

    legendInd = [];
    % loop over all conditions
    for i = 1:numCond
        % plot individual flies
        plot(allCondLagsT{i}, allCellsK{i}, ...
            'LineWidth',0.5, 'Color', c(i,:)');

         % plot mean across flies
        legendInd(i) = plot(allCondLagsT{i}, allCondMeanK{i}, ...
            'LineWidth',2, 'Color', c(i,:));
    end

%     % line at t = 0
    line([0 0], yScaleK, 'LineWidth', 1, 'Color', 'k');

    % label axes
    ylim(yScaleK);
    ylabel('Kernel');

    xlim([allCondLagsT{1}(1) allCondLagsT{1}(end)]);
    xlabel('Time (s)');

    legend(legendInd,condNames);


    % plot peak times
    figure;
    hold on;

    xVec = 1:numCond;

    for i = 1:numCond
        % plot individual cells as points
        numCells = length(allCellsPeakT{i});
        % get repmat xVec for this number of cells
        repXVec = repmat(xVec(i),numCells,1);

        plot(repXVec, allCellsPeakT{i}, ...
            '.','Marker','.','LineWidth',0.5, 'Color', c(1,:));

        % plot mean +/- SEM across flies
        errorbar(xVec(i),allCondMeanPeakT(i),allCondSEMsPeakT(i),...
            '_','LineWidth', 2, 'CapSize', 0, 'Color', c(2,:));
    end

    xScale = xlim;
    xScale(1) = xScale(1) - (0.5 * (xVec(end)-xVec(1)));
    xScale(2) = xScale(2) + (0.5 * (xVec(end)-xVec(1)));
    xlim(xScale);

    xticks(xVec);
    xticklabels(condNames);

    ylim(yScalePeakT);

    ylabel('Peak Time (s)');
end