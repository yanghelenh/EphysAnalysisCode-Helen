% plotEphysBehXCorr.m
%
% Function to plot cross-correlation of ephys parameter (output of
%  getXCorrEphysContParam_cell()). Plots the cross-correlation itself as 
%  well as the peak time, across conditions. Plots one line/pt per cell and
%  mean +/- SEM across cells
% Select output files through GUI
%
% INPUTS:
%   datDir - folder containing getXCorrEphysContParam_cell() output files
%   numCond - number of conditions
%   condNames - cell array of length numCond, string for name of each 
%       condition   
%   yScaleXC - limits of y axis for plot of cross-correlation
%   yScalePeakT - limits of y axis for plot of peak times
%
% OUTPUTS:
%   none, but generates plots
%
% CREATED: 9/1/23 - HHY
%
% UPDATED:
%   9/1/23 - HHY
%   9/6/23 - HHY - add time of center of mass plot
%
function allCellsPeakT = plotEphysBehXCorr(datDir, numCond, condNames, yScaleXC, ...
    yScalePeakT)

    % preallocate 
    % cell array where each element will be numCells x numTPts matrix of
    %  cross-correlations
    allCellsXC = cell(numCond,1);
    % cell array where each element will be numCells length vector of peak
    %  times
    allCellsPeakT = cell(numCond,1);
    % cell array where each element will be numCells length vector of
    %  center of mass times
    allCellsCOMtime = cell(numCond,1);

    % cell array where each element will be 1 x numTPts vector of mean/SEM
    %  across flies
    allCondMeanXC = cell(numCond,1);
    allCondSEMsXC = cell(numCond,1);
    allCondMeanPeakT = zeros(numCond,1);
    allCondSEMsPeakT = zeros(numCond,1);
    allCondMeanCOMtime = zeros(numCond,1);
    allCondSEMsCOMtime = zeros(numCond,1);
    allCondLagsT = cell(numCond,1);

    % loop across number of conditions, get data files for each cell
    % compute mean across flies
    for i = 1:numCond
        [outputFNames, outputPath] = uigetfile('*.mat', ...
            'Select getXCorrEphysContParam_cell() files', ...
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
        thisCondXC = [];
        thisCondPeakT = zeros(numCells,1);
        % time of center of mass
        thisCondCOMtime = zeros(numCells,1);

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
            load(outputFullPath, 'xCorr', 'lagsT', 'peakT');

            % save this peak time
            thisCondPeakT(j) = peakT;

            % save this cross-corr
            thisCondXC = [thisCondXC, xCorr];

            % calculate time of center of mass (time at which half max in
            %  cumulative sum occurs)
            % remove NaNs
            nanLog = isnan(xCorr);
            noNanXCorr = xCorr;
            noNanXCorr(nanLog) = [];
            noNanLagsT = lagsT;
            noNanLagsT(nanLog) = [];

            % minimum subtracted xCorr, so not dealing with negatives
            minSubXCorr = noNanXCorr - min(noNanXCorr);
            % cumulative sum
            xCorrCumSum = cumsum(minSubXCorr);
            % half max of cumulative sum
            halfMax = xCorrCumSum(end)/2;
            % find indices that bracket half-max value
            ind1 = find(halfMax >= xCorrCumSum, 1, 'last');
            ind2 = find(halfMax < xCorrCumSum, 1, 'first');

            % linear interpolation b/w 2 indices that bracket half-max
            %  value, best estimate of time of half max
            tHalfMax = interp1([xCorrCumSum(ind1), xCorrCumSum(ind2)], ...
                [noNanLagsT(ind1), noNanLagsT(ind2)], halfMax, 'linear');

            % save this time of center of mass
            thisCondCOMtime(j) = tHalfMax;

        end

        % crosscorr, peakT, lagsT across cond
        allCellsXC{i} = thisCondXC;
        allCellsPeakT{i} = thisCondPeakT;
        allCellsCOMtime{i} = thisCondCOMtime;
        allCondLagsT{i} = lagsT;

        % get mean and SEM
        thisMean = zeros(size(thisCondXC,1),1);
        thisSEM = zeros(size(thisCondXC,1),1);
        for j = 1:size(thisCondXC,1)
            thisRow = thisCondXC(j,:);
            thisRow = thisRow(~isnan(thisRow));
            thisMean(j) = mean(thisRow);
            thisSEM(j) = std(thisRow) / sqrt(length(thisRow));
        end

        allCondMeanXC{i} = thisMean;
        allCondSEMsXC{i} = thisSEM;

%         allCondMeanXC{i} = mean(thisCondXC,2);
%         allCondSEMsXC{i} = std(thisCondXC,[],2) / sqrt(numCells);
        allCondMeanPeakT(i) = mean(thisCondPeakT);
        allCondSEMsPeakT(i) = std(thisCondPeakT) / sqrt(numCells);

        allCondMeanCOMtime(i) = mean(thisCondCOMtime);
        allCondSEMsCOMtime(i) = std(thisCondCOMtime) / sqrt(numCells);
    end

    % plot autocorr
    figure;
    c = colormap('lines');
    hold on;

    legendInd = [];
    % loop over all conditions
    for i = 1:numCond
        % plot individual flies
        plot(allCondLagsT{i}, allCellsXC{i}, ...
            'LineWidth',0.5, 'Color', c(i,:)');

         % plot mean across flies
        legendInd(i) = plot(allCondLagsT{i}, allCondMeanXC{i}, ...
            'LineWidth',2, 'Color', c(i,:));
    end

    % line at t = 0
    line([0 0], yScaleXC, 'LineWidth', 1, 'Color', 'k');

    % label axes
    ylim(yScaleXC);
    ylabel('Cross-correlation');

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


    % plot center of mass times
    figure;
    hold on;

    xVec = 1:numCond;

    for i = 1:numCond
        % plot individual cells as points
        numCells = length(allCellsCOMtime{i});
        % get repmat xVec for this number of cells
        repXVec = repmat(xVec(i),numCells,1);

        plot(repXVec, allCellsCOMtime{i}, ...
            '.','Marker','.','LineWidth',0.5, 'Color', c(1,:));

        % plot mean +/- SEM across flies
        errorbar(xVec(i),allCondMeanCOMtime(i),allCondSEMsCOMtime(i),...
            '_','LineWidth', 2, 'CapSize', 0, 'Color', c(2,:));
    end

    xScale = xlim;
    xScale(1) = xScale(1) - (0.5 * (xVec(end)-xVec(1)));
    xScale(2) = xScale(2) + (0.5 * (xVec(end)-xVec(1)));
    xlim(xScale);

    xticks(xVec);
    xticklabels(condNames);

    ylim(yScalePeakT);

    ylabel('Time of Center of Mass (s)');
end