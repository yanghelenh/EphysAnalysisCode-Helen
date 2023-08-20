% plotOptoFicTrac_allFlies.m
%
% Function that takes output files from extractFicTracOpto_fly() and
%  generates a plot of the specified FicTrac param. Plots 
%  mean for each fly as well as across flies. User selects one dur, one ND,
%  one FicTrac parameter, and whether to plot raw values or change in the
%  parameter
% Select output files through GUI
%
% INPUTS:
%   datDir - directory with output files
%   whichParam - which FicTrac parameter to plot
%   whichDur - which duration conditions to plot
%   whichND - which NDs to plot
%   plotChange - boolean for whether to plot change or raw values (true for
%       change)
%   yScale - scale for plots, as [min max]
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 8/18/23 - HHY
%
% UPDATED:
%   8/18/23 - HHY
%
function plotOptoFicTrac_allFlies(datDir, whichParam, whichDur,...
    whichND, plotChange, yScale)

    % prompt user to select output files from extractFicTracOpto_fly()
    [outputFNames, outputPath] = uigetfile('*.mat', 'Select Step Param files', ...
        datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(outputFNames))
        numFlies = length(outputFNames);
    else
        numFlies = 1;
    end

    % initialize
    allFliesMeans = [];
    allFliesSEMs = [];    


    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];

        % load data
        load(outputFullPath, 'fictracOptoMean', 'fictracOptoSEM', ...
            'changeFictracOptoMean', 'changeFictracOptoSEM',...
            'durs', 'NDs', 'trialTimes');

        thisDurInd = find(durs == whichDur);
        thisNDInd = find(NDs == whichND);
        thisTrialTimes = trialTimes{thisDurInd};

        if (plotChange)

        thisMean = changeFictracOptoMean(thisDurInd, thisNDInd).(whichParam);
        thisSEM = changeFictracOptoSEM(thisDurInd, thisNDInd).(whichParam);

        else
            thisMean = fictracOptoMean(thisDurInd, thisNDInd).(whichParam);
            thisSEM = fictracOptoSEM(thisDurInd, thisNDInd).(whichParam);
        end


        % save means and SEMs for this fly
        allFliesMeans = cat(2, allFliesMeans, thisMean);
        allFliesSEMs = cat(2, allFliesSEMs, thisSEM);
    end

    % compute mean, median, SEM across all flies
    meanAllFlies = zeros(size(allFliesMeans,1),1);
    medianAllFlies = zeros(size(allFliesMeans,1),1);
    SEMAllFlies = zeros(size(allFliesMeans,1),1);

    for i = 1:size(allFliesMeans,1)
        thisRow = allFliesMeans(i,:);
        thisMean = mean(thisRow(~isnan(thisRow)));
        thisMedian = median(thisRow(~isnan(thisRow)));
        thisStd = std(thisRow(~isnan(thisRow)));
        thisN = length(thisRow(~isnan(thisRow)));

        meanAllFlies(i) = thisMean;
        medianAllFlies(i) = thisMedian;
        SEMAllFlies(i) = thisStd / sqrt(thisN);
    end


    % initialize figure
    figure;
    c = colormap('lines');

    hold on;

    % plot indiv flies
    plot(thisTrialTimes, allFliesMeans,...
        'LineWidth',0.5, 'Color', c(1,:));

    % plot mean
    plot(thisTrialTimes, meanAllFlies, 'Color', c(2,:), ...
        'LineWidth', 2);


    ylim(yScale);
    xScale = [thisTrialTimes(1) thisTrialTimes(end)];
    xlim(xScale)

    % plot lines for opto
    line([0,0], yScale, 'Color','k', 'LineWidth', 1);
    line([whichDur whichDur], yScale, 'Color','k', 'LineWidth', 1);

    titleStr = sprintf('%s, ND=%.1f, dur=%.1f', whichParam, whichND, ...
        whichDur);
    sgtitle(titleStr);
end
