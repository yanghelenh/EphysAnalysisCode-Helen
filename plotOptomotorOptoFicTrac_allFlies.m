% plotOptomotorOptoFicTrac_allFlies.m
%
% Function that takes outputs from extractOptomotorFicTracOptoCond_fly()
%  and plots the average optomotor response over time
% Plots mean across flies
% Has option to plot mean for each fly
% Select one duration, one FicTrac parameter, one or more NDs
% If more than 1 ND, different NDs plotted in different colors
% Select output files through GUI
% 
% INPUTS:
%   datDir - directory with output files
%   whichParam - which FicTrac parameter to plot
%   whichVel - which velocity condition to plot
%   whichNDs - which ND(s) to plot
%   plotIndiv - boolean for whether to plot means of individual flies
%   plotMean - boolean for whether to plot mean across flies
%   yScale - scale for plot, as [min max]
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED:
%   4/10/24 - HHY
%
% UPDATED:
%   4/10/24 - HHY
%
function plotOptomotorOptoFicTrac_allFlies(datDir, whichParam, whichVel, ...
    whichNDs, plotIndiv, yScale)

    % prompt user to select output files from 
    %  extractOptomotorFicTracOptoCond_fly()
    [outputFNames, outputPath] = uigetfile('*.mat', ...
        'Select extractOptomotorFicTracOptoCond_fly() files', ...
        datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(outputFNames))
        numFlies = length(outputFNames);
    else
        numFlies = 1;
    end

    % number of amplitude conditions
    numNDs = length(whichNDs);

    % preallocate - mean for each fly
    allFliesMeans = [];

    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];

        % load data
        load(outputFullPath, 'changeFictracOptoMean', 'vels', 'NDs', ...
            'trialTimes', 'optomotorDur', 'optoStartTime', 'optoEndTime');

        thisVelInd = find(vels == whichVel);
        thisTrialTimes = trialTimes{thisVelInd};

        % preallocate - mean for this fly
        thisFlyMeans = zeros(length(thisTrialTimes), numNDs);

        % loop through all NDs
        for j = 1:numNDs
            thisNDInd = find(NDs == whichNDs(j));

            % save mean for this ND
            thisFlyMeans(:,j) = ...
                changeFictracOptoMean(thisVelInd,thisNDInd).(whichParam);
        end

        % save means for this fly
        allFliesMeans = cat(3, allFliesMeans, thisFlyMeans);
    end

    % compute mean and SEM across all flies
    meanAllFlies = zeros(size(allFliesMeans,1),size(allFliesMeans, 2));
    SEMAllFlies = zeros(size(allFliesMeans,1),size(allFliesMeans, 2));
    for i = 1:size(allFliesMeans,1)
        for j = 1:size(allFliesMeans, 2)
            thisRow = allFliesMeans(i,j,:);
            thisMean = mean(thisRow(~isnan(thisRow)));
            thisStd = std(thisRow(~isnan(thisRow)));
            thisN = length(thisRow(~isnan(thisRow)));

            meanAllFlies(i,j) = thisMean;
            SEMAllFlies(i,j) = thisStd / sqrt(thisN);
        end
    end

    % initialize figure
    figure;
    c = colormap('lines');

    % initialize - handles for mean lines
    meanLineHdl = zeros(1,numNDs);
    % initialize - text for legend for each line
    legendStr = cell(1,numNDs);

    for j = 1:numNDs
        % plot indiv flies
        if plotIndiv
            plot(thisTrialTimes, squeeze(allFliesMeans(:,j,:)),...
                'LineWidth',0.5, 'Color', c(j,:));
        end

        hold on;

        % plot mean
        meanLineHdl(j) = plot(thisTrialTimes, meanAllFlies(:,j), ...
            'Color', c(j,:), 'LineWidth', 2);


        ylim(yScale);
        xScale = [thisTrialTimes(1) thisTrialTimes(end)];
        xlim(xScale)

        % plot lines for optomotor and opto stim
        % optomotor start and end, as dashed line
        line([0,0], yScale, 'Color','k', 'LineWidth', 1, ...
            'LineStyle', '--'); % start
        line([optomotorDur,optomotorDur], yScale, 'Color','k', ...
            'LineWidth', 1, 'LineStyle', '--'); % end
        % optogenetic stimulation start and end, as dotted line
        line([optoStartTime,optoStartTime], yScale, 'Color','k', ...
            'LineWidth', 1, 'LineStyle', ':'); % start
        line([optoEndTime,optoEndTime], yScale, 'Color','k', ...
            'LineWidth', 1, 'LineStyle', ':'); % end

        % generate legend string for this ND
        legendStr{j} = sprintf('ND=%.1f', whichNDs(j));
    end

    % plot legend
    legend(meanLineHdl,legendStr);

    % label axes
    ylabel(whichParam);
    xlabel('time (s)');

    % title string
    titleStr = sprintf('%s, Vel=%.1f', whichParam, whichVel);
    sgtitle(titleStr);

end