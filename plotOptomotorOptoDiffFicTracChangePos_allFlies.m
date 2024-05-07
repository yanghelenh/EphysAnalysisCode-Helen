% plotOptomotorOptoDiffFicTracChangePos_allFlies.m
%
% Function that takes outputs from extractOptomotorFicTracOptoCond_fly()
%  and plots the change in FicTrac position parameters. Specifically, for 
%  each fly, plots the average across trials of the difference in 
%  cumulative position/angle between two time periods, each defined by a 
%  start time point and an end time point relative to the optomotor 
%  stimulus start. For each fly, option for having the mean of each 
%  optogenetic condition subtracted from the mean of the no stimulation 
%  condition, and the mean across flies is plotted.
% Select output files through GUI
%
% INPUTS:
%   datDir - directory with output files
%   whichParam - which FicTrac parameter to plot
%   whichVel - which optomotor stimulus velocity to plot
%   whichNDs - which NDs to plot
%   plotTimes1 - [start time pt, end time pt] vector of time points relative
%       to the optomotor stimulus start over which to calculate difference
%       in cumulative angle/position, first time period
%   plotTimes2 - [start time pt, end time pt] vector of time points relative
%       to the optomotor stimulus start over which to calculate difference
%       in cumulative angle/position, second time period
%   plotDiff - boolean for whether to plot difference from no stim
%       condition
%   yScale - 2 element vector for y-axis limits
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 4/30/24 - HHY
%
% UPDATED:
%   4/30/24 - HHY
%
function allFliesMeans = plotOptomotorOptoDiffFicTracChangePos_allFlies(datDir, whichParam, ...
    whichVel, whichNDs, plotTimes1, plotTimes2, plotDiff, yScale)

    % prompt user to select output files from extractFicTracOpto_fly()
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

    % preallocate - mean and SEM for each fly
    allFliesMeans = zeros(numFlies, numNDs);
    allFliesSEMs = zeros(numFlies, numNDs);

    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];

        % load data
        load(outputFullPath, 'fictracOpto', 'vels', 'trialTimes');

        thisVelInd = find(vels == whichVel);
        thisTrialTimes = trialTimes{thisVelInd};

        % get indices for start and end time points of comparison, first
        %  time period
        startInd1 = find(plotTimes1(1) <= thisTrialTimes, 1, 'first');
        endInd1 = find(plotTimes1(2) >= thisTrialTimes, 1, 'last');

        % get indices for start and end time points of comparison, second
        %  time period
        startInd2 = find(plotTimes2(1) <= thisTrialTimes, 1, 'first');
        endInd2 = find(plotTimes2(2) >= thisTrialTimes, 1, 'last');

        % preallocate - mean and SEM for this fly
        thisFlyMeans = zeros(1, numNDs);
        thisFlySEMs = zeros(1, numNDs);

        % get difference in cumulative position/angle between start and end
        %  points for no stim condition (ND = -1)
        noStimLog = fictracOpto(thisVelInd).whichND == -1;
        noStimMean = mean((fictracOpto(thisVelInd).(whichParam)(endInd2, noStimLog) - ...
                fictracOpto(thisVelInd).(whichParam)(startInd2, noStimLog)) - ...
                (fictracOpto(thisVelInd).(whichParam)(endInd1, noStimLog) - ...
                fictracOpto(thisVelInd).(whichParam)(startInd1, noStimLog))); 


        % loop through all NDs
        for j = 1:numNDs
            thisNDLog = fictracOpto(thisVelInd).whichND == whichNDs(j);

            % get difference in cumulative position/angle between start and
            %  end time points
            diffAllTrials = (fictracOpto(thisVelInd).(whichParam)(endInd2, thisNDLog) - ...
                fictracOpto(thisVelInd).(whichParam)(startInd2, thisNDLog)) - ...
                (fictracOpto(thisVelInd).(whichParam)(endInd1, thisNDLog) - ...
                fictracOpto(thisVelInd).(whichParam)(startInd1, thisNDLog));

            % save mean and SEM for this ND
            if plotDiff
                thisFlyMeans(j) = mean(diffAllTrials) - noStimMean;
            else
                thisFlyMeans(j) = mean(diffAllTrials);
            end
            thisFlySEMs(j) = std(diffAllTrials) / sqrt(length(diffAllTrials));
        end

        % save means and SEMs for this fly
        allFliesMeans(i,:) = thisFlyMeans;
        allFliesSEMs(i,:) = thisFlySEMs;
    end

    % get mean and SEM across flies
    for i = 1:numNDs
        thisNDVals = allFliesMeans(:,i);
        meanAllFlies(i) = mean(thisNDVals(~isnan(thisNDVals)));
        semAllFlies(i) = std(thisNDVals(~isnan(thisNDVals))) / ...
            sqrt(length(thisNDVals(~isnan(thisNDVals))));
%     meanAllFlies = mean(allFliesMeans,1);
%     semAllFlies = std(allFliesMeans,[],1) / sqrt(numFlies);
    end


    % plot
    figure;
    c = colormap('lines');

    % x axis 
    xVec = 1:numNDs;

    % x labels
    tickLabels = cell(numNDs,1);
    for i = 1:numNDs
        tickLabels{i} = sprintf('ND=%.1f', whichNDs(i));
    end

    % plot individual flies
    for i = 1:numFlies
        plot(xVec, allFliesMeans(i,:), ...
            'Marker', '.','LineWidth',0.5, 'Color', c(1,:));
    
        hold on;
    end

    % plot mean across flies
    plot(xVec, meanAllFlies, ...
        'Marker','_','LineWidth', 2, 'Color', c(2,:));

    ylim(yScale);
    xScale = xlim;
    xScale(1) = xScale(1) - (0.3 * (xVec(end)-xVec(1)));
    xScale(2) = xScale(2) + (0.3 * (xVec(end)-xVec(1)));
    xlim(xScale);

    % line at 0
    line(xScale,[0,0],'Color','k', 'LineWidth', 2);

    % tick labels
    xticks(xVec);
    xticklabels(tickLabels);

    xlabel('ND');
    ylabel(whichParam);
    title(sprintf('%s, Vel=%.1f, change (t=%.1f to %.1f) - (t=%.1f to %.1f)',...
        whichParam, whichVel, plotTimes2(1), plotTimes2(2), ...
        plotTimes1(1), plotTimes1(2)));
end