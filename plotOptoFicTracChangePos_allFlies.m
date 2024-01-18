% plotOptoFicTracChangePos_allFlies.m
%
% Function that plots effect of opto stim on FicTrac position
%  parameters. Specifically, for each fly, plots the average across trials
%  of the difference in cumulative position/angle between a start time 
%  point and an end time point relative to the opto stim start. For
%  each fly, the mean of each opto condition is subtracted from the mean of
%  the no stimulation condition, and the mean across flies is plotted.
% Operates on the output of extractFicTracOpto_fly(). Select output files
%  through GUI
%
% INPUTS:
%   datDir - directory with output files
%   whichParam - which FicTrac parameter to plot
%   whichND - which NDs to plot
%   plotTimes - [start time pt, end time pt] vector of time points relative
%       to current injection start over which to calculate difference in
%       cumulative angle/position
%   yScale - 2 element vector for y-axis limits
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 8/23/23 - HHY
%
% UPDATED:
%   8/23/23 - HHY
%
function allFliesMeans = plotOptoFicTracChangePos_allFlies(datDir, whichParam, ...
    whichND, plotTimes, yScale)

    % prompt user to select output files from extractFicTracOpto_fly()
    [outputFNames, outputPath] = uigetfile('*.mat', 'Select FicTracOpto files', ...
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
    numNDs = length(whichND);

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
        load(outputFullPath, 'fictracOpto', 'trialTimes');

        % get indices for start and end time points of comparison
        startInd = find(plotTimes(1) <= trialTimes{1}, 1, 'first');
        endInd = find(plotTimes(2) >= trialTimes{1}, 1, 'last');

        % preallocate - mean and SEM for this fly
        thisFlyMeans = zeros(1, numNDs);
        thisFlySEMs = zeros(1, numNDs);

        % get difference in cumulative position/angle between start and end
        %  points for no stim condition (ND = -1)
        noStimLog = fictracOpto.whichND == -1;
        noStimMean = mean(fictracOpto.(whichParam)(endInd, noStimLog) - ...
                fictracOpto.(whichParam)(startInd, noStimLog));


        % loop through all NDs
        for j = 1:numNDs
            thisNDLog = fictracOpto.whichND == whichND(j);

            % get difference in cumulative position/angle between start and
            %  end time points
            diffAllTrials = (fictracOpto.(whichParam)(endInd, thisNDLog) - ...
                fictracOpto.(whichParam)(startInd, thisNDLog));

            % save mean and SEM for this ND
            thisFlyMeans(j) = mean(diffAllTrials) - noStimMean;
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
        tickLabels{i} = sprintf('ND=%.1f', whichND(i));
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
    title(sprintf('%s difference t=%.1f - t=%.1f', whichParam, ...
        plotTimes(2), plotTimes(1)));
end