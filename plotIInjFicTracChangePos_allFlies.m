% plotIInjFicTracChangePos_allFlies.m
%
% Function that plots effect of current injection on FicTrac position
%  parameters. Specifically, for each fly, plots the average across trials
%  of the difference in cumulative position/angle between a start time 
%  point and an end time point relative to the current injection start. For
%  each fly, the mean of each IInj condition is subtracted from the mean of
%  the no stimulation condition, and the mean across flies is plotted.
% Operates on the output of extractFicTracIInj_fly(). Select output files
%  through GUI
%
% INPUTS:
%   datDir - directory with output files
%   whichParam - which FicTrac parameter to plot
%   whichAmp - which amps to plot
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
function allFliesMeans = plotIInjFicTracChangePos_allFlies(datDir, whichParam, ...
    whichAmp, plotTimes, yScale)

    % prompt user to select output files from extractFicTracIInj_fly()
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

    % number of amplitude conditions
    numAmps = length(whichAmp);

    % preallocate - mean and SEM for each fly
    allFliesMeans = zeros(numFlies, numAmps);
    allFliesSEMs = zeros(numFlies, numAmps);

    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];

        % load data
        load(outputFullPath, 'fictracIInj', 'trialTimes');

        % get indices for start and end time points of comparison
        startInd = find(plotTimes(1) <= trialTimes{1}, 1, 'first');
        endInd = find(plotTimes(2) >= trialTimes{1}, 1, 'last');

        % preallocate - mean and SEM for this fly
        thisFlyMeans = zeros(1, numAmps);
        thisFlySEMs = zeros(1, numAmps);

        % get difference in cumulative position/angle between start and end
        %  points for no stim condition (amp = 0)
        noStimLog = fictracIInj.whichAmp == 0;
        noStimMean = mean(fictracIInj.(whichParam)(endInd, noStimLog) - ...
                fictracIInj.(whichParam)(startInd, noStimLog));


        % loop through all amps
        for j = 1:numAmps
            thisAmpLog = fictracIInj.whichAmp == whichAmp(j);

            % get difference in cumulative position/angle between start and
            %  end time points
            diffAllTrials = (fictracIInj.(whichParam)(endInd, thisAmpLog) - ...
                fictracIInj.(whichParam)(startInd, thisAmpLog));

            % save mean and SEM for this amp
            thisFlyMeans(j) = mean(diffAllTrials) - noStimMean;
            thisFlySEMs(j) = std(diffAllTrials) / sqrt(length(diffAllTrials));
        end

        % save means and SEMs for this fly
        allFliesMeans(i,:) = thisFlyMeans;
        allFliesSEMs(i,:) = thisFlySEMs;
    end

    % get mean and SEM across flies
    meanAllFlies = mean(allFliesMeans,1);
    semAllFlies = std(allFliesMeans,[],1) / sqrt(numFlies);


    % plot
    figure;
    c = colormap('lines');

    % x axis 
    xVec = 1:numAmps;

    % x labels
    tickLabels = cell(numAmps,1);
    for i = 1:numAmps
        tickLabels{i} = sprintf('Amp=%d pA', whichAmp(i));
    end

    % plot individual flies
    plot(xVec, allFliesMeans, ...
        'Marker', '.','LineWidth',0.5, 'Color', c(1,:));

    hold on;

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

    xlabel('IInj amplitude');
    ylabel(whichParam);
    title(sprintf('%s difference t=%.1f - t=%.1f', whichParam, ...
        plotTimes(2), plotTimes(1)));
end