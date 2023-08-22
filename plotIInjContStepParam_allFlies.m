% plotIInjContStepParam_allFlies.m
%
% Function that takes output files from extractContStepParamsIInj_fly() and
%  generates a plot of the specified legStep param for each leg. Plots 
%  mean for each fly as well as across flies. User selects one dur, one amp,
%  one step parameter
% Select output files through GUI
%
% INPUTS:
%   datDir - directory with output files
%   whichParam - which step parameter to plot
%   whichDur - which duration conditions to plot
%   whichAmp - which amps to plot
%   yScale - scale for plots, as [min max]
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 8/21/23 - HHY
%
% UPDATED:
%   8/21/23 - HHY
%
function plotIInjContStepParam_allFlies(datDir, whichParam, whichDur,...
    whichAmp, yScale)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % circular step parameters
    circStepParams = {'stepDirections'};

    % prompt user to select output files from extractContStepParamsIInj_fly()
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
        load(outputFullPath, 'legStepsContIInjMean', 'legStepsContIInjSEM', ...
            'durs', 'amps', 'trialTimes');

        thisDurInd = find(durs == whichDur);
        thisAmpInd = find(amps == whichAmp);
        thisTrialTimes = trialTimes{thisDurInd};

        thisMean = legStepsContIInjMean(thisDurInd, thisAmpInd).(whichParam);
        thisSEM = legStepsContIInjSEM(thisDurInd, thisAmpInd).(whichParam);


        % save means and SEMs for this fly
        allFliesMeans = cat(3, allFliesMeans, thisMean);
        allFliesSEMs = cat(3, allFliesSEMs, thisSEM);
    end


    % compute mean, median, SEM across all flies
    meanAllFlies = zeros(size(allFliesMeans,1),size(allFliesMeans, 2));
    medianAllFlies = zeros(size(allFliesMeans,1),size(allFliesMeans, 2));
    SEMAllFlies = zeros(size(allFliesMeans,1),size(allFliesMeans, 2));
    % if circular
    if(any(strcmpi(whichParam, circStepParams)))
        for i = 1:size(allFliesMeans,1)
            for j = 1:size(allFliesMeans, 2)
                thisRow = squeeze(allFliesMeans(i,j,:));
                thisMean = rad2deg(circ_mean(deg2rad(thisRow(~isnan(thisRow)))));
                thisMedian = rad2deg(circ_median(deg2rad(thisRow(~isnan(thisRow)))));
                thisStd = rad2deg(circ_std(deg2rad(thisRow(~isnan(thisRow)))));
                thisN = length(thisRow(~isnan(thisRow)));

                meanAllFlies(i,j) = thisMean;
                SEMAllFlies(i,j) = thisStd / sqrt(thisN);
            end
        end
    % not circular
    else
        for i = 1:size(allFliesMeans,1)
            for j = 1:size(allFliesMeans, 2)
                thisRow = allFliesMeans(i,j,:);
                thisMean = mean(thisRow(~isnan(thisRow)));
                thisMedian = median(thisRow(~isnan(thisRow)));
                thisStd = std(thisRow(~isnan(thisRow)));
                thisN = length(thisRow(~isnan(thisRow)));

                meanAllFlies(i,j) = thisMean;
                medianAllFlies(i,j) = thisMedian;
                SEMAllFlies(i,j) = thisStd / sqrt(thisN);
            end
        end

    end

%     % for plotting, if circular parameter and stance, wrap to 360
%     if (strcmpi(whichPhase, 'stance') && any(strcmpi(legStepParam, circStepParams)))
%         allFliesNormMeans = wrapTo360(allFliesNormMeans);
%         meanAllFlies = wrapTo360(meanAllFlies);
%         %medianAllFlies = wrapTo360(medianAllFlies);
%     end

    % initialize figure
    figure;
    c = colormap('lines');

    % loop across all legs, one subplot per leg
    for j = 1:6
        subplot(3,2,subInd(j));
        hold on;

        % plot indiv flies
        plot(thisTrialTimes, squeeze(allFliesMeans(:,j,:)),...
            'LineWidth',0.5, 'Color', c(1,:));

        % plot mean
        plot(thisTrialTimes, meanAllFlies(:,j), 'Color', c(2,:), ...
            'LineWidth', 2);


        ylim(yScale);
        xScale = [thisTrialTimes(1) thisTrialTimes(end)];
        xlim(xScale)

        % plot lines for opto
        line([0,0], yScale, 'Color','k', 'LineWidth', 1);
        line([whichDur whichDur], yScale, 'Color','k', 'LineWidth', 1);
    end

    titleStr = sprintf('%s, Amp=%d pA, dur=%.1f s', whichParam, whichAmp, ...
        whichDur);
    sgtitle(titleStr);
end
