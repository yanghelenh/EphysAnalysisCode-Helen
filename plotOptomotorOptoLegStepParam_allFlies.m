% plotOptomotorOptoLegStepParam_allFlies.m
%
% Function that takes output files from 
%  extractOptomotorLegStepParamsOptoCond_fly() and generates a plot of the 
%  specified legStep param for each leg. Plots mean +/- SEM for each fly as 
%  well as across flies. Different conditions in different columns. Points 
%  from the same fly are connected by lines. 
% Slight modification of plotOptoLegStepParam_allFlies(), has vels instead
%  of durs
% Select output files through GUI
%
% INPUTS:
%   datDir - directory with output files
%   whichParam - which step parameter to plot
%   vels - which velocity conditions to plot
%   NDs - which NDs to plot
%   yScale - scale for plots, as [min max]
%   plotAvg - 'mean' or 'median' for which type of average across flies
%   plotIndiv - boolean for whether to plot individual flies
%   plotDiff - boolean for whether to plot absolute value of step parameter
%       or difference from no stim condition
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 1/29/24 - HHY
%
% UPDATED:
%   1/29/24 - HHY
%
function plotOptomotorOptoLegStepParam_allFlies(...
    datDir, whichParam, whichPhase,vels, NDs, yScale, plotAvg, ...
    plotIndiv, plotDiff)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % circular step parameters
    circStepParams = {'stepDirections'};

    % prompt user to select output files from saveLegStepParamByCond_fly()
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
        load(outputFullPath, 'legStepsOptoMeans', 'legStepsOptoSEM', ...
            'condKeyVels', 'condKeyNDs');

        if (strcmpi(whichPhase, 'stance'))
%             if(any(strcmpi(whichParam, circStepParams)))
%                 thisMean = wrapTo360(legStepsOptoMeans.stance.(whichParam));
%             else
                thisMean = legStepsOptoMeans.stance.(whichParam);
%             end
            thisSEM = legStepsOptoSEM.stance.(whichParam);
        elseif (strcmpi(whichPhase, 'swing'))
            thisMean = legStepsOptoMeans.swing.(whichParam);
            thisSEM = legStepsOptoSEM.swing.(whichParam);
        end
        
        % if we're plotting difference from no stim condition
        if (plotDiff)
            noStimInd = find(condKeyNDs == -1);
            for j = 1:length(noStimInd)
                thisVel = condKeyVels(noStimInd(j));
                sameVelInd = find(condKeyVels == thisVel);

                noStimMean = thisMean(noStimInd(j),:);
                for k = 1:length(sameVelInd)
%                     thisMean(sameDurInd(k),:) = thisMean(sameDurInd(k),:) - ...
%                         noStimMean;

                    if(any(strcmpi(whichParam, circStepParams)))
                        thisMean(sameVelInd(k),:) = wrapTo180(...
                            thisMean(sameVelInd(k),:) - ...
                            noStimMean); 
                    else
                        thisMean(sameVelInd(k),:) = thisMean(sameVelInd(k),:) - ...
                            noStimMean;
                    end
                end
            end
        end


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

    % get which conditions to plot (of specified durs and NDs)
    plotCondInd = [];
    for i = 1:length(condKeyVels)
        if (any(condKeyVels(i) == vels) && any(condKeyNDs(i) == NDs))
            plotCondInd = [plotCondInd; i];
        end
    end

    xVec = 1:length(plotCondInd);

    % get labels
    tickLabels = cell(size(plotCondInd));
    for i = 1:length(plotCondInd)
        tickLabels{i} = sprintf('ND=%.1f', ...
            condKeyNDs(plotCondInd(i)));
    end


    % initialize figure
    figure;
    c = colormap('lines');

    % loop across all legs, one subplot per leg
    for j = 1:6
        subplot(3,2,subInd(j));
        hold on;

        % if plot mean across flies
        if (strcmpi(plotAvg,'mean') && plotIndiv)
%             errorbar(xVec', meanAllFlies(:,j), SEMAllFlies(:,j), ...
%                 '_','LineWidth', 2, 'CapSize', 0, 'Color', c(2,:));
            plot(xVec', meanAllFlies(plotCondInd,j), ...
                'Marker','_','LineWidth', 2, 'Color', c(2,:));
        elseif (strcmpi(plotAvg,'mean') && ~plotIndiv)
            errorbar(xVec', meanAllFlies(plotCondInd,j), SEMAllFlies(plotCondInd,j), ...
                '_','LineWidth', 2, 'CapSize', 0, 'Color', c(2,:));
        elseif (strcmpi(plotAvg,'median') && plotIndiv)
            plot(xVec', medianAllFlies(plotCondInd,j), ...
                '_','LineWidth', 2, 'Color', c(2,:));
        elseif (strcmpi(plotAvg,'median') && ~plotIndiv)
            errorbar(xVec', meanAllFlies(plotCondInd,j), SEMAllFlies(plotCondInd,j), ...
                '_','LineWidth', 2, 'CapSize', 0, 'Color', c(2,:));
        end

        % if plot indiv flies
        if (plotIndiv)
            for i = 1:size(allFliesMeans,3)
                if (strcmpi(plotAvg, 'none'))
                    errorbar(xVec', allFliesMeans(plotCondInd,j,i), allFliesSEMs(plotCondInd,j,i), ...
                        'Marker', '.','LineWidth',1, 'CapSize', 0, 'Color', c(1,:));
                else
                    plot(xVec', allFliesMeans(plotCondInd,j,i), ...
                        'Marker', '.','LineWidth',0.5, 'Color', c(1,:));
                end
            end
        end

        ylim(yScale);
        xScale = xlim;
        xScale(1) = xScale(1) - (0.1 * (xVec(end)-xVec(1)));
        xScale(2) = xScale(2) + (0.1 * (xVec(end)-xVec(1)));
        xlim(xScale);

        if (plotDiff)
            line(xScale,[0,0],'Color','k', 'LineWidth', 2);
        end

        xticks(xVec);
        if (j == 3 || j == 6)
            xticklabels(tickLabels);
        else
            xticklabels(repmat({''},1,length(plotCondInd)));
        end
    end

    
    sgtitle(sprintf('%s, %s, Vel=%.1f',whichParam, whichPhase, vels));
end
