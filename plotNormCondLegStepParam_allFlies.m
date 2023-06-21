% plotNormCondLegStepParam_allFlies.m
%
% Function that takes output files from saveLegStepParamByCond_fly() and
%  generates a plot of the legStep param for each leg. Plots mean +/- SEM 
%  for each fly. Different conditions in different columns. Points from the
%  same fly are connected by lines. 
% Select output files through GUI
%
% INPUTS:
%   datDir - directory with output files
%   yScale - scale for plots, as [min max]
%   plotAvg - 'mean', 'median', or 'neither' for whether to plot mean,
%       median, or neither across all flies
%   plotIndiv - boolean for whether to plot individual flies
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 2/14/23 - HHY
%
% UPDATED:
%   2/14/23 - HHY
%   5/18/23 - HHY - add plotting of average across flies
%   6/14/23 - HHY - deal with circular variables appropriately
%
function plotNormCondLegStepParam_allFlies(datDir, yScale, plotAvg, ...
    plotIndiv)

    % all the step parameters where values need to be * -1 for left turns
    flipStepParams = {'stepVelY', 'stepAEPY', 'stepPEPY', 'stepFtLat',...
        'stepFtYaw', 'stepDirections'};
    % all step parameters that are circular variables - need to use
    %  circular stats - 6/14/23 - HHY
    circStepParams = {'stepDirections'};

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

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
    allFliesNormMeans = [];
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
        load(outputFullPath, 'stepValMeans', 'stepValSEM', 'legStepParam',...
            'stim', 'flipLegsLR', 'whichPhase');

        % normalize means for this fly
        normStepValMeans = stepValMeans(2:end,:) - repmat(stepValMeans(1,:),...
            length(stepValMeans(2:end,1)),1);
        % if circular parameter, wrap again to +/-180 deg
        if(any(strcmpi(legStepParam, circStepParams)))
            normStepValMeans = wrapTo180(normStepValMeans);
        end

        % if flipped LR and step direction, need to invert
        if (flipLegsLR)
            if (any(strcmpi(legStepParam, flipStepParams)))
                normStepValMeans = -1 * normStepValMeans;
            end
        end

        % number of conditions
        numCats = size(stepValMeans,1) - 1;

        % save means and SEMs for this fly
        allFliesNormMeans = cat(3, allFliesNormMeans, normStepValMeans);
        allFliesSEMs = cat(3, allFliesSEMs, stepValSEM(2:end,:));
        
        % x vector for plotting (for number of conditions)
        xVec = 1:numCats;
    end

    % compute mean, median, SEM across all flies
    meanAllFlies = zeros(size(allFliesNormMeans,1),size(allFliesNormMeans, 2));
    medianAllFlies = zeros(size(allFliesNormMeans,1),size(allFliesNormMeans, 2));
    SEMAllFlies = zeros(size(allFliesNormMeans,1),size(allFliesNormMeans, 2));
    % if circular
    if(any(strcmpi(legStepParam, circStepParams)))
        for i = 1:size(allFliesNormMeans,1)
            for j = 1:size(allFliesNormMeans, 2)
                thisRow = squeeze(allFliesNormMeans(i,j,:));
                thisMean = rad2deg(circ_mean(deg2rad(thisRow(~isnan(thisRow)))));
                thisMedian = rad2deg(circ_median(deg2rad(thisRow(~isnan(thisRow)))));
                thisStd = rad2deg(circ_std(deg2rad(thisRow(~isnan(thisRow)))));
                thisN = length(thisRow(~isnan(thisRow)));

                meanAllFlies(i,j) = thisMean;
                medianAllFlies(i,j) = thisMedian;
                SEMAllFlies(i,j) = thisStd / sqrt(thisN);

%         meanAllFlies = rad2deg(circ_mean(deg2rad(allFliesNormMeans), [], 3));
%         %medianAllFlies = rad2deg(circ_median(deg2rad(allFliesNormMeans), 3));
%         SEMAllFlies = rad2deg(circ_std(deg2rad(allFliesNormMeans),[],[],3)) / ...
%             sqrt(size(allFliesNormMeans,3));
            end
        end
    % not circular
    else
        for i = 1:size(allFliesNormMeans,1)
            for j = 1:size(allFliesNormMeans, 2)
                thisRow = allFliesNormMeans(i,j,:);
                thisMean = mean(thisRow(~isnan(thisRow)));
                thisMedian = median(thisRow(~isnan(thisRow)));
                thisStd = std(thisRow(~isnan(thisRow)));
                thisN = length(thisRow(~isnan(thisRow)));

                meanAllFlies(i,j) = thisMean;
                medianAllFlies(i,j) = thisMedian;
                SEMAllFlies(i,j) = thisStd / sqrt(thisN);
            end
        end
%         meanAllFlies = mean(allFliesNormMeans, 3);
%         medianAllFlies = median(allFliesNormMeans, 3);
%         SEMAllFlies = std(allFliesNormMeans,[],3) / sqrt(size(allFliesNormMeans,3));
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

        % if plot mean across flies
        if (strcmpi(plotAvg,'mean') && plotIndiv)
%             errorbar(xVec', meanAllFlies(:,j), SEMAllFlies(:,j), ...
%                 '_','LineWidth', 2, 'CapSize', 0, 'Color', c(2,:));
            plot(xVec', meanAllFlies(:,j), ...
                'Marker','_','LineWidth', 2, 'Color', c(2,:));
        elseif (strcmpi(plotAvg,'mean') && ~plotIndiv)
            errorbar(xVec', meanAllFlies(:,j), SEMAllFlies(:,j), ...
                '_','LineWidth', 2, 'CapSize', 0, 'Color', c(2,:));
        elseif (strcmpi(plotAvg,'median') && plotIndiv)
            plot(xVec', medianAllFlies(:,j), ...
                '_','LineWidth', 2, 'Color', c(2,:));
        elseif (strcmpi(plotAvg,'median') && ~plotIndiv)
            errorbar(xVec', meanAllFlies(:,j), SEMAllFlies(:,j), ...
                '_','LineWidth', 2, 'CapSize', 0, 'Color', c(2,:));
        end

        % if plot indiv flies
        if (plotIndiv)
            for i = 1:size(allFliesNormMeans,3)
                if (strcmpi(plotAvg, 'none'))
                    errorbar(xVec', allFliesNormMeans(:,j,i), allFliesSEMs(:,j,i), ...
                        'Marker', '.','LineWidth',1, 'CapSize', 0, 'Color', c(1,:));
                else
                    plot(xVec', allFliesNormMeans(:,j,i), ...
                        'Marker', '.','LineWidth',0.5, 'Color', c(1,:));
                end
            end
        end

        ylim(yScale);
        xScale = xlim;
        xScale(1) = xScale(1) - (0.1 * (xVec(end)-xVec(1)));
        xScale(2) = xScale(2) + (0.1 * (xVec(end)-xVec(1)));
        xlim(xScale);


        line(xScale,[0,0],'Color','k', 'LineWidth', 2);

        xticks(xVec);
        xticklabels(repmat({''},1,numCats));
    end

    

    % generate labels for x axis
    xLabels = cell(size(xVec));

    if (strcmpi(stim.whichStim, 'iInj'))
        % generate key for mapping b/w indices and amps and durs
        iInjCatAmps = zeros(1,numCats);
        iInjCatDurs = zeros(1,numCats);
        % counter index into vectors
        counter = 1; 
    
        % assign amps to indices
        for i = 1:length(stim.amps)
            for j = 1:length(stim.durs)
                iInjCatAmps(counter) = stim.amps(i);
                iInjCatDurs(counter) = stim.durs(j);
    
                counter = counter + 1;
            end
        end
    elseif (strcmpi(stim.whichStim, 'opto'))
        % generate key for mapping b/w indices and NDs and durs
        optoCatNDs = ones(1,numCats) * -1;
        optoCatDurs = ones(1,numCats) * -1;
        % counter index into vectors
        counter = 1; 

        % assign NDs, durs to indices
        for i = 1:length(stim.NDs)
            for j = 1:length(stim.durs)
                optoCatNDs(counter) = stim.NDs(i);
                optoCatDurs(counter) = stim.durs(j);

                counter = counter + 1;
            end
        end
    end

    for i = 1:numCats
        if (strcmpi(stim.whichStim, 'iInj'))
            xLabels{i} = sprintf('%d pA, %.1f s', iInjCatAmps(i), ...
                iInjCatDurs(i));
        elseif (strcmpi(stim.whichStim, 'opto'))
            xLabels{i} = sprintf('ND=%.1f, %.1f s', optoCatNDs(i), ...
                optoCatDurs(i));
        end
    end

%     xticklabels(xLabels);


%     sgtitle(legStepParam);
end
