% plotOpto_AEPPEP_allFlies.m
%
% Function that takes output from extractLegStepParamsOpto_fly() and
%  plots the AEP/PEP in the XY plane, either as the actual values or as the
%  difference from the no stimulation condition
% Has option to plot mean across flies (X) and/or mean for each individual 
%  fly (dot).
% Select output files through GUI. One output file per fly
% Modification of: plotNormAEPPEP_allFlies()
%
% INPUTS:
%   datDir - directory with output files
%   whichParam - which step parameter to plot ('AEP' or 'PEP')
%   whichPhase - which step phase ('stance' or 'swing')
%   durs - which duration conditions to plot
%   NDs - which NDs to plot
%   xScale - x scale for plots. If plotting difference, as [min max]. If
%       plotting actual values, as [min max; min max; min max] for front,
%       mid, hind leg pairs
%   yScale - y scale for plots. If plotting difference, as [min max]. If
%       plotting actual values, as [min max; min max; min max] for front,
%       mid, hind leg pairs
%   plotIndiv - boolean for whether to plot individual flies
%   plotDiff - boolean for whether to plot absolute value of step parameter
%       or difference from no stim condition
%
% OUTPUTS: none, but generates figure
%
% CREATED: 4/26/24 - HHY
%
% UPDATED:
%   4/27/24 - HHY
%
function plotOpto_AEPPEP_allFlies(datDir, whichParam, whichPhase, ...
    durs, NDs, xScale, yScale, plotIndiv, plotDiff)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % prompt user to select output files from extractLegStepParamsOpto_fly()
    [outputFNames, outputPath] = uigetfile('*.mat', ['...' ...
        'Select Step Param files'], datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(outputFNames))
        numFlies = length(outputFNames);
    else
        numFlies = 1;
    end

    % initialize
    allFliesXMeans = [];
    allFliesYMeans = [];

    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];

        % load data
        load(outputFullPath, 'legStepsOptoMeans', 'condKeyDurs', ...
            'condKeyNDs');

        if (strcmpi(whichPhase, 'stance'))
            if (strcmpi(whichParam, 'AEP'))
                thisMeanX = legStepsOptoMeans.stance.stepAEPX;
                thisMeanY = legStepsOptoMeans.stance.stepAEPY;
            elseif (strcmpi(whichParam, 'PEP'))
                thisMeanX = legStepsOptoMeans.stance.stepPEPX;
                thisMeanY = legStepsOptoMeans.stance.stepPEPY;
            else
                disp('Invalid step parameter. AEP or PEP only.')
                return;
            end
        elseif (strcmpi(whichPhase, 'swing'))
            if (strcmpi(whichParam, 'AEP'))
                thisMeanX = legStepsOptoMeans.swing.stepAEPX;
                thisMeanY = legStepsOptoMeans.swing.stepAEPY;
            elseif (strcmpi(whichParam, 'PEP'))
                thisMeanX = legStepsOptoMeans.swing.stepPEPX;
                thisMeanY = legStepsOptoMeans.swing.stepPEPY;
            else
                disp('Invalid step parameter. AEP or PEP only.')
                return;
            end
        else
            disp('Invalid step phase. Swing or stance only.')
            return;
        end
        
        % if we're plotting difference from no stim condition
        %  find difference in X and Y means
        if (plotDiff)
            noStimInd = find(condKeyNDs == -1);
            for j = 1:length(noStimInd)
                thisDur = condKeyDurs(noStimInd(j));
                sameDurInd = find(condKeyDurs == thisDur);

                noStimMeanX = thisMeanX(noStimInd(j),:);
                noStimMeanY = thisMeanY(noStimInd(j),:);

                for k = 1:length(sameDurInd)
                    thisMeanX(sameDurInd(k),:) = ...
                        thisMeanX(sameDurInd(k),:) - noStimMeanX;
                    thisMeanY(sameDurInd(k),:) = ...
                        thisMeanY(sameDurInd(k),:) - noStimMeanY;            
                end
            end
        end


        % save means for this fly
        allFliesXMeans = cat(3, allFliesXMeans, thisMeanX);
        allFliesYMeans = cat(3, allFliesYMeans, thisMeanY);
    end

    % compute mean, SEM across all flies
    meanXAllFlies = zeros(size(allFliesXMeans,1),size(allFliesXMeans, 2));
    meanYAllFlies = zeros(size(allFliesYMeans,1),size(allFliesYMeans, 2));
    SEMXAllFlies = zeros(size(allFliesXMeans,1),size(allFliesXMeans, 2));
    SEMYAllFlies = zeros(size(allFliesYMeans,1),size(allFliesYMeans, 2));

    % X
    for i = 1:size(allFliesXMeans,1)
        for j = 1:size(allFliesXMeans, 2)
            thisRow = allFliesXMeans(i,j,:);
            thisMean = mean(thisRow(~isnan(thisRow)));
            thisStd = std(thisRow(~isnan(thisRow)));
            thisN = length(thisRow(~isnan(thisRow)));

            meanXAllFlies(i,j) = thisMean;
            SEMXAllFlies(i,j) = thisStd / sqrt(thisN);
        end
    end

    % Y
    for i = 1:size(allFliesYMeans,1)
        for j = 1:size(allFliesYMeans, 2)
            thisRow = allFliesYMeans(i,j,:);
            thisMean = mean(thisRow(~isnan(thisRow)));
            thisStd = std(thisRow(~isnan(thisRow)));
            thisN = length(thisRow(~isnan(thisRow)));

            meanYAllFlies(i,j) = thisMean;
            SEMYAllFlies(i,j) = thisStd / sqrt(thisN);
        end
    end

    % get which conditions to plot (of specified durs and NDs)
    plotCondInd = [];
    for i = 1:length(condKeyDurs)
        if (any(condKeyDurs(i) == durs) && any(condKeyNDs(i) == NDs))
            plotCondInd = [plotCondInd; i];
        end
    end

    % get condition labels for legend
    legendLabels = cell(size(plotCondInd));
    for i = 1:length(plotCondInd)
        legendLabels{i} = sprintf('ND=%.1f, dur=%.1f s', ...
            condKeyNDs(plotCondInd(i)), condKeyDurs(plotCondInd(i)));
    end

    % initialize figure
    figure;
    c = colormap('lines');

    % loop across all legs, one subplot per leg
    for j = 1:6
        % handle for mean plot
        meanHndl = [];

        % which subplot
        subplot(3,2,subInd(j));
        hold on;

        % loop through all conditions
        for k = 1:length(plotCondInd)
            thisInd = plotCondInd(k);

            % plot indiv
            if (plotIndiv)
                % plot mean
                meanHndl(k) = plot(meanYAllFlies(thisInd,j), ...
                    meanXAllFlies(thisInd,j),...
                    'Marker','x','LineStyle', 'none', 'Color', c(k,:));
                % plot indiv
                plot(squeeze(allFliesYMeans(thisInd,j,:)), ...
                    squeeze(allFliesXMeans(thisInd,j,:)), 'Marker', '.', ...
                    'LineStyle','none', 'Color', c(k,:));
            % no indiv, plot mean +/- SEM instead
            else
                meanHndl(k) = errorbar(meanYAllFlies(thisInd,j),...
                    meanXAllFlies(thisInd,j), SEMXAllFlies(thisInd,j), ...
                    SEMXAllFlies(thisInd,j), SEMYAllFlies(thisInd,j), ...
                    SEMYAllFlies(thisInd,j), 'Marker', 'x', ...
                    'LineStyle', 'none', 'CapSize', 0, 'Color', c(k,:));
            end
        end

        % scale
        if (plotDiff)

            xlim(xScale);
            ylim(yScale);

            % x and y axis lines
            line([0 0],ylim, 'Color','k', 'LineWidth', 1);
            line(xlim,[0,0],'Color','k', 'LineWidth', 1);

        else
            % scale depends on which legs
            if (j==1 || j==4) % front legs
                xlim(xScale(1,:));
                ylim(yScale(1,:));
            elseif (j==2 || j==5) % mid legs
                xlim(xScale(2,:));
                ylim(yScale(2,:)); 
            else % hind legs
                xlim(xScale(3,:));
                ylim(yScale(3,:));
            end

        end

        % reverse y axis (x values) so head (neg vals) is at top
        set(gca, 'YDir','reverse');

        % x and y on same scale
        pbaspect([1 1 1]);

        % add legend to bottom right
        if (j==3)
            legend(meanHndl,legendLabels);
        end
    end

    % figure title
    if (plotDiff)
        sgtitle(['Difference in ' whichParam ' ' whichPhase]);
    else
        sgtitle([whichParam ' ' whichPhase]);
    end
end