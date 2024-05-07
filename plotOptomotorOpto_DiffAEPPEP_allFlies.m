% plotOptomotorOpto_DiffAEPPEP_allFlies.m
%
% Function that takes output from 
%  extractOptomotorLegStepParamsOptoCond_fly() and plots the AEP/PEP in the
%  XY plane, as the difference between steps from 2 different time periods
%  (different output files, matched by name).
% Has option to plot mean across flies (X) and/or mean for each individual 
%  fly (dot).
% Select output files through GUI. One output file per fly
% Modification of: plotOpto_AEPPEP_allFlies()
%
% INPUTS:
%   datDir - directory with output files
%   whichParam - which step parameter to plot ('AEP' or 'PEP')
%   whichPhase - which step phase ('stance' or 'swing')
%   vels - which optomotor velocity conditions to plot
%   NDs - which NDs to plot
%   xScale - x scale for plots, as [min max]
%   yScale - y scale for plots, as [min max]
%   plotIndiv - boolean for whether to plot individual flies
%
% OUTPUTS: none, but generates figure
%
% CREATED: 4/30/24 - HHY
%
% UPDATED:
%   4/30/24 - HHY
%
function plotOptomotorOpto_DiffAEPPEP_allFlies(datDir, whichParam, ...
    whichPhase, vels, NDs, xScale, yScale, plotIndiv)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % prompt user to select output files from saveLegStepParamByCond_fly()
    [outputFNames1, outputPath1] = uigetfile('*.mat', ...
        'Select set 1 Step Param files', ...
        datDir, 'MultiSelect', 'on');

    [outputFNames2, outputPath2] = uigetfile('*.mat', ...
        'Select set 2 Step Param files', ...
        datDir, 'MultiSelect', 'on');


    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(outputFNames1))
        numFlies = length(outputFNames1);
    else
        numFlies = 1;
    end

    % initialize
    allFliesXMeans = [];
    allFliesYMeans = [];


    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(outputFNames1))
            outName = outputFNames1{i};
        else
            outName = outputFNames1;
        end
        
        outputFullPath1 = [outputPath1 outName];

        % load data, set 1
        load(outputFullPath1, 'legStepsOptoMeans', 'legStepsOptoSEM', ...
            'condKeyVels', 'condKeyNDs');

        if (strcmpi(whichPhase, 'stance'))
            if (strcmpi(whichParam, 'AEP'))
                thisMeanX1 = legStepsOptoMeans.stance.stepAEPX;
                thisMeanY1 = legStepsOptoMeans.stance.stepAEPY;
            elseif (strcmpi(whichParam, 'PEP'))
                thisMeanX1 = legStepsOptoMeans.stance.stepPEPX;
                thisMeanY1 = legStepsOptoMeans.stance.stepPEPY;
            else
                disp('Invalid step parameter. AEP or PEP only.')
                return;
            end
        elseif (strcmpi(whichPhase, 'swing'))
            if (strcmpi(whichParam, 'AEP'))
                thisMeanX1 = legStepsOptoMeans.swing.stepAEPX;
                thisMeanY1 = legStepsOptoMeans.swing.stepAEPY;
            elseif (strcmpi(whichParam, 'PEP'))
                thisMeanX1 = legStepsOptoMeans.swing.stepPEPX;
                thisMeanY1 = legStepsOptoMeans.swing.stepPEPY;
            else
                disp('Invalid step parameter. AEP or PEP only.')
                return;
            end
        else
            disp('Invalid step phase. Swing or stance only.')
            return;
        end

        outputFullPath2 = [outputPath2 outName];
        
        % load data, set 2
        load(outputFullPath2, 'legStepsOptoMeans', 'legStepsOptoSEM', ...
            'condKeyVels', 'condKeyNDs');

        if (strcmpi(whichPhase, 'stance'))
            if (strcmpi(whichParam, 'AEP'))
                thisMeanX2 = legStepsOptoMeans.stance.stepAEPX;
                thisMeanY2 = legStepsOptoMeans.stance.stepAEPY;
            elseif (strcmpi(whichParam, 'PEP'))
                thisMeanX2 = legStepsOptoMeans.stance.stepPEPX;
                thisMeanY2 = legStepsOptoMeans.stance.stepPEPY;
            else
                disp('Invalid step parameter. AEP or PEP only.')
                return;
            end
        elseif (strcmpi(whichPhase, 'swing'))
            if (strcmpi(whichParam, 'AEP'))
                thisMeanX2 = legStepsOptoMeans.swing.stepAEPX;
                thisMeanY2 = legStepsOptoMeans.swing.stepAEPY;
            elseif (strcmpi(whichParam, 'PEP'))
                thisMeanX2 = legStepsOptoMeans.swing.stepPEPX;
                thisMeanY2 = legStepsOptoMeans.swing.stepPEPY;
            else
                disp('Invalid step parameter. AEP or PEP only.')
                return;
            end
        else
            disp('Invalid step phase. Swing or stance only.')
            return;
        end

        % take difference
        thisDiffMeanX = thisMeanX2 - thisMeanX1;
        thisDiffMeanY = thisMeanY2 - thisMeanY1;


        % save means for this fly
        allFliesXMeans = cat(3, allFliesXMeans, thisDiffMeanX);
        allFliesYMeans = cat(3, allFliesYMeans, thisDiffMeanY);
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
    for i = 1:length(condKeyVels)
        if (any(condKeyVels(i) == vels) && any(condKeyNDs(i) == NDs))
            plotCondInd = [plotCondInd; i];
        end
    end

    % get condition labels for legend
    legendLabels = cell(size(plotCondInd));
    for i = 1:length(plotCondInd)
        legendLabels{i} = sprintf('ND=%.1f, vel=%d', ...
            condKeyNDs(plotCondInd(i)), condKeyVels(plotCondInd(i)));
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
        xlim(xScale);
        ylim(yScale);

        % x and y axis lines
        line([0 0],ylim, 'Color','k', 'LineWidth', 1);
        line(xlim,[0,0],'Color','k', 'LineWidth', 1);

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
    sgtitle(['Difference in ' whichParam ' ' whichPhase]);
end