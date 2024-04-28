% plotBallLegStepParamIndpt_diffAEPPEP_allFlies.m
%
% Function that takes output from saveBallLegStepParamCond_indpt() and
%  plots the difference in AEP or PEP for each leg in the XY plane
%  between 2 conditions.
% Different conditions are different output files. Matched by name. One
%  output file per fly. Has option to plot mean across flies (X) and/or
%  mean for each individual fly (dot).
% Select output files through GUI. Select set 1 then set 2. Data plotted as
%  set2 - set1
%
% INPUTS:
%   datDir - directory with output files
%   whichParam - which step parameter to plot ('AEP' or 'PEP')
%   whichPhase - which step phase ('stance' or 'swing')
%   xScale - x scale for plots. If plotting difference, as [min max]. If
%       plotting actual values, as [min max; min max; min max] for front,
%       mid, hind leg pairs
%   yScale - scale for plots. If plotting difference, as [min max]. If
%       plotting actual values, as [min max; min max; min max] for front,
%       mid, hind leg pairs
%   plotIndiv - boolean for whether to plot individual flies 
%
% OUTPUTS: none, but generates figure
%
% CREATED: 4/26/24 - HHY
%
% UPDATED:
%   4/27/24 - HHY
%
function plotBallLegStepParamIndpt_diffAEPPEP_allFlies(datDir, ...
    whichParam, whichPhase, xScale, yScale, plotIndiv)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % prompt user to select output files from 
    %  saveBallLegStepParamCond_indpt()
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
        
        % files have same names, different paths
        outputFullPath1 = [outputPath1 outName];
        outputFullPath2 = [outputPath2 outName];

        if (strcmpi(whichPhase, 'stance'))
            % load data, both sets, stance
            load(outputFullPath1, 'selStanceParams', 'legIDs');
            set1Params = selStanceParams;
            load(outputFullPath2, 'selStanceParams', 'legIDs');
            set2Params = selStanceParams;
        elseif (strcmpi(whichPhase, 'swing'))
            % load data, both sets, swing
            load(outputFullPath1, 'selSwingParams', 'legIDs');
            set1Params = selSwingParams;
            load(outputFullPath2, 'selSwingParams', 'legIDs');
            set2Params = selSwingParams;
        else
            disp('Invalid step phase. Swing or stance only.')
            return;
        end

        % param values for AEP/PEP specifically
        if (strcmpi(whichParam, 'AEP'))
            thisParamXSet2 = set2Params.stepAEPX;
            thisParamYSet2 = set2Params.stepAEPY;

            thisParamXSet1 = set1Params.stepAEPX;
            thisParamYSet1 = set1Params.stepAEPY;
        elseif (strcmpi(whichParam, 'PEP'))
            thisParamXSet2 = set2Params.stepPEPX;
            thisParamYSet2 = set2Params.stepPEPY;

            thisParamXSet1 = set1Params.stepPEPX;
            thisParamYSet1 = set1Params.stepPEPY;
        else
            disp('Invalid step parameter. AEP or PEP only.')
            return;
        end

        % compute means for this fly, for each leg
        % preallocate column vector for mean for each leg
        thisMeanX = zeros(length(legIDs.ind),1);
        thisMeanY = zeros(length(legIDs.ind),1);

        % loop through all legs
        for j = 1:length(legIDs.ind)
            thisLegInd = legIDs.ind(j);

            % set 1 param values for this leg
            thisLegXParam1 = ...
                thisParamXSet1(set1Params.stepWhichLeg == thisLegInd);
            thisLegYParam1 = ...
                thisParamYSet1(set1Params.stepWhichLeg == thisLegInd);

            % set 2 param values for this leg
            thisLegXParam2 = ...
                thisParamXSet2(set2Params.stepWhichLeg == thisLegInd);
            thisLegYParam2 = ...
                thisParamYSet2(set2Params.stepWhichLeg == thisLegInd);

            % get difference of means, also eliminate any NaNs
            thisMeanX(j) = mean(thisLegXParam2(~isnan(thisLegXParam2))) - ...
                mean(thisLegXParam1(~isnan(thisLegXParam1)));
            thisMeanY(j) = mean(thisLegYParam2(~isnan(thisLegYParam2))) - ...
                mean(thisLegYParam1(~isnan(thisLegYParam1)));
        end

        % save means for this fly in matrix for all flies
        allFliesXMeans = cat(2, allFliesXMeans, thisMeanX);
        allFliesYMeans = cat(2, allFliesYMeans, thisMeanY);
    end

    % compute mean, SEM across all flies
    meanXAllFlies = mean(allFliesXMeans, 2);
    meanYAllFlies = mean(allFliesYMeans, 2);
    SEMXAllFlies = std(allFliesXMeans, [], 2) / sqrt(numFlies);
    SEMYAllFlies = std(allFliesYMeans, [], 2) / sqrt(numFlies);


    % initialize figure
    figure;
    c = colormap('lines');

    % loop across all legs, one subplot per leg
    for j = 1:6
        % which subplot
        subplot(3,2,subInd(j));
        hold on;

        % plot indiv
        if (plotIndiv)
            % plot mean
            plot(meanYAllFlies(j), ...
                meanXAllFlies(j),...
                'Marker','x','LineStyle', 'none', 'Color', c(1,:));
            % plot indiv
            plot(allFliesYMeans(j,:), ...
                allFliesXMeans(j,:), 'Marker', '.', ...
                'LineStyle','none', 'Color', c(1,:));
        % no indiv, plot mean +/- SEM instead
        else
            errorbar(meanYAllFlies(j),...
                meanXAllFlies(j), SEMXAllFlies(j), ...
                SEMXAllFlies(j), SEMYAllFlies(j), ...
                SEMYAllFlies(j), 'Marker', 'x', ...
                'LineStyle', 'none', 'CapSize', 0, 'Color', c(1,:));
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

    end

    % figure title
    sgtitle(['Difference in ' whichParam ' ' whichPhase]);
end