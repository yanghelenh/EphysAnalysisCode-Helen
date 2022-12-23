% plotNormAEPPEP_allFlies.m
%
% Function that takes output files from saveAEPPEPbyCond_fly() and
%  generates a plot of the AEP or PEP for each leg, normalized to the
%  unmanipulated condition. Plots mean +/- SEM for each fly. Subplots for 
%  each leg (left legs on left, front legs on top)
% Select output files through GUI
%
% INPUTS:
%   datDir - directory with output files
%   xyScale - scale for plots, as [min max]
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 12/2/22 - HHY
%
% UPDATED:
%   12/2/22 - HHY
%   12/6/22 - HHY - 
%
function plotNormAEPPEP_allFlies(datDir, xyScale)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % prompt user to select output files from saveAEPPEPbyCond_fly()
    [outputFNames, outputPath] = uigetfile('*.mat', 'Select AEP PEP files', ...
        datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(outputFNames))
        numFlies = length(outputFNames);
    else
        numFlies = 1;
    end

    % initialize figures
    aepFig = figure;
    pepFig = figure;
    c = colormap('lines');

    for i = 1:numFlies

        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath outName];

        % load data
        load(outputFullPath, 'AEPxMeans', 'AEPxSEM', 'AEPyMeans',...
            'AEPySEM', 'PEPxMeans', 'PEPxSEM', 'PEPyMeans',...
            'PEPySEM', 'flipLegsLR');

        % normalize means for this fly
        normAEPxMeans = AEPxMeans(2:end,:) - repmat(AEPxMeans(1,:),...
            length(AEPxMeans(2:end,1)),1);
        normAEPyMeans = AEPyMeans(2:end,:) - repmat(AEPyMeans(1,:),...
            length(AEPyMeans(2:end,1)),1);

        normPEPxMeans = PEPxMeans(2:end,:) - repmat(PEPxMeans(1,:),...
            length(PEPxMeans(2:end,1)),1);
        normPEPyMeans = PEPyMeans(2:end,:) - repmat(PEPyMeans(1,:),...
            length(PEPyMeans(2:end,1)),1);

        % if flipped legs LR, also need to invert y axis for difference
        if (flipLegsLR)
            normAEPyMeans = -1 * normAEPyMeans;
            normPEPyMeans = -1 * normPEPyMeans;
        end

        numCond = length(AEPxMeans(2:end,1));

        % plot this fly's data into AEP figure
        figure(aepFig);

        for j = 1:6
            subplot(3,2,subInd(j));
            hold on;

            for k = 1:numCond
                errorbar(normAEPyMeans(k,j), normAEPxMeans(k,j), ...
                    AEPxSEM(k+1,j), AEPxSEM(k+1,j), ...
                    AEPySEM(k+1,j), AEPySEM(k+1,j), ...
                    'Marker','x', 'LineStyle','none','Color',c(k,:), ...
                    'LineWidth',1, 'CapSize', 0);
            end

            xlim(xyScale);
            ylim(xyScale);

            % reverse y axis (x values) so head (neg vals) is at top
            set(gca, 'YDir','reverse');

            line([0 0],ylim, 'Color','k', 'LineWidth', 2);
            line(xlim,[0,0],'Color','k', 'LineWidth', 2);

            pbaspect([1 1 1]);
        end

        % plot this fly's data into PEP figure
        figure(pepFig);

        for j = 1:6
            subplot(3,2,subInd(j));
            hold on;

            for k = 1:numCond
                errorbar(normPEPyMeans(k,j), normPEPxMeans(k,j), ...
                    PEPxSEM(k+1,j), PEPxSEM(k+1,j), ...
                    PEPySEM(k+1,j), PEPySEM(k+1,j), ...
                    'Marker','x', 'LineStyle','none','Color',c(k,:), ...
                    'LineWidth',1, 'CapSize', 0);
            end

            xlim(xyScale);
            ylim(xyScale);

            % reverse y axis (x values) so head (neg vals) is at top
            set(gca, 'YDir','reverse');

            line([0 0],ylim, 'Color','k', 'LineWidth', 2);
            line(xlim,[0,0],'Color','k', 'LineWidth', 2);

            pbaspect([1 1 1]);
        end
    end

    figure(aepFig);
    sgtitle('AEP');

    figure(pepFig);
    sgtitle('PEP');


end
