% plotSpikerateAutocorr.m
%
% Function to plot autocorrelation of spike rate (output of
%  getSpikerateAutocorr_cell()). Plots the autocorrelation itself as well
%  as the FWHM, across conditions. Plots one line/pt per cell and mean +/-
%  SEM across cells
% Select output files through GUI
%
% INPUTS:
%   datDir - folder containing getSpikerateAutocorr_cell() output files
%   numCond - number of conditions
%   condNames - cell array of length numCond, string for name of each 
%       condition   
%   yScaleAC - limits of y axis for plot of autocorr
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 9/1/23 - HHY
%
% UPDATED:
%   9/1/23 - HHY
%
function plotSpikerateAutocorr(datDir, numCond, condNames, yScaleAC)

    % preallocate 
    % cell array where each element will be numCells x numTPts matrix of
    %  autocorrelations
    allCellsAC = cell(numCond,1);
    % cell array where each element will be numCells length vector of FWHM
    allCellsFWHM = cell(numCond,1);

    % cell array where each element will be 1 x numTPts vector of mean/SEM
    %  across flies
    allCondMeanAC = cell(numCond,1);
    allCondSEMsAC = cell(numCond,1);
    allCondMeanFWHM = zeros(numCond,1);
    allCondSEMsFWHM = zeros(numCond,1);
    allCondLagsT = cell(numCond,1);

    % loop across number of conditions, get data files for each cell
    % compute mean across flies
    for i = 1:numCond
        [outputFNames, outputPath] = uigetfile('*.mat', ...
            'Select getSpikerateAutocorr_cell() files', ...
            datDir, 'MultiSelect', 'on');

        % if only 1 file selected, not cell array; make sure loop still
        %  works 
        % num flies is number of cells
        if (iscell(outputFNames))
            numCells = length(outputFNames);
        else
            numCells = 1;
        end

        % preallocate
        thisCondAC = [];
        thisCondFWHM = zeros(numCells,1);

        % loop through all cells
        for j = 1:numCells
            % handle whether it's a cell array or not
            if (iscell(outputFNames))
                outName = outputFNames{j};
            else
                outName = outputFNames;
            end
            
            outputFullPath = [outputPath outName];

            % load data file
            load(outputFullPath, 'autoCorr', 'lagsT');

            % get FWHM for this cell
            thisCondFWHM(j) = fwhm(lagsT, autoCorr);

            % save this autocorr
            thisCondAC = [thisCondAC, autoCorr];

        end

        % autocorr, FWHM, lagsT across cond
        allCellsAC{i} = thisCondAC;
        allCellsFWHM{i} = thisCondFWHM;
        allCondLagsT{i} = lagsT;

        % get mean and SEM
        allCondMeanAC{i} = mean(thisCondAC,2);
        allCondSEMsAC{i} = std(thisCondAC,[],2) / sqrt(numCells);
        allCondMeanFWHM(i) = mean(thisCondFWHM);
        allCondSEMsFWHM(i) = std(thisCondFWHM) / sqrt(numCells);
    end

    % plot autocorr
    figure;
    c = colormap('lines');
    hold on;

    legendInd = [];
    % loop over all conditions
    for i = 1:numCond
        % plot individual flies
        plot(allCondLagsT{i}, allCellsAC{i}, ...
            'LineWidth',0.5, 'Color', c(i,:)');

         % plot mean across flies
        legendInd(i) = plot(allCondLagsT{i}, allCondMeanAC{i}, ...
            'LineWidth',2, 'Color', c(i,:));
    end

    % label axes
    ylim(yScaleAC);
    ylabel('Autocorrelation');

    xlim([allCondLagsT{1}(1) allCondLagsT{1}(end)]);
    xlabel('Time (s)');

    legend(legendInd,condNames);


    % plot FWHM
    figure;
    hold on;

    xVec = 1:numCond;

    for i = 1:numCond
        % plot individual cells as points
        numCells = length(allCellsFWHM{i});
        % get repmat xVec for this number of cells
        repXVec = repmat(xVec(i),numCells,1);

        plot(repXVec, allCellsFWHM{i}, ...
            '.','LineWidth',0.5, 'Color', c(1,:));

        % plot mean +/- SEM across flies
        errorbar(xVec(i),allCondMeanFWHM(i),allCondSEMsFWHM(i),...
            '_','LineWidth', 2, 'CapSize', 0, 'Color', c(2,:));
    end

    xScale = xlim;
    xScale(1) = xScale(1) - (0.5 * (xVec(end)-xVec(1)));
    xScale(2) = xScale(2) + (0.5 * (xVec(end)-xVec(1)));
    xlim(xScale);

    xticks(xVec);
    xticklabels(condNames);

    ylabel('FWHM (s)');
end