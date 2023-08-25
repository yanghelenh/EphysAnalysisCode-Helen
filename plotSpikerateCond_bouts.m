% plotSpikerateCond_bouts.m
%
% Function to plot output of saveSpikerateCond_bouts(), pooled across
%  multiple flies and conditions. 
% User specifies number of conditions, and then selects files (1 per fly) 
%  for each condition through GUI. Plots mean for each fly and mean across 
%  flies for each condition (different colors for different conditions). 
%  Thin lines for individual flies; thick lines for means across flies.
%
% INPUTS:
%   datDir - full path to directory with output files of
%       saveSpikerateCond_bouts()
%   numCond - number of conditions to plot (must be matched in num time
%       points)
%   condNames - cell array of length numCond, for names of each condition
%   yScale - 2 element vector for y-axis limits
%
% OUTPUTS:
%   none, but produces plot
%
% CREATED: 8/24/23 - HHY
%
% UPDATED:
%   8/24/23 - HHY
%
function plotSpikerateCond_bouts(datDir, numCond, condNames)
    
    % preallocate 
    % cell array where each element will be numFlies x numTPts matrix of
    %  mean/SEM for each fly
    allFlyMeans = cell(numCond,1);
    allFlySEMs = cell(numCond,1);

    % cell array where each element will be 1 x numTPts vector of mean/SEM
    %  across flies
    allCondMeans = cell(numCond,1);
    allCondSEMs = cell(numCond,1);

    % cell arrays where each element is cell array of output file
    %  names/paths
    fileNames = cell(numCond,1);
    filePaths = cell(numCond,1);

    % loop across number of conditions, get data files for each fly
    % compute mean for each fly and mean across flies
    for i = 1:numCond
        [outputFNames, outputPath] = uigetfile('*.mat', ...
            'Select saveSpikerateCond_bouts files', ...
            datDir, 'MultiSelect', 'on');

        % if only 1 file selected, not cell array; make sure loop still
        %  works 
        % num flies is number of files
        if (iscell(outputFNames))
            numFlies = length(outputFNames);
        else
            numFlies = 1;
        end

        % preallocate
        thisCondMean = [];
        thisCondSEM = [];

        % loop through all flies
        for j = 1:numFlies
            % handle whether it's a cell array or not
            if (iscell(outputFNames))
                outName = outputFNames{j};
            else
                outName = outputFNames;
            end
            
            outputFullPath = [outputPath outName];

            % load data file
            load(outputFullPath, 'allSpikerate', 't', 'numBouts');

            % get mean subtracted
            subMean = mean(mean(allSpikerate(10:40,:),2));

            % compute mean for this fly
            thisMean = mean(allSpikerate, 2) - subMean; 
            thisSEM = std(allSpikerate, [], 2) / sqrt(numBouts);

            % save this mean and SEM
            thisCondMean = [thisCondMean, thisMean];
            thisCondSEM = [thisCondSEM, thisSEM];
        end

        % get mean and SEM across flies for this condition
        allCondMeans{i} = mean(thisCondMean,2);
        allCondSEMs{i} = std(thisCondMean, [], 2) / ...
            sqrt(size(thisCondMean,2));

        % save mean and SEM for each fly
        allFlyMeans{i} = thisCondMean;
        allFlySEMs{i} = thisCondSEM;
    end

    % plotting
    figure;
    c = colormap('lines');
    hold on;

    legendInd = [];
    % loop through conditions
    for i = 1:numCond
        % plot individual flies
        plot(t, allFlyMeans{i}, 'LineWidth',0.5, 'Color', c(i,:)');

        % plot mean across flies
        legendInd(i) = plot(t, allCondMeans{i}, ...
            'LineWidth',2, 'Color', c(i,:));
    end

    xlabel('Time relative to yaw peak (s)');
    ylabel('Spike rate (Hz)');
    legend(legendInd,condNames);
end