% plotSpikesPhaseDelays.m
%
% Function to plot amplitude of modulation of spike timing by phase, across
%  temporal offsets
% Fits 1 period sine curve to values for each delay, plots amplitude (2*a)
%  equation: a*sin(x+c) + d
% Data from getSpikesStepPhase_cell()
% 
% INPUTS:
%   datDir - directory with output files
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 9/18/23 - HHY
%
% UPDATED:
%   9/18/23 - HHY
%
function plotSpikesComboPhaseDelays(datDir)

     % prompt user to select output files from saveLegStepParamByCond_fly()
    [outputFNames, outputPath] = uigetfile('*.mat', ...
        'Select getSpikesStepPhase_cell() files', ...
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
    allAmps = [];

    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath filesep outName];

        load(outputFullPath, 'normPhaseSpikes','phaseBinEdges',...
            'tDelay', 'allFitMdl');

        % get phase bin mids, for plotting
        phaseBinStarts = phaseBinEdges(1:(end-1));
        phaseBinEnds = phaseBinEdges(2:end);
        phaseBinMids = (phaseBinStarts + phaseBinEnds)/2;
        phaseBinMids = rad2deg(phaseBinMids);

        % convert to radians
        phaseBinMidsRads = deg2rad(phaseBinMids);

        % preallocate
        thisAmps = zeros(length(tDelay),1); % num delays 

        % loop over all time offsets
        for j = 1:length(tDelay)

            thisAmps(j) = abs(2*allFitMdl{j}.a);
        end

        % add this fly's amps to all
        allAmps = cat(2,allAmps,thisAmps);
    end

    % compute mean and SEMacross flies
    meanAmpsAllFlies = mean(allAmps,2);
    semAmpsAllFlies = std(allAmps, [], 2) / sqrt(numFlies);


    % plotting
    figure;
    c = colormap('lines');

    hold on;

    % plot individual flies
    for k = 1:numFlies
        plot(tDelay, allAmps(:,k), ...
            'LineWidth', 1, 'Color', c(k,:));
    end

    % plot mean across flies
%     plot(tDelay, meanAmpsAllFlies, ...
%         'LineWidth', 2, 'Color', 'k');

    % label x axis
    xlabel('Temporal offset (s)');
    % label y axis
    ylabel('Modulation amplitude (Hz)');

    yScale = ylim;
    % line at t = 0
    line([0 0], yScale, 'Color', 'k', 'LineWidth', 1);

    xlim([tDelay(1) tDelay(end)]);

end