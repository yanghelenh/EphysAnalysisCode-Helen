% plotSpikesStepPhase.m
%
% Function to plot timing of spikes relative to leg phase, at different
%  time offsets between ephys and behavior. I.e., plots the output of
%  getSpikesStepPhase_cell(), but combining multiple flies.
% One figure per delay; 6 subplots, 1 for each leg per figure
%
% INPUTS:
%   datDir - directory with output files
%   whichDelays - which delays between ephys and behavior to plot
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 9/3/23 - HHY
%
% UPDATED:
%   9/3/23 - HHY
%
function plotSpikesStepPhase(datDir, whichDelays)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

    % initialize
    allFliesSpikes = [];

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

    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(outputFNames))
            outName = outputFNames{i};
        else
            outName = outputFNames;
        end
        
        outputFullPath = [outputPath filesep outName];

        load(outputFullPath, 'normPhaseSpikes','phaseBinEdges',...
            'tDelay', 'spikeCounts', 'phaseCounts');

        % normalize by subtracting by value at bin 36
        for j = 1:length(tDelay)
            for k = 1:6
                normPhaseSpikes(j,:,k) = normPhaseSpikes(j,:,k) - ...
                    normPhaseSpikes(j,36,k);
            end
        end

        % get running matrix of all flies normPhaseSpikes
        allFliesSpikes = cat(4,allFliesSpikes, normPhaseSpikes);
%         allFliesSpikes = cat(4,allFliesSpikes, spikeCounts);
%         allFliesSpikes = cat(4,allFliesSpikes, phaseCounts);
    end

    % compute mean and SEM across flies
    meanSpikes = mean(allFliesSpikes,4);
    allFliesSEMs = std(allFliesSpikes, [], 4) / sqrt(numFlies);

    % get phase bin mids, for plotting
    phaseBinStarts = phaseBinEdges(1:(end-1));
    phaseBinEnds = phaseBinEdges(2:end);
    phaseBinMids = (phaseBinStarts + phaseBinEnds)/2;


    % plotting, loop over all tDelays
    for i = 1:length(whichDelays)
        thisDelayInd = find(whichDelays(i) == tDelay);

        figure;
        c = colormap('lines');

        for j = 1:6
            
            subplot(3,2,subInd(j));
            hold on;

            % plot individual flies
            for k = 1:numFlies
                plot(phaseBinMids, squeeze(allFliesSpikes(thisDelayInd,:,j, k)), ...
                    'LineWidth', 1, 'Color', c(1,:));
            end

            % plot mean across flies
            plot(phaseBinMids, squeeze(meanSpikes(thisDelayInd,:,j)), ...
                'LineWidth', 2, 'Color', c(2,:));

            % label x axis
            if (j == 3 || j == 6)
                xlabel('Step phase (deg)');
            end

            xlim([0 360]);
        end

        sgtitle(sprintf('Delay = %.2f',whichDelays(i)));
    end


end