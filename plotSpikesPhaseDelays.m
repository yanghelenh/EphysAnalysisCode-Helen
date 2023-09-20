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
function plotSpikesPhaseDelays(datDir)

    % legs to subplot indices
    % puts left legs on left, and front legs on top
    subInd = [2 4 6 1 3 5]; 

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
            'tDelay', 'spikeCounts', 'phaseCounts');

        % get phase bin mids, for plotting
        phaseBinStarts = phaseBinEdges(1:(end-1));
        phaseBinEnds = phaseBinEdges(2:end);
        phaseBinMids = (phaseBinStarts + phaseBinEnds)/2;

        % convert to radians
        phaseBinMidsRads = deg2rad(phaseBinMids);

        % preallocate
        thisAmps = zeros(length(tDelay),6); % num delays x legs

        % loop over all time offsets
        for j = 1:length(tDelay)
            % loop over all legs
            for k = 1:6
                thisNormPhaseSpikes = squeeze(normPhaseSpikes(j,:,k));

                % get median value as starting point for offset
                medVal = median(thisNormPhaseSpikes);

                % get sine fit to this delay and leg
                fitMdl = fit2PiPerSine(phaseBinMidsRads, ...
                    thisNormPhaseSpikes, [rand(), rand(), medVal]);

                % get amplitude from sine fit
                curAmp = abs(2*fitMdl.a);

                % add to thisAmps matrix
                thisAmps(j,k) = curAmp;
            end
        end

        % add this fly's amps to all
        allAmps = cat(3,allAmps,thisAmps);
    end

    % compute mean and SEMacross flies
    meanAmpsAllFlies = mean(allAmps,3);
    semAmpsAllFlies = std(allAmps, [], 3) / sqrt(numFlies);


    % plotting
    figure;
    c = colormap('lines');

    % loop over all legs
    for j = 1:6
        
        subplot(3,2,subInd(j));
        hold on;

        % plot individual flies
        for k = 1:numFlies
            plot(tDelay, squeeze(allAmps(:,j,k)), ...
                'LineWidth', 1, 'Color', c(1,:));
        end

        % plot mean across flies
        plot(tDelay, squeeze(meanAmpsAllFlies(:,j)), ...
            'LineWidth', 2, 'Color', c(2,:));

        % label x axis
        if (j == 3 || j == 6)
            xlabel('Temporal offset (s)');
        end

        yScale = ylim;
        % line at t = 0
        line([0 0], yScale, 'Color', 'k', 'LineWidth', 1);

        xlim([tDelay(1) tDelay(end)]);
    end

    sgtitle('Amplitude of modulation (Hz)');

end