% plotSpikerateMoveNotMove.m
%
% Function to plot the mean spike rate during moving and not moving periods
%  (output of getMoveNotMoveSpikerate()). 
%
% INPUTS:
%   datDir - folder containing getXCorrEphysContParam_cell() output files
%   yScale - scale for y-axis
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 9/2/23 - HHY
%
% UPDATED:
%   9/2/23 - HHY
%
function diffSpikerate = plotSpikerateMoveNotMove(datDir)

    % get file through GUI
    [outputFNames, outputPath] = uigetfile('*.mat', ...
        'Select getMoveNotMoveSpikerate() file', ...
        datDir, 'MultiSelect', 'off');

    fullFilePath = [outputPath filesep outputFNames];

    load(fullFilePath, 'meanSpikerate');

    % get mean, SEM across flies
    meanAllFlies = mean(meanSpikerate, 1);
    SEMAllFlies = std(meanSpikerate, [], 1) / sqrt(size(meanAllFlies,1));

    % get difference b/w moving and not moving (move - not move)
    diffSpikerate = meanSpikerate(:,1) - meanSpikerate(:,2);
    meanDiff = mean(diffSpikerate);
    semDiff = std(diffSpikerate) / sqrt(length(diffSpikerate));

    % x vector for plotting move and not move 
    xVec = 1:2;
    xVec = xVec';

    % plot move and not move
    figure;
    c = colormap('lines');

    hold on;

    % plot for each individual fly
    plot(xVec, meanSpikerate, ...
        'Marker', 'x','LineWidth',0.5, 'Color', c(1,:));

    % plot mean across flies
    errorbar(xVec, meanAllFlies, SEMAllFlies, ...
        '_','LineWidth', 2, 'CapSize', 0, 'Color', c(2,:));

    xScale = xlim;
    xScale(1) = xScale(1) - (0.5 * (xVec(end)-xVec(1)));
    xScale(2) = xScale(2) + (0.5 * (xVec(end)-xVec(1)));
    xlim(xScale);

    % label x-axis
    xticks(xVec);
    xticklabels({'Moving', 'Not Moving'});

    % label y-axis
    ylabel('Spike Rate (Hz)');


    % plot difference (move - not move)
    figure;
    hold on;

    % plot for each individual fly
    plot(ones(length(diffSpikerate),1), diffSpikerate, ...
        '.', 'Marker', 'x','LineWidth',0.5, 'Color', c(1,:));

    % plot mean across flies
    errorbar(1, meanDiff, semDiff, ...
        '_','LineWidth', 2, 'CapSize', 0, 'Color', c(2,:));

    % line at y = 0
    line(xlim, [0 0], 'LineWidth', 1, 'Color', 'k');

    % label x-axis
    xticks(1);
    xticklabels('Move - Not Move');

    ylabel('Spike Rate Difference (Hz)');

end