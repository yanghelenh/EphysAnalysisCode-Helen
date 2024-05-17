% plotAvgEphysVYawVelBinFwdAcc.m
%
% Function to plot output of getCorrelationEphysFwdAccYawVel_cond_all() or
%  getSlideWinEphysFwdAccYawVel_cond_all(). 
% Plots yaw velocity on x, normalized spike rate on y, binned by forward
%  acceleration (i.e. different lines for different forward acceleration
%  bins). Mean spike rate for that yaw velocity+forward acceleration bin
% Each forward acceleration bin in a different color. Thin lines are
%  individual flies (if plotted) and thick lines are means across flies
%
% INPUTS:
%   datDir - path to folder containing output files
%   indivFlies - plot individual flies or not boolean
%   numBinsYaw - number of bins of yaw velocity
%   yawVelRange - range of yaw velocity values to plot and bin across. 
%     As 2 element vector for start and end
%   fwdAccBins - matrix specifying forward acceleration bins where each row
%     is the start (col 1) and end (col 2) of one bin, and the number of
%     rows is the number of bins
%   minNumVals - minimum number of values that must be in bin for it to
%       have a value
%   yScale - scale of y-axis
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 5/15/24 - HHY
%
% UPDATED:
%   5/15/24 - HHY
%
function plotAvgEphysVYawVelBinFwdAcc(datDir, indivFlies, numBinsYaw, ...
    yawVelRange, fwdAccBins, minNumVals, yScale)

    % prompt user to select getCorrelationEphysFwdAccYawVel_cond_all() or
    %  getSlideWinEphysFwdAccYawVel_cond_all() files
    [corrEphysFNames, corrEphysPath] = uigetfile('*.mat', ...
        'Select output file', datDir, 'MultiSelect', 'on');

    % if only 1 file selected, not cell array; make sure loop still
    %  works 
    % num flies is number of files
    if (iscell(corrEphysFNames))
        numFlies = length(corrEphysFNames);
    else
        numFlies = 1;
    end

    % get boundaries of yaw velocity bins
    yawBinSize = (yawVelRange(2) - yawVelRange(1)) / numBinsYaw;
    yawBinEdges = yawVelRange(1):yawBinSize:yawVelRange(2);
    yawBinStarts = yawBinEdges(1:(end-1));
    yawBinEnds = yawBinEdges(2:end);
    yawBinMids = (yawBinStarts + yawBinEnds)/2;

    % number of forward acceleration bins
    numFwdAccBins = size(fwdAccBins,1);

    % preallocate
    allFliesAvg = zeros(numFlies, length(yawBinMids), numFwdAccBins);
    allFliesErr = zeros(numFlies, length(yawBinMids), numFwdAccBins);


    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(corrEphysFNames))
            outName = corrEphysFNames{i};
        else
            outName = corrEphysFNames;
        end
   

        % load data from cond_bout file
        fullFilePath = [corrEphysPath filesep corrEphysFNames{i}];
    
        % load variables
        load(fullFilePath, 'ephysValsNorm', 'fwdAccVals', ...
            'yawVelVals', 'ephysParam', 'tDelay');

        % hack to flip L/R for 220907_fly01
        if (contains(outName, '220907_fly01'))
            yawVelVals = yawVelVals * -1;
        end


        % find mean and SEM for each bin
        % loop through all forward acceleration bins
        for j = 1:numFwdAccBins
            % get all values that fall into this forward acceleration bin
            thisFwdBinLog = (fwdAccVals >= fwdAccBins(j,1)) & ...
                (fwdAccVals < fwdAccBins(j,2));

            % loop through all yaw velocity bins
            for k = 1:numBinsYaw
                % get all values that fall into this yaw velocity bin
                thisYawBinLog = (yawVelVals >= yawBinStarts(k)) & ...
                (yawVelVals < yawBinEnds(k));

                % logical for values that fall into both fwd and yaw bin
                thisBinLog = thisFwdBinLog & thisYawBinLog;

                % get ephys values for this bin
                thisEphys = ephysValsNorm(thisBinLog);

                % remove NaNs
                thisEphys(isnan(thisEphys)) = [];

                % number of values in this bin
                numVals = length(thisEphys);

                % if number of values in this bin exceeds minimum, compute
                %  mean for bin and add to output matrix
                if (numVals > minNumVals)
                    % get mean and SEM for this bin
                    allFliesAvg(i,k,j) = mean(thisEphys);
                    allFliesErr(i,k,j) = std(thisEphys) / ...
                        sqrt(length(thisEphys));
                else % otherwise, NaN
                    allFliesAvg(i,k,j) = nan;
                    allFliesErr(i,k,j) = nan;
                end

            end
        end
    end


    % mean and error across flies, preallocate
    totAvg = zeros(numBinsYaw, numFwdAccBins);
    totErr = zeros(numBinsYaw, numFwdAccBins);

    % get mean and SEM across flies
    for i = 1:numFwdAccBins
        for j = 1:numBinsYaw
            thisBin = allFliesAvg(~isnan(allFliesAvg(:,j,i)),j,i);
            totAvg(j,i) = mean(thisBin);
            totErr(j,i) = std(thisBin) / sqrt(length(thisBin));
        end
    end

    % generate figure
    h = zeros(1,numFwdAccBins);
    figure;

    c = colormap('lines');

    hold on;


    for i = 1:numFwdAccBins
        if (indivFlies)
        % plot individual flies, if flagged
        plot(yawBinMids,squeeze(allFliesAvg(:,:,i)), 'Color', c(i,:));

        % plot mean across flies - thick line
        h(i) = plot(yawBinMids, totAvg(:,i), 'Color', c(i,:), ...
            'LineWidth', 2);
        else
            h(i) = plot_err_patch_v2(yawBinMids, totAvg(:,i), ...
                totErr(:,i), c(i,:)*0.5, c(i,:));
        end
    end

    xlim(yawVelRange);
    ylim(yScale);

    % generate strings for legend (different fwd acceleration bins)
    legendStr = {};
    for i = 1:numFwdAccBins
        legendStr{i} = sprintf('FwdAcc %.1f to %.1f mm/s^2', ...
            fwdAccBins(i,1), fwdAccBins(i,2));
    end

    % ephys param label
    yLblStr = sprintf('%s, delay = %d ms', ephysParam, tDelay * 1000);
    ylabel(yLblStr);

    xlabel('Yaw Velocity (deg/s)');

    legend(h, legendStr);
end