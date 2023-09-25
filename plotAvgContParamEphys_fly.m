% plotAvgContParamEphys_fly.m
%
% Function to plot output of getCorrelationEphysContParam_cond_all() as
%  average across flies and/or individual flies
% Spike rate on x, behavior on y
%
% INPUTS:
%   datDir - path to folder containing savePCAinterpStepParam() output files
%   avg - 'mean' or 'median' for which type of average to plot
%   indivFlies - plot individual flies or not boolean
%   numBins - number of bins in x (spike rate)
%   xRange - range of x values to plot and bin across. As 2 element vector 
%       for start and end
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 9/24/23 - HHY
%
% UPDATED:
%   9/24/23 - HHY
%
function plotAvgContParamEphys_fly(datDir, avg, indivFlies, numBins, ...
    xRange, minNumVals, isCirc, yScale)

    % prompt user to select getCorrelationEphysContParam() files
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

    % get boundaries of bins
    binSize = (xRange(2) - xRange(1)) / numBins;
    binEdges = xRange(1):binSize:xRange(2);
    binStarts = binEdges(1:(end-1));
    binEnds = binEdges(2:end);
    binMids = (binStarts + binEnds)/2;

    % preallocate
    allFliesAvg = zeros(numFlies, length(binMids));
    allFliesErr = zeros(numFlies, length(binMids));


    for i = 1:numFlies
        % handle whether it's a cell array or not
        if (iscell(corrEphysFNames))
            outName = corrEphysFNames{i};
        else
            outName = corrEphysFNames;
        end
        
        outputFullPath = [corrEphysPath outName];

        % load data from cond_bout file
        fullFilePath = [corrEphysPath filesep corrEphysFNames{i}];
    
        % load variables
        load(fullFilePath, 'ephysVals', 'ephysValsNorm', 'behVals1D', ...
            'ephysParam', 'behParams', 'legs', 'tDelay');

%         ephysValsNorm = ephysVals;

        if (contains(outName, '220907_fly01_yawVel') || ...
                contains(outName, '220907_fly01_stepDirection'))
            behVals1D = behVals1D * -1;
        end

%         thisLog = (behVals1D >= xRange(1) ) & (behVals1D <= xRange(2));
% 
%         corrcoef(behVals1D(thisLog), ephysValsNorm(thisLog))


        % loop through all bins, find mean and SEM
        for j = 1:numBins
            % get logical for which samples fall into this bin
            thisBinLog = (ephysValsNorm >= binStarts(j)) & ...
                (ephysValsNorm < binEnds(j));
            % get ephys values for these samples
            thisBehVal = behVals1D(thisBinLog);

            thisBehVal(isnan(thisBehVal)) = [];
            numVals = length(thisBehVal);

            if (numVals > minNumVals)
                if (strcmpi(avg,'mean'))
                    if ~(isCirc)
                    % get mean and SEM for this bin
                        allFliesAvg(i,j) = mean(thisBehVal);
                        allFliesErr(i,j) = std(thisBehVal) / sqrt(length(thisBehVal));
                    else
                        thisBehVal = deg2rad(thisBehVal);
                        allFliesAvg(i,j) = wrapTo180(rad2deg(circ_mean(thisBehVal)));
                        allFliesErr(i,j) = rad2deg(circ_std(thisBehVal)) / sqrt(length(thisBehVal));
                    end
                elseif (strcmpi(avg,'median'))
                    if ~(isCirc)
                        % get median and MAD for this bin
                        allFliesAvg(i,j) = median(thisBehVal);
                        allFliesErr(i,j) = mad(thisBehVal,1);
                    else
                        % get median and MAD for this bin
                        allFliesAvg(i,j) = circ_median(thisBehVal);
                        allFliesErr(i,j) = mad(thisBehVal,1);
                    end
                end
            else
                allFliesAvg(i,j) = nan;
                allFliesErr(i,j) = nan;
            end
        end
    end



    % generate figure
    figure;

    hold on;

    totAvg = zeros(1,length(binMids));
    totErr = zeros(1, length(binMids));
    % get average and error across flies
    if (strcmpi(avg,'mean'))

        for i = 1:size(allFliesAvg,2)
%             thisBin = allFliesAvg(~isnan(allFliesAvg(:,i)),i);
            thisBin = allFliesAvg(:,i);
            if ~(isCirc)
                totAvg(i) = mean(thisBin);
                totErr(i) = std(thisBin) / sqrt(length(thisBin));
            else
                thisBin = deg2rad(thisBin);
                totAvg(i) = wrapTo180(rad2deg(circ_mean(thisBin)));
                totErr(i) = rad2deg(circ_std(thisBin)) / sqrt(length(thisBin));
            end
        end
    elseif (strcmpi(avg,'median'))
        for i = 1:size(allFliesAvg,2)
%             thisBin = allFliesAvg(~isnan(allFliesAvg(:,i)),i);
            thisBin = allFliesAvg(:,i);
            totAvg(i) = median(thisBin);
            totErr(i) = mad(thisBin);
        end
    end

    % plot average across flies
%     plot_err_patch_v2(binMids, totAvg, totErr,...
%         [0 0.4470 0.7410],[0.3010 0.7450 0.9330]);

    plot(binMids,totAvg,'Color','k','LineWidth',2);

    % plot individual flies if flagged, black lines
    if (indivFlies)
        plot(binMids,allFliesAvg');
    end

    % ephys param label
    ephysLblStr = sprintf('%s, delay = %d ms', ephysParam, tDelay * 1000);
    xlabel(ephysLblStr);

    behLblStr = [];
    if iscell(behParams)
        for i = 1:length(behParams)
            if ~i==1
                behLblStr = [behLblStr ', ' behParams{i} ' ' legs{i}];
            else
                behLblStr = [behParams{i} ' ' legs{i}];
            end
        end
    else
        behLblStr = [behParams ' ' legs];
    end
    ylabel(behLblStr);

    xlim(xRange);
    ylim(yScale);
end