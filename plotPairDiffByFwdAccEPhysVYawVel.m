% plotPairDiffByFwdAccEPhysVYawVel.m
%
% Function to plot output of getCorrelationEphysFwdAccYawVel_cond_all() or
%  getSlideWinEphysFwdAccYawVel_cond_all() as difference in average ephys
%  value between two forward acceleration bins, in specific yaw velocity 
%  bins
% I.e. for each fly, average ephys value in one forward acceleration + yaw
%  velocity bin - average ephys value in other forward acceleration + yaw
%  velocity bin. Specify 2 forward acceleration bins and 1 or more yaw
%  velocity bins
%
% INPUTS:
%   datDir - path to folder containing output files
%   yawVelBins - matrix specifying yaw velocity bins where each row
%     is the start (col 1) and end (col 2) of one bin, and the number of
%     rows is the number of bins
%   fwdAccBins - matrix specifying forward acceleration bins where each row
%     is the start (col 1) and end (col 2) of one bin; must be exactly 2
%     rows
%   minNumVals - minimum number of values that must be in bin for it to
%       have a value
%   indivFlies - plot individual flies or not boolean
%   yScale - scale of y-axis
%
% OUTPUTS:
%   none, but makes plot
%
% CREATED: 5/16/24 - HHY
%
% UPDATED:
%   5/16/24 - HHY
%
function allFliesDiff = plotPairDiffByFwdAccEPhysVYawVel(datDir, yawVelBins, ...
    fwdAccBins, minNumVals, indivFlies, yScale)

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

    numYawVelBins = size(yawVelBins,1);
    numFwdAccBins = size(fwdAccBins,1);

    if(numFwdAccBins ~=2)
        disp('Number of forward acceleration bins must be 2');
        return;
    end

    % preallocate
    allFliesMeans = zeros(numFlies,numYawVelBins,numFwdAccBins);
    allFliesDiff = zeros(numFlies, numYawVelBins);


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
        load(fullFilePath, 'ephysVals', 'ephysValsNorm', 'fwdAccVals', ...
            'yawVelVals', 'ephysParam', 'tDelay');

        % hack to flip L/R for 220907_fly01
        if (contains(outName, '220907_fly01'))
            yawVelVals = yawVelVals * -1;
        end

        % loop through all yaw velocity bins
        for j = 1:numYawVelBins
            % get all values that fall into this yaw velocity bin
            thisYawBinLog = (yawVelVals >= yawVelBins(j,1)) & ...
                (yawVelVals < yawVelBins(j,2));

            meanEachFwd = nan(numFwdAccBins,1);
            for k = 1:numFwdAccBins
                thisFwdBinLog = (fwdAccVals >= fwdAccBins(k,1)) & ...
                    (fwdAccVals < fwdAccBins(k,2));
                
                % logical for values that fall into both fwd and yaw bin
                thisBinLog = thisFwdBinLog & thisYawBinLog;

                % get ephys values for this bin
                thisEphys = ephysValsNorm(thisBinLog);
%                 thisEphys = ephysVals(thisBinLog);

                % remove NaNs
                thisEphys(isnan(thisEphys)) = [];

                % number of values in this bin
                numVals = length(thisEphys);

                % get mean value in this bin
                if (numVals >= minNumVals)
                    meanEachFwd(k) = mean(thisEphys);
                end
            end

            % get difference between forward acceleration bins for this yaw
            %  velocity (2nd - 1st)
            allFliesDiff(i,j) = meanEachFwd(2) - meanEachFwd(1);

            allFliesMeans(i,j,1) = meanEachFwd(1);
            allFliesMeans(i,j,2) = meanEachFwd(2);
        end
    end

    meanAllFlies = nan(numYawVelBins,1);
    semAllFlies = nan(numYawVelBins,1);
    % get mean and SEM across flies
    for i = 1:numYawVelBins
        thisBin = allFliesDiff(~isnan(allFliesDiff(:,i)),i);
        meanAllFlies(i) = mean(thisBin);
        semAllFlies(i) = std(thisBin) / sqrt(numFlies);
    end

    % find mean and SEM difference across flies
%     meanAllFlies = mean(allFliesDiff,1);
%     semAllFlies = std(allFliesDiff,1) / sqrt(numFlies);

    % mean and error across flies, preallocate
    totAvg = zeros(numYawVelBins, numFwdAccBins);
    totErr = zeros(numYawVelBins, numFwdAccBins);

    % get mean and SEM across flies
    for i = 1:numFwdAccBins
        for j = 1:numYawVelBins
            thisBin = allFliesMeans(~isnan(allFliesMeans(:,j,i)),j,i);
            totAvg(j,i) = mean(thisBin);
            totErr(j,i) = std(thisBin) / sqrt(length(thisBin));
        end
    end


    % x vector
    xVec = 1:numYawVelBins;
    xVec = xVec';

    
    % plot
    figure;
    c = colormap('lines');
    hold on;

    if (indivFlies)
        % plot indivdual flies
        plot(xVec,allFliesDiff', 'Marker', '.','LineWidth',0.5, ...
            'Color', c(1,:));
        % plot mean
        plot(xVec,meanAllFlies,'Marker', '_', 'Color', c(2,:));
    else
        % plot mean
        errorbar(xVec, meanAllFlies, semAllFlies, ...
            '_','LineWidth', 2, 'CapSize', 0, 'Color', c(2,:));
    end

    xScale = xlim;
    xScale(1) = xScale(1) - (0.5 * (xVec(end)-xVec(1)));
    xScale(2) = xScale(2) + (0.5 * (xVec(end)-xVec(1)));
    xlim(xScale);

    % line at y=0
    line(xlim, [0 0], 'LineWidth', 1, 'Color', 'k');

    if (size(yScale,1)~=1)
        ylim(yScale(1,:));
    else
        ylim(yScale);
    end


    % label x-axis
    xticks(xVec);
    yawBinLbl = {};
    for i = 1:numYawVelBins
        yawBinLbl{i} = sprintf('%d to %d', yawVelBins(i,1), ...
            yawVelBins(i,2));
    end

    xticklabels(yawBinLbl);

    % label axes
    xlabel('Yaw velocity bins (deg/s)');
    ylabel('Difference in normalized spike rate');

    ttlStr = sprintf('Forward Acceleration %d to %d - %d to %d mm/s^2', ...
        fwdAccBins(2,1), fwdAccBins(2,2), fwdAccBins(1,1), fwdAccBins(1,2));
    title(ttlStr);


    % plot
    if (numYawVelBins == 1)
        % x vector
        xVec = 1:numFwdAccBins;
        xVec = xVec';
    
        figure;
        c = colormap('lines');
        hold on;
    
        if (indivFlies)
            % plot indivdual flies
            for i = 1:numFlies
                plot(xVec,squeeze(allFliesMeans(i,:,:)), 'Marker', '.','LineWidth',0.5, ...
                    'Color', c(1,:));
            end
            % plot mean
            plot(xVec,totAvg,'Marker', '_', 'Color', c(2,:));
        else
            % plot mean
            errorbar(xVec, totAvg, totErr, ...
                '_','LineWidth', 2, 'CapSize', 0, 'Color', c(2,:));
        end
    
        xScale = xlim;
        xScale(1) = xScale(1) - (0.5 * (xVec(end)-xVec(1)));
        xScale(2) = xScale(2) + (0.5 * (xVec(end)-xVec(1)));
        xlim(xScale);
    
        ylim(yScale(2,:));
    
        % label x-axis
        xticks(xVec);
        fwdBinLbl = {};
        for i = 1:numFwdAccBins
            fwdBinLbl{i} = sprintf('%d to %d', fwdAccBins(i,1), ...
                fwdAccBins(i,2));
        end
    
        xticklabels(fwdBinLbl);
    
        % label axes
        xlabel('Forward Acceleration bins (mm/s^2)');
        ylabel('Normalized spike rate');
    end
end