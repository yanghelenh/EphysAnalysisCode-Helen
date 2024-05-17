% plotHeatmapEphysVYawVelBinFwdAcc.m
%
% Function to plot output of getCorrelationEphysFwdAccYawVel_cond_all() or
%  getSlideWinEphysFwdAccYawVel_cond_all() as heatmap. 
% Plots yaw velocity on x, forward acceleration on y, normalized spike rate
%  in color. Mean spike rate for that yaw velocity+forward acceleration bin
% When plot individual, separate heatmap for each fly
%
% INPUTS:
%   datDir - path to folder containing output files
%   indivFlies - plot individual flies or not boolean
%   numBinsYaw - number of bins of yaw velocity
%   yawVelRange - range of yaw velocity values to plot and bin across 
%     As 2 element vector for start and end
%   numBinsFwd - number of bins of forward acceleration
%   fwdAccRange - range of forward acceleration values to plot and bin
%     across
%   minNumVals - minimum number of values that must be in bin for it to
%       have a value
%   zScale - color scale, as 2 element vector
%
% OUTPUTS:
%   none, but generates plot
%
% CREATED: 5/15/24 - HHY
%
% UPDATED:
%   5/15/24 - HHY
%
function plotHeatmapEphysVYawVelBinFwdAcc(datDir, indivFlies, numBinsYaw, ...
    yawVelRange, numBinsFwd, fwdAccRange, minNumVals, zScale)

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

    % get boundaries of forward acceleration bins
    fwdBinSize = (fwdAccRange(2) - fwdAccRange(1)) / numBinsFwd;
    fwdBinEdges = fwdAccRange(1):fwdBinSize:fwdAccRange(2);
    fwdBinStarts = fwdBinEdges(1:(end-1));
    fwdBinEnds = fwdBinEdges(2:end);
    fwdBinMids = (fwdBinStarts + fwdBinEnds)/2;


    % preallocate
    allFliesAvg = zeros(numFlies, numBinsYaw, numBinsFwd);
    allFliesErr = zeros(numFlies, numBinsYaw, numBinsFwd);


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
        for j = 1:numBinsFwd
            % get all values that fall into this forward acceleration bin
            thisFwdBinLog = (fwdAccVals >= fwdBinStarts(j)) & ...
                (fwdAccVals < fwdBinEnds(j));

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

    % mean across flies, preallocate
    totAvg = zeros(numBinsYaw, numBinsFwd);

    % get mean and SEM across flies
    for i = 1:numBinsFwd
        for j = 1:numBinsYaw
            thisBin = allFliesAvg(~isnan(allFliesAvg(:,j,i)),j,i);
            totAvg(j,i) = mean(thisBin);
        end
    end

    % color scale params
    colorRes = 256; % number of different colors in heatmap, must be even
    nanVal = [0.5 0.5 0.5]; % value for bins without enough data
    numTicks = 11; % number of ticks for colorbar

    % get color scale for heatmap
    colorScale = parula(colorRes);
    % index 1 corresponds to color for NaN
    colorScale = [nanVal; colorScale];
    minColorInd = 2; % since ind 1 is NaN color, min ind for a color is 2

    % min to max of zScale to color bar indices
    zRange = zScale(2) - zScale(1);
    % scale factor b/w color resolution and zRange
    zcScFctr = colorRes / zRange;

    % scale factors for colorbar
    % exclude 1 from showing on colorbar
    colorbarLims = [minColorInd colorRes]; 

    % where ticks are, on zScale; will become tick labels
    tickLabels = zScale(1):(zRange/(numTicks - 1)):zScale(2);
    % where ticks are, in indicies
    tickLocs = zcScFctr .* tickLabels + 1;


    % plot heatmap for mean across flies
    figure;

    % generate heat map in colorbar indicies for plotting
    % conversion b/w z and color
    heatmapPlotMat = zcScFctr .* totAvg' + 1;

    % values that exceed zScale boundaries in either direction are set
    %  to min/max color values
    heatmapPlotMat(heatmapPlotMat < minColorInd) = minColorInd;
    heatmapPlotMat(heatmapPlotMat > colorRes) = colorRes;

    % NaNs in heatmap set to 1, to plot as color for bins without
    %  enough data
    heatmapPlotMat(isnan(heatmapPlotMat)) = 1;

    colormap(colorScale);
    
    imgHandle = imagesc(yawBinMids, fwdBinMids, heatmapPlotMat);
    % invert axes so smaller values on y axis toward bottom
    ax = gca;
    ax.YDir = 'normal';
    % this property needs to be set here and not as option in imagesc
    %  initial call (doesn't work there)
    imgHandle.CDataMapping = 'direct';
    axis square;
    
    % axis labels
    xlabel('Yaw Velocity (deg/s)');
    ylabel('Fwd Acceleration (mm/s2)');
    
    % colorbar, scaled and labeled appropriately for z-axis data
    colorbarHandle = colorbar('LimitsMode', 'manual', 'Limits',...
        colorbarLims, 'TicksMode', 'manual', 'Ticks', tickLocs, ...
        'TickLabelsMode', 'manual', 'TickLabels', tickLabels);
    colorbarHandle.Label.String = ephysParam;


    % plot heatmap for individual flies
    if (indivFlies)
        for i = 1:numFlies
            figure;

            % get this fly's name
            if (iscell(corrEphysFNames))
                outName = corrEphysFNames{i};
            else
                outName = corrEphysFNames;
            end
            flyName = outName(1:19);

            % this fly's matrix for plotting
            thisMat = squeeze(allFliesAvg(i,:,:))';
        
            % generate heat map in colorbar indicies for plotting
            % conversion b/w z and color
            heatmapPlotMat = zcScFctr .* thisMat + 1;
        
            % values that exceed zScale boundaries in either direction are set
            %  to min/max color values
            heatmapPlotMat(heatmapPlotMat < minColorInd) = minColorInd;
            heatmapPlotMat(heatmapPlotMat > colorRes) = colorRes;
        
            % NaNs in heatmap set to 1, to plot as color for bins without
            %  enough data
            heatmapPlotMat(isnan(heatmapPlotMat)) = 1;
        
            colormap(colorScale);
            
            imgHandle = imagesc(yawBinMids, fwdBinMids, heatmapPlotMat);
            % invert axes so smaller values on y axis toward bottom
            ax = gca;
            ax.YDir = 'normal';
            % this property needs to be set here and not as option in imagesc
            %  initial call (doesn't work there)
            imgHandle.CDataMapping = 'direct';
            axis square;
            
            % axis labels
            xlabel('Yaw Velocity (deg/s)');
            ylabel('Fwd Acceleration (mm/s^2)');
            
            % colorbar, scaled and labeled appropriately for z-axis data
            colorbarHandle = colorbar('LimitsMode', 'manual', 'Limits',...
                colorbarLims, 'TicksMode', 'manual', 'Ticks', tickLocs, ...
                'TickLabelsMode', 'manual', 'TickLabels', tickLabels);
            colorbarHandle.Label.String = ephysParam;

            title(flyName);
        end
    end

end