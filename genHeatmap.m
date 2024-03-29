% heatmapMoveCondData.m
%
% Function to generate heat map of movement conditioned data generated by
%  moveCondPairData(). Similar to scatterMoveCondData() but as heat map,
%  plots the value of the dependent variable z against two independent
%  variables, x and y. Each value in heatmap is a mean of the values that
%  contributed to it.
% Specify x, y, and z data by field name. Can run any combination of
%  fictrac and imaging data, though some aren't very meaningful.
% Plots either multiple flies on the same heatmap or individual flies each
%  as own figure.
% Plots either only data during moving bouts or all data.
% Can specify temporal offset between z data and x,y data. 
% Specify number of bins as well as min and max values for all variables.
%
% INPUTS:
%   condPairData - array of structs, output of moveCondPairData()
%   xDataName - string specifying name of data to put on x-axis;
%       independent variable 1
%   yDataName - string specifying name of data to put on y-axis;
%       independent variable 2
%   zDataName - string specifying name of data to represent as color on
%       heat map; dependent variable; if this is 'counts', heat map is
%       instead count of number of values in each bin (relates x and y data)
%   xScale - 3 element vector of [xmin, xmax, xNumBins]
%   yScale - 3 element vector of [ymin, ymax, yNumBins]
%   zScale - z axis limits, as 2 element vector
%   minNumVals - minimum number of values that need to be in a bin for it
%       to be plotted in heat map, per fly
%   offset - offset between zData and x/y Data, positive values are zData
%       after x/yData, negative values are zData before x/yData; in units
%       of samples
%   samePlot - binary for whether to put all flies on same plot (1) or to
%       generate a separate figure for each fly (0)
%   plotNotMove - binary for whether to plot not moving data points; (1)
%       for yes, (0) for no; nothing distinguishes moving from not moving
%       points in heat map if (1)
%   degPerMM - if not empty, will convert any fictrac variables using mm to
%       deg, using this conversion factor
%   ttl - plot title
%
% OUTPUTS:
%   f - handle to figure(s)
%   heatmapMat - matrix of heat map values (not color indicies) that went
%       into plots
%   countsMat - number of values per bin
%
% CREATED: 9/13/19 - HHY
%
% UPDATED: 9/19/19 - HHY
%   9/27/19 - HHY - allow acceleration as a fictrac variable, change how
%       units handled to deal appropriately
%

function [f, heatmapMat, countsMat] = genHeatmap(xDat, yDat, zDat,...
    xDataName, yDataName, zDataName, xScale, yScale, zScale, minNumVals,...
    offsets, degPerMM, ttl)

    colorRes = 256; % number of different colors in heatmap, must be even
    colorbarWhite = zScale(2)/2; % using redblue colorbar, what value is white
    nanVal = [0.5 0.5 0.5]; % value for bins without enough data
    numTicks = 11; % number of ticks for colorbar
    
    % fictrac behavioral variables
    behVars = {'fwdVel', 'slideVel', 'yawAngVel', 'yawAngSpd', ...
        'totAngSpd', 'fwdAcc', 'slideAcc', 'yawAngAcc', 'totAngAccMag'};
    behVarsUnits = {'mm/s', 'mm/s', 'deg/s', 'deg/s', 'deg/s', ...
        'mm/s^2', 'mm/s^2', 'deg/s^2', 'deg/s^2'};
    % imaging variables
    ephysVars = {'spikeRate', 'medFiltV'};
    
    % preallocate matrix that will become heat map, dimension 3 is number
    %  of flies, rows = y, col = x
    heatmapMat = zeros(yScale(3), xScale(3), length(offsets));
    % preallocate matrix that counts number of samples that go into each
    %  bin
    countsMat = zeros(size(heatmapMat));
    
    % (xmax-xmin)/xNumBins
    xBinWidth = (xScale(2) - xScale(1)) / xScale(3);
    % start values of x bins
    xBinStarts = xScale(1):xBinWidth:(xScale(2) - xBinWidth);
    % end values of x bins
    xBinEnds = (xScale(1) + xBinWidth):xBinWidth:xScale(2);
    % midpoints of x bins (for plotting)
    xBinMids = (xBinStarts + xBinEnds) / 2;
    
    % same for y bins
    yBinWidth = (yScale(2) - yScale(1)) / yScale(3);
    yBinStarts = yScale(1):yBinWidth:(yScale(2) - yBinWidth);
    yBinEnds = (yScale(1) + yBinWidth):yBinWidth:yScale(2);
    yBinMids = (yBinStarts + yBinEnds) / 2;
    
        % get x data, assumes xDataName is field of ephysSpikes or fictrac
    if (strcmpi(ephysVars{1}, xDataName))
        xUnits = 'Hz';
    elseif (strcmpi(ephysVars{2}, xDataName))
        xUnits = 'mV';
    else
        % find which behavioral variable it is
        xBehVarInd = find(strcmpi(behVars, xDataName));
        % if the x data is in mm and the user desires a conversion
        %  to degrees
        if ~isempty(degPerMM) && ...
                ~isempty(strfind(behVarsUnits{xBehVarInd}, 'mm'))
            xDat = xDat .* degPerMM;
            xUnits = 'deg/s';
        else
            xUnits = behVarsUnits{xBehVarInd};
        end
    end
    % get y data, assumes yDataName is field of img or fictrac
    if (strcmpi(ephysVars{1}, yDataName))
        yUnits = 'Hz';
    elseif (strcmpi(ephysVars{2}, yDataName))
        yUnits = 'mV';
    else
        yBehVarInd = find(strcmpi(behVars, yDataName));
        % if the y data is in mm and the user desires a conversion
        %  to degrees
        if ~isempty(degPerMM) && ...
                ~isempty(strfind(behVarsUnits{yBehVarInd}, 'mm'))
            yDat = yDat .* degPerMM;
            yUnits = 'deg/s';
        else
            yUnits = behVarsUnits{yBehVarInd};
        end
    end
    zIsCount = 0; % boolean for whether z-data is count
    % get z data, assumes zDataName is field of ephysSpikes or fictrac
    if (strcmpi(ephysVars{1}, zDataName))
        zUnits = 'Hz';
    elseif (strcmpi(ephysVars{2}, zDataName))
        zUnits = 'mV';
    elseif (any(strcmpi(behVars, zDataName)))
        zBehVarInd = find(strcmpi(behVars, zDataName));
        % if the z data is in mm and the user desires a conversion
        %  to degrees
        if ~isempty(degPerMM) && ...
                strfind(behVarsUnits{zBehVarInd}, 'mm')
            zDat = zDat .* degPerMM;
            zUnits = 'deg/s';
        else
            zUnits = behVarsUnits{zBehVarInd};
        end
    elseif (strcmpi('counts', zDataName))
        zDat = ones(size(yDat));
        zIsCount = 1;
    end
    
    origXDat = xDat;
    origYDat = yDat;
    origZDat = zDat;
    
    % loop through all offsets
    for i = 1:length(offsets)
        xDat = origXDat;
        yDat = origYDat;
        zDat = origZDat;
        
        % introduce offset to data
        xDat = circshift(xDat, offsets(i));
        yDat = circshift(yDat, offsets(i));
        % delete appropriate number of elements from xDat, yDat, zDat, so
        %  wrapped elements aren't used
        if (offsets(i) < 0) % negative offsets, remove from end
            xDat = xDat(1:(end + offsets(i)));
            yDat = yDat(1:(end + offsets(i)));
            zDat = zDat(1:(end + offsets(i)));
        elseif (offsets(i) > 0) % positive offsets, remove from front
            xDat = xDat((offsets(i) + 1):end);
            yDat = yDat((offsets(i) + 1):end);
            zDat = zDat((offsets(i) + 1):end);
        end        
        
        % loop through all data points, put each one in appropriate bin
        for j = 1:length(xDat)
            % point has all valid values (not NaN)
            if (~isnan(xDat(j)) && ~isnan(yDat(j)) && ~isnan(zDat(j)))
                % find index of bin, going from both directions of start
                %  and end of bin
                whichXBinStart = find(xDat(j) >= xBinStarts, 1, 'last');
                whichXBinEnd = find(xDat(j) < xBinEnds, 1, 'first');
            
                % if either is empty, xDat point exceeds x limits
                if (~isempty(whichXBinStart) && ~isempty(whichXBinEnd))
                    whichXBin = whichXBinStart;
                else
                    whichXBin = [];
                end
                
                % find index of bin, going from both directions of start
                %  and end of bin
                whichYBinStart = find(yDat(j) >= yBinStarts, 1, 'last');
                whichYBinEnd = find(yDat(j) < yBinEnds, 1, 'first');
                
                % if either is empty, yDat point exceeds y limits
                if (~isempty(whichYBinStart) && ~isempty(whichYBinEnd))
                    whichYBin = whichYBinStart;
                else
                    whichYBin = [];
                end
                
                % put data point in correct bin
                if (~isempty(whichXBin) && ~isempty(whichYBin))
                    % sum in value, will later use for averaging
                    heatmapMat(whichYBin, whichXBin, i) = ...
                        heatmapMat(whichYBin, whichXBin, i) + zDat(j);
                    % update counter, will use for averaging and filtering
                    countsMat(whichYBin, whichXBin, i) = ...
                        countsMat(whichYBin, whichXBin, i) + 1;
                end
               
            end 
        end

    end
    

    % only when data to be plotted is not counts; get heat map values to
    %  plot, scaling
    if ~zIsCount
        % get mean value for each heatmap bin
        % NaN for every bin that has fewer than min number of data points
        %  allowed
        countsMat(countsMat < minNumVals) = nan;
        heatmapMat = heatmapMat ./ countsMat;
    
        % zScale
        colorScale = parula(colorRes);
        % shift so only 1 middle value is white
        colorScale(2:(colorRes/2 + 1), :) = colorScale(1:(colorRes/2),:);
        % index 1 corresponds to color for NaN
        colorScale(1,:) = nanVal;
        minColorInd = 2; % since ind 1 is NaN color, min ind for a color is 2

        % max difference from value set to white
        maxAmp = max(abs(zScale - colorbarWhite));
        % use maxAmp to determine scale factor between zScale and colorScale:
        % maxAmp corresponds to 1/2 of colorScale; c per z
        zcScFctr = ((colorRes - minColorInd)/2) / maxAmp;

        % generate heat maps in colorbar indicies for plotting
        % conversion b/w z and color
        heatmapPlotMat = zcScFctr .* (heatmapMat - colorbarWhite) + ...
            (colorRes/2 + 1);

        % values that exceed zScale boundaries in either direction are set
        %  to min/max color values
        heatmapPlotMat(heatmapPlotMat < minColorInd) = minColorInd;
        heatmapPlotMat(heatmapPlotMat > colorRes) = colorRes;

        % NaNs in heatmap set to 1, to plot as color for bins without
        %  enough data
        heatmapPlotMat(isnan(heatmapPlotMat)) = 1;

        % scale factors for colorbar
        % exclude 1 from showing on colorbar
        colorbarLims = [minColorInd colorRes]; 

        % where ticks are, on zScale; will become tick labels
        tickLabels = zScale(1):((zScale(2)-zScale(1))/(numTicks - 1)):zScale(2);
        % where ticks are, in indicies
        tickLocs = zcScFctr .* (tickLabels - colorbarWhite) + ...
            (colorRes/2 + 1);
    
    % plot counts; set up color scale
    else
        % zScale
        colorScale = parula(colorRes);
        minColorInd = 1; % shift corresponding to 1 vs 0 indexing
        
        zcScFctr = (colorRes - minColorInd) / (zScale(2)-zScale(1));
        
            % conversion b/w z and color
            heatmapPlotMat = zcScFctr .* (countsMat - zScale(1)) + ...
                minColorInd;

        
        % scale factors for colorbar
        colorbarLims = [minColorInd colorRes]; 

        % where ticks are, on zScale; will become tick labels
        tickLabels = zScale(1):((zScale(2)-zScale(1))/(numTicks - 1)):zScale(2);
        % where ticks are, in indicies
        tickLocs = zcScFctr .* (tickLabels - zScale(1)) + minColorInd;
    end
    
    % plot heatmap(s)
    for i = 1:size(heatmapPlotMat,3)
        f(i) = figure;
        colormap(colorScale);
        
        imgHandle = imagesc(xBinMids, yBinMids, heatmapPlotMat(:,:,i));
        % invert axes so smaller values on y axis toward bottom
        ax = gca;
        ax.YDir = 'normal';
        % this property needs to be set here and not as option in imagesc
        %  initial call (doesn't work there)
        imgHandle.CDataMapping = 'direct';
        axis square;
        
        % axis labels
        xlabel(sprintf('%s (%s)', xDataName, xUnits));
        ylabel(sprintf('%s (%s)', yDataName, yUnits));
        
        % colorbar, scaled and labeled appropriately for z-axis data
        colorbarHandle = colorbar('LimitsMode', 'manual', 'Limits',...
            colorbarLims, 'TicksMode', 'manual', 'Ticks', tickLocs, ...
            'TickLabelsMode', 'manual', 'TickLabels', tickLabels);
        if ~zIsCount
            colorbarHandle.Label.String = sprintf('%s (%s)', ...
                zDataName, zUnits);
        else
            colorbarHandle.Label.String = zDataName;
        end
        
        % plot title
        if ~zIsCount
            ttlStr = sprintf('%s offset %d %s vs. %s and %s', ttl, ...
                offsets(i), zDataName, xDataName, yDataName);

        else
            ttlStr = sprintf('%s %s %s vs. %s', ttl, ...
                zDataName, yDataName, xDataName);
        end
        
        title(ttlStr); 
    end
    
end