% xcorrWGaps.m
%
% Function to compute cross-correlation when there are non-random gaps in
%  the data.
% Operates on 2 input vectors. Gaps in vectors should be NaNs. Input
%  vectors must have same sampling rate.
% Uses corr function to calculate Pierson's correlation for each lag, after
%  removing missing data points.
%
% INPUTS:
%   x - first input vector
%   y - second input vector
%   numLags - number of lags to calculate, should be odd number (0 plus
%       equal number of lags on each side of 0)
%
% OUTPUTS:
%   allCorr - all correlations, of length numLags  
%   lags - all lags used, in samples, matched to allCorr. Negative lags are
%       later time points in x correlated with earlier time points in y.
%
% CREATED: 8/31/23 - HHY
%
% UPDATED:
%   8/31/23 - HHY
%
function [allCorr, lags] = xcorrWGaps(x, y, numLags)

    % convert x and y to columns
    if (~iscolumn(x))
        x = x';
    end
    if (~iscolumn(y))
        y = y';
    end

    % get number of lags to each side of 0 (will add lags if numLags is
    %  even)
    numLags1Side = ceil((numLags-1)/2);

    % preallocate 
    allCorr = zeros(numLags1Side * 2 + 1, 1);

    % get lags
    negLags = fliplr(1:numLags1Side) * -1;
    posLags = 1:numLags1Side;

    lags(1:numLags1Side) = negLags;
    lags(numLags1Side + 1) = 0;
    lags((numLags1Side + 2):(numLags1Side * 2 + 1)) = posLags;
    lags = lags'; % as column vector

    % loop through negative lags
    for i = 1:numLags1Side
        thisLag = negLags(i) * -1;

        % get x and y with this lag
        % for negative lag, remove elements at front of x from
        %  consideration, so later time points in x aligned to earlier ones
        %  in y
        thisX = x;
        thisX(1:thisLag) = [];

        thisY = y;
        thisY((end-thisLag+1):end) = [];

        % remove any NaNs
        % logical of NaNs for x and y
        xNaNLog = isnan(thisX);
        yNaNLog = isnan(thisY);

        % logical true for NaN at index for either x or y
        allNaNLog = xNaNLog | yNaNLog;

        % remove NaNs
        thisX(allNaNLog) = [];
        thisY(allNaNLog) = [];

        if (~isempty(thisX) && ~isempty(thisY))
            % get correlation b/w x and y
            allCorr(i) = corr(thisX, thisY);
        else
            allCorr(i) = nan;
        end
    end

    % loop through positive lags
    for i = 1:numLags1Side
        thisLag = posLags(i);

        % get x and y with this lag
        % for positive lag, remove elements at front of y from
        %  consideration, so later time points in y aligned to earlier ones
        %  in x
        thisY = y;
        thisY(1:thisLag) = [];

        thisX = x;
        thisX((end-thisLag+1):end) = [];

        % remove any NaNs
        % logical of NaNs for x and y
        xNaNLog = isnan(thisX);
        yNaNLog = isnan(thisY);

        % logical true for NaN at index for either x or y
        allNaNLog = xNaNLog | yNaNLog;

        % remove NaNs
        thisX(allNaNLog) = [];
        thisY(allNaNLog) = [];

        if (~isempty(thisX) && ~isempty(thisY))
            % get correlation b/w x and y
            allCorr(i + numLags1Side + 1) = corr(thisX, thisY);
        else
            allCorr(i + numLags1Side + 1) = nan;
        end
    end  

    % correlation for lag = 0
    thisX = x;
    thisY = y;

    % remove any NaNs
    % logical of NaNs for x and y
    xNaNLog = isnan(thisX);
    yNaNLog = isnan(thisY);

    % logical true for NaN at index for either x or y
    allNaNLog = xNaNLog | yNaNLog;

    % remove NaNs
    thisX(allNaNLog) = [];
    thisY(allNaNLog) = [];

    if (~isempty(thisX) && ~isempty(thisY))
        allCorr(numLags1Side + 1) = corr(thisX, thisY);
    else
        allCorr(numLags1Side + 1) = nan;
    end
end