% fwhm.m
%
% Function to compute full width at half max (FWHM) given x and y. X is
%  what width will be computed on; y defines peak.
%
% INPUTS:
%   x - vector that defines width of peak
%   y - vector that defines height of peak
%   valInd - indices into x/y in which to find peak
%   minOrMax - 'max' or 'min' for whether to find fwhm of a maximum or
%       minimum
%
% OUTPUTS:
%   fwhmOut - FWHM for this peak, in units of x
%
% CREATED: 9/1/23 - HHY
%
% UPDATED:
%   9/1/23 - HHY
%   10/6/23 - HHY - allow only range of indices to be considered as
%   location for max; also works on min
%
function fwhmOut = fwhm(x, y, valInd, minOrMax)
    
    if (strcmpi(minOrMax,'max'))
        % get index of max
        [maxVal, maxInd] = max(y(valInd));
        maxInd = valInd(1) + maxInd - 1;
    
        % get half-max value
        halfMaxVal = maxVal / 2;
    
        % going from max index backwards, find first index that falls below
        %  halfMaxVal
        backInd = maxInd:-1:1;
        yLHMBackInd = backInd(find(y(backInd) < halfMaxVal, 1, 'first'));
    
        % going from max index forwards, find first index that falls below
        %  halfMaxVal
        fwdInd = maxInd:length(y);
        yLHMFwdInd = fwdInd(find(y(fwdInd) < halfMaxVal, 1, 'first'));
    
    
    %     % get start and end indices of y that are greater than the half-max val
    %     yGHMStartInd = find(y > halfMaxVal,1,'first');
    %     yGHMEndInd = find(y > halfMaxVal, 1, 'last');
    
        % correct if start and end are at edges, use ends of x vector
        if (isempty(yLHMBackInd) || yLHMBackInd == 1)
            xValStart = x(1);
        else
    %         % use linear interpolation to get x value of end
    %         xValStart = interp1(y((yGHMStartInd-1):yGHMStartInd),...
    %             x((yGHMStartInd-1):yGHMStartInd), halfMaxVal, 'linear');
            % use linear interpolation to get x value of end
            xValStart = interp1(y(yLHMBackInd:(yLHMBackInd+1)),...
                x(yLHMBackInd:(yLHMBackInd+1)), halfMaxVal, 'linear');
        end
    
        if (isempty(yLHMFwdInd) || yLHMFwdInd == length(x))
            xValEnd = x(end);
        else
    %         % use linear interpolation to get x value of end
    %         xValEnd = interp1(y(yGHMEndInd:(yGHMEndInd + 1)),...
    %             x(yGHMEndInd:(yGHMEndInd + 1)), halfMaxVal, 'linear');
            xValEnd = interp1(y((yLHMFwdInd-1):yLHMFwdInd),...
                x((yLHMFwdInd-1):yLHMFwdInd), halfMaxVal, 'linear');
        end
    else
        % get index of min
        [minVal, minInd] = min(y(valInd));
        minInd = valInd(1) + minInd - 1;

        % get max value
        maxVal = max(y);
    
        % get half-min value
        halfMinVal = (maxVal - minVal) / 2;
    
        % going from min index backwards, find first index that goes above
        %  halfMinVal
        backInd = minInd:-1:1;
        yLHMBackInd = backInd(find(y(backInd) > halfMinVal, 1, 'first'));
    
        % going from min index forwards, find first index that goes above
        %  halfMinVal
        fwdInd = minInd:length(y);
        yLHMFwdInd = fwdInd(find(y(fwdInd) > halfMinVal, 1, 'first'));
    
    
    %     % get start and end indices of y that are greater than the half-max val
    %     yGHMStartInd = find(y > halfMaxVal,1,'first');
    %     yGHMEndInd = find(y > halfMaxVal, 1, 'last');
    
        % correct if start and end are at edges, use ends of x vector
        if (isempty(yLHMBackInd) || yLHMBackInd == 1)
            xValStart = x(1);
        else
    %         % use linear interpolation to get x value of end
    %         xValStart = interp1(y((yGHMStartInd-1):yGHMStartInd),...
    %             x((yGHMStartInd-1):yGHMStartInd), halfMaxVal, 'linear');
            % use linear interpolation to get x value of end
            xValStart = interp1(y(yLHMBackInd:(yLHMBackInd+1)),...
                x(yLHMBackInd:(yLHMBackInd+1)), halfMinVal, 'linear');
        end
    
        if (isempty(yLHMFwdInd) || yLHMFwdInd == length(x))
            xValEnd = x(end);
        else
    %         % use linear interpolation to get x value of end
    %         xValEnd = interp1(y(yGHMEndInd:(yGHMEndInd + 1)),...
    %             x(yGHMEndInd:(yGHMEndInd + 1)), halfMaxVal, 'linear');
            xValEnd = interp1(y((yLHMFwdInd-1):yLHMFwdInd),...
                x((yLHMFwdInd-1):yLHMFwdInd), halfMinVal, 'linear');
        end
    end

    % fwhm is difference between x val for start and end
    fwhmOut = xValEnd - xValStart;
end