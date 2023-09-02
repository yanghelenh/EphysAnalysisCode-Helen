% fwhm.m
%
% Function to compute full width at half max (FWHM) given x and y. X is
%  what width will be computed on; y defines peak.
%
% INPUTS:
%   x - vector that defines width of peak
%   y - vector that defines height of peak
%
% OUTPUTS:
%   fwhmOut - FWHM for this peak, in units of x
%
% CREATED: 9/1/23 - HHY
%
% UPDATED:
%   9/1/23 - HHY
%
function fwhmOut = fwhm(x, y)
    % get half-max value
    halfMaxVal = max(y) / 2;

    % get start and end indices of y that are greater than the half-max val
    yGHMStartInd = find(y > halfMaxVal,1,'first');
    yGHMEndInd = find(y > halfMaxVal, 1, 'last');

    % correct if start and end are at edges, use ends of x vector
    if (yGHMStartInd == 1 || isempty(yGHMStartInd))
        xValStart = x(1);
    else
        % use linear interpolation to get x value of end
        xValStart = interp1(y((yGHMStartInd-1):yGHMStartInd),...
            x((yGHMStartInd-1):yGHMStartInd), halfMaxVal, 'linear');
    end

    if (yGHMEndInd == length(x) || isempty(yGHMEndInd))
        xValEnd = x(end);
    else
        % use linear interpolation to get x value of end
        xValEnd = interp1(y(yGHMEndInd:(yGHMEndInd + 1)),...
            x(yGHMEndInd:(yGHMEndInd + 1)), halfMaxVal, 'linear');
    end

    % fwhm is difference between x val for start and end
    fwhmOut = xValEnd - xValStart;
end