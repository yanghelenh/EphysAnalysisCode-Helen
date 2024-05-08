% madRmOutliers.m
%
% Function that removes outliers from vector where outliers are greater
%  than x MADs away from the median, where x is a user specified parameter
%
% INPUTS:
%   inVec - input vector
%   numMADs - number of MADs away to consider as outlier
%
% OUTPUTS:
%   outVec - output vector, which is input vector with outliers removed
%
% CREATED: 5/7/24 - HHY
%
% UPDATED:
%   5/7/14 - HHY
%
function outVec = madRmOutliers(inVec, numMADs)
    % input has to be vector, throw error otherwise
    if isvector(inVec)
        % numMADs has to be a postive number, throw error otherwise
        if (numMADs < 0)
            error('Threshold has to be a positive number of MADs');
        else
            % get MAD and median
            thisMAD = mad(inVec,1);
            thisMedian = median(inVec);

            % upper and lower threshold
            upThresh = thisMedian + thisMAD * numMADs;
            lowThresh = thisMedian - thisMAD * numMADs;

            % logical for which values of input are valid (between lower
            %  and upper thresholds)
            validLog = (inVec >= lowThresh) & (inVec <= upThresh);

            % get output
            outVec = inVec(validLog);
        end
    else
        error('Input has to be a vector');
    end
end