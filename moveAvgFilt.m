% moveAvgFilt.m
%
% Function to perform moving average using specified time window to average
%  over.
%
% INPUTS:
%   in - input time trace to perform moving average on
%   sampRate - sampling rate of input
%   avgWindow - time window to average over, in seconds
%
% OUTPUTS:
%   out - the filtered input
%
% CREATED: 3/7/19 HHY
% UPDATED: 3/7/19 HHY
%

function out = moveAvgFilt(in, sampRate, avgWindow)
    
    % averaging window, in samples
    windowSize = round(avgWindow * sampRate);
    
    % filter parameters for filtfilt
    b = (1/windowSize)*ones(1,windowSize);
    a = 1;
    
    % perform moving average, use filtfilt so no shift introduced
    out = filtfilt(b,a, in);
end