% lowPassFiler.m
%
% LOWPASSFILTER Filters the data using a lowpass butterworth filter of 
%  specified cut out
%   
% INPUTS:
%   data - trace to be filtered
%   lowPassCutOff - value (Hz) that will be the top limit of the filter
%   sampleRate - rate data is sampled at to allow correct conversion in to
%     Hz.
% 
% OUTPUT:
%   out - filtered version of data
% 
% CREATED: Yvette Fisher 1/2018
% UPDATED: 1/31/19 HHY
%

function [out] = lowPassFilter(data, lowPassCutOff, sampleRate)

% low pass filter the data trace :
% fprintf('\nLow-pass filtering at %d Hz\n',lowPassCutOff);

% build 2nd order butter function
[b,a] = butter(2 , lowPassCutOff / (sampleRate/2), 'low');

% filter data using butterworth function
out = filtfilt(b, a, data);

end

