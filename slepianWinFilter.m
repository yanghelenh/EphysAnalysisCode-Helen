% slepianWinFilter.m
%
% Implements low pass filtering by convolving input with 0th order Slepian
%  window with user specified half bandwidth (in Hz) and sample rate 
%  (in Hz)
% 
% Based off of py.filter_tools.lpf_linear_filter in a2lib from S. Holtz,
%  but corrected so input amplitude isn't changed by filtering
%
% See dpss documentation in matlab
%
% INPUT:
%   in - input to filter
%   bandwidth - bandwidth of Slepian window, in Hz
%   sampRate - sample rate (in Hz) of input (and Slepian window)
%
% OUTPUT:
%   out - input filtered by convolution with Slepian window
%
% CREATED: 8/27/19 - HHY
%
% UPDATED: 8/27/19 - HHY
%

function [out] = slepianWinFilter(in, bandwidth, sampRate)

    % length of Slepian window: 1/halfBandwidth * sampRate
    winLen = round(1/(bandwidth * 0.5) * sampRate);
    
    % set NW = 1, return only 1st sequence
    slepianWin = dpss(winLen, 1, 1);
    
    % normalize Slepian window
    normSlepWin = slepianWin ./ (sum(slepianWin));
    
    % filter input
    out = conv(in, normSlepWin', 'same');
end