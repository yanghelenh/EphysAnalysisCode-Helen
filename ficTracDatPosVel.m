% ficTracDatPosVel.m
%
% Extracts FicTrac position and velocity from unwrapped position output
%  from .dat file
%  
% INPUTS:
%   unwrappedPos - unwrapped position in radians
%   sampleRate - Rate the data was aquired at (samples/second)
%   lowPassFilterCutOff - frequency that the position signal will be low
%     pass filtered at (Hz)
%   maxFlyVelocity - max value of realistic fly movement (deg/s) 
%
% OUTPUT:
%   velocityOut -array containing ball's instentanous velocity (degree/sec)
%   accumulatedPositionOut - array containing the filtered and unwraped
%     position signal
% 
% CREATED: 10/4/21 - HHY
% UPDATED: 10/4/21 - HHY
% 

function [velocityOut, accumulatedPositionOut] = ficTracDatPosVel(...
    unwrappedPos, sampleRate, lowPassFilterCutOff, maxFlyVelocity)

% low pass filter the position array
filteredPosition = lowPassFilter(unwrappedPos, lowPassFilterCutOff, ...
    sampleRate);
% low pass filter the position array again to be more aggressive
filteredPosition = lowPassFilter(filteredPosition, lowPassFilterCutOff,...
    sampleRate);


% transform from radians into degrees, send to user
accumulatedPositionOut = (filteredPosition / (2*pi)) * 360;

% take derivative and adjust for sample rate to solve for deg/s
%velocityOut = diff( accumulatedPositionOut ) .* sampleRate ; % degees/sec
velocityOut = gradient(accumulatedPositionOut) .* sampleRate ; % degees/sec

%low pass filter the velocity signal
%velocityOut = lowPassFilter(velocityOut, lowPassFilterCutOff, sampleRate);

% remove velocity values that are too large to be possible for the fly
velocityOut = replaceValuesOutsideThresholdBound(velocityOut, ...
    maxFlyVelocity);

end

