% ficTracSignalDecoding.m
%
% FICTRACSIGNALDECODING takes a fictrac position value and extracts ball 
%  velocity. This function will take a FicTrac output signal as aquired by 
%  the DAQ as an analog signal and solve for the ball's angular velocity in
%  the given dimention. To do this the signal is then UNWRAPPED to handle 
%  the abrupt transitions caused by when the ball rotates completely and 
%  the signal resets (0->10 volts or 10 volts -> 0 transistions). Then the
%  signal will be further CLEANED to remove extra position values 
%  surrounding those signal reset time points. Then the position signal is 
%  LOW PASS FILTERED two times!! to remove noise/ sampling of fictrac.
%  Then the velocity of the ball with be solved for in degree/s by using 
%  the gradient function, and taking into consideration the sample rate the
%  data was collected at. Velocity values above what a resonable fly would 
%  turn the ball at (maxFlyVelocity) are discarded
%  
% INPUTS:
%   ficTracBallPosition - array containing data from 0-10 volts relating
%     to the ball's position
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
% CREATED: 1/2018 YEF
% UPDATED: 11/24/18 HHY
% 

function [velocityOut, accumulatedPositionOut] = ficTracSignalDecoding(...
    ficTracBallPosition, sampleRate, lowPassFilterCutOff, maxFlyVelocity)

FICTRAC_MAX_VOLTAGE = 10;  % volts

% transfrom ficTrac signal into radians  
posRadians = ficTracBallPosition .* 2 .* pi ./ FICTRAC_MAX_VOLTAGE; 

% upwrap position signal
unwrappedPos = unwrap(posRadians);

% find indexes where the unwrapping happened (tolerance = pi)
upwrappedIndexes = find (abs(diff(posRadians)) > pi); 

NUM_SAMPLES_FROM_WRAP_TO_REPLACE = 2;
% handle edge case so we don't fill off the edge of the trace
upwrappedIndexes = upwrappedIndexes(upwrappedIndexes > ...
    NUM_SAMPLES_FROM_WRAP_TO_REPLACE & ...
    upwrappedIndexes < (length(unwrappedPos) - ...
    NUM_SAMPLES_FROM_WRAP_TO_REPLACE) ); 

cleanedPos = unwrappedPos;
% replace potentially problematic indexes with Nan
for i = 1:length(upwrappedIndexes)
    index_start = upwrappedIndexes(i) -  NUM_SAMPLES_FROM_WRAP_TO_REPLACE; 
    index_end = upwrappedIndexes(i) +  NUM_SAMPLES_FROM_WRAP_TO_REPLACE; 
    
    cleanedPos(index_start:index_end) = NaN;
end

% replace NaN values with the last preceding value that was a real number
nanIDX = find(isnan(cleanedPos)); % find NaN indexes
% replace with preceeding value
while(~isempty(nanIDX))
    cleanedPos(nanIDX) = cleanedPos(nanIDX - 1);
    
    % find any remaining NaN
    nanIDX  = find(isnan(cleanedPos));
end

% low pass filter the position array
filteredPosition = lowPassFilter(cleanedPos, lowPassFilterCutOff, ...
    sampleRate);
% low pass filter the position array again to be more aggressive
filteredPosition = lowPassFilter(filteredPosition, lowPassFilterCutOff,...
    sampleRate);

 
% plotting to check how well unwrapping, cleaning and filtering worked
% can be commented out once you are happy with the parameters
%{
figure('Position',[50, 50, 1000, 400]);  set(gcf, 'Color', 'w');
ax(1) = subplot(4, 1, 1);
plot( posRadians ); hold on;
ylabel ('rad');
ax(2) = subplot (4, 1, 2);
plot( unwrappedPos ); hold on;
ax(3) = subplot (4, 1, 3);
plot( cleanedPos ); hold on;
ax(4) = subplot (4, 1, 4);
plot( filteredPosition ); hold on;
linkaxes(ax,'x');
%}

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

% plotting to check degree calulation and velocity
%{
figure('Position',[50, 50, 1000, 400]);  set(gcf, 'Color', 'w');
bx(1) = subplot(3, 1, 1);
plot( filteredPosition ); hold on;
bx(1) = subplot(3, 1, 2);
plot( accumulatedPositionOut ); hold on;
ylabel ('deg');
bx(2) = subplot (3, 1, 3);
plot( velocityOut ); hold on;
linkaxes(bx,'x');
%}

end

