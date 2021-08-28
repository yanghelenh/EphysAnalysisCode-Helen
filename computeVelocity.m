% computeVelocity.m
%
% Simple function to lowpass filter position and compute velocity

function vel = computeVelocity(pos, sampleRate)
    [b,a] = butter(2 , 0.5, 'low');

    % filter data using butterworth function
    filtPos = filtfilt(b, a, pos);
    filtPos = filtfilt(b, a, filtPos);
    % transform from radians into degrees
    filtPosDeg = (filtPos / (2*pi)) * 360;
    % velocity
    vel = gradient(filtPosDeg) .* sampleRate;
end