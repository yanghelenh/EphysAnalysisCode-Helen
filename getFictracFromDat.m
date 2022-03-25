% getFictracFromDat.m
%
% Function to extract FicTrac data from .dat file output by FicTrac instead
%  of from DAQ output
% 
% INPUTS:
%   datFilePath - full path to .dat file
%   daqData - struct of data from experimental DAQ, processed by
%       preprocessUserDaq()
%   daqTime - 
%
% OUTPUTS:
%   fictrac - struct of outputs with following fields:
%       yawAngVel - yaw angular velocity, in deg/sec
%       yawAngPosWrap - yaw angular position (heading), wrapped between 0 
%           and 360 deg
%       fwdVel - forward velocity, in mm/sec
%       fwdCumPos - distance traveled in the forward direction, with trial
%           start being 0 mm (signed)
%       slideVel - slide/lateral velocity, in mm/sec
%       slideCumPos - distance traveled in the lateral direction, with
%           trial start being 0 mm (signed)
%       xPos - x position, in mm, if fly were walking on xy plane
%       yPos - y position, in mm, if fly were walking on xy plane
%       t - time in seconds of each timepoint of the above vectors
%
% CREATED: 10/4/21 - HHY
%
% UPDATED:
%   10/4/21 - HHY
% 
function fictrac = getFictracFromDat(datFilePath, daqData, daqTime)

    % frame start indicies, strobe signal - find falling edges
    frameStarts = find(diff(daqData.ficTracCamFrames) < -0.1);
    
    % fictrac vid frame times
    ficTracVidFrameTimes = daqTime(frameStarts + 1);

    
    % save into fictrac output struct
    fictrac.t = ficTracVidFrameTimes;
    
    sampRate = 1/median(diff(ficTracVidFrameTimes));

    % read in file
    fileID = fopen(datFilePath, 'r');

    % format spec for this file, 23 comma separated values
    formatSpec = '%d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, \n';

    % read in whole file
    sizeOutput = [23 Inf];

    fileOutput = fscanf(fileID, formatSpec, sizeOutput);

    fclose(fileID);
    
    % get outputs
    yawPos = unwrap(fileOutput(17, :)); % need to unwrap heading
    fwdPos = fileOutput(20,:);
    slidePos = fileOutput(21, :);

    % constants
    % lowpass filter cutoff (approximately half of FicTrac sample rate)
    LOWPASS_FILTER_CUTOFF = 40; % in Hz 
    MAX_YAW_VELOCITY = 2500; % deg/sec
    MAX_FWD_VELOCITY = 2500; % deg/sec
    MAX_SLIDE_VELOCITY = 2500; % deg/sec
    BALL_DIAM = 6.46; % diameter of ball, in mm
    
    % yaw/heading    
    [yawAngVel, yawAngPos] = ficTracDatPosVel(yawPos,...
        sampRate, LOWPASS_FILTER_CUTOFF, MAX_YAW_VELOCITY);
    
    % wrap yaw angular position to 360 deg instead of it being cumulative
    yawAngPosWrap = wrapTo360(yawAngPos);
    
    % conversion factor between degrees and mm
    circum = BALL_DIAM * pi; % circumference of ball, in mm
    mmPerDeg = circum / 360; % mm per degree of ball
    
    % forward direction
    [fwdAngVel, fwdAngPos] = ficTracDatPosVel(fwdPos,...
        sampRate, LOWPASS_FILTER_CUTOFF, MAX_FWD_VELOCITY);
    fwdVel = fwdAngVel .* mmPerDeg; % velocity in mm/sec
    % cumulative forward position in mm, where start of trial is at 0
    fwdCumPos = (fwdAngPos - fwdAngPos(1)) .* mmPerDeg; 
    
    % slide direction    
    [slideAngVel, slideAngPos] = ficTracDatPosVel(slidePos,...
        sampRate, LOWPASS_FILTER_CUTOFF, MAX_SLIDE_VELOCITY);
    
    slideVel = slideAngVel .* mmPerDeg; % velocity in mm/sec
    % cumulative slide position in mm, where start of trial is at 0
    slideCumPos = (slideAngPos-slideAngPos(1)) .* mmPerDeg;     
    
    % position incorporating heading - as if fly were walking on x-y plane,
    %  x-y coordinates at each time point
    % start with fly at (0,0) and facing 0 deg
    zeroedYawAngPos = yawAngPos - yawAngPos(1); 
    
    % movement in x (in degrees) at each time point
    xChangePos = (fwdAngVel ./ sampRate) .* sind(zeroedYawAngPos) + ...
        (slideAngVel ./ sampRate) .* sind(zeroedYawAngPos + 90);  

    % x position in mm (i.e. x-coordinate of fly's position at each time 
    %  point), starts at 0
    xPos = (cumsum(xChangePos) - xChangePos(1)) .* mmPerDeg;
   
    % movement in y (in degrees) at each time point
    yChangePos = (fwdAngVel ./ sampRate) .* cosd(zeroedYawAngPos) + ...
        (slideAngVel ./ sampRate) .* cosd(zeroedYawAngPos + 90);

    % y position in mm (i.e. y-coordinate of fly's position at each time 
    %  point), starts at 0
    yPos = (cumsum(yChangePos) - yChangePos(1)) .* mmPerDeg;
    
    % put all output variables into output fictrac struct
    fictrac.yawAngVel = yawAngVel;
    fictrac.yawAngPosWrap = yawAngPosWrap;
    fictrac.fwdVel = fwdVel;
    fictrac.fwdCumPos = fwdCumPos;
    fictrac.slideVel = slideVel;
    fictrac.slideCumPos = slideCumPos;
    fictrac.xPos = xPos;
    fictrac.yPos = yPos;
    
end