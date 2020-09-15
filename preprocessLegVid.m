% preprocessLegVid.m
%
% Function for preprocessing raw voltage signals from the leg tracking
%  camera (strobe signal with falling edge as frame start). Extracts leg
%  video frame times. Checks that number of frames camera says it has 
%  captured matches number of trigger pulses sent. 
% Returns these parameters as outputs in leg struct.
% Adapted from function of same name in 2PAnalysisCode-Helen
%
%
% INPUTS:
%   legVidPath - full path to leg video mp4
%   daqData - struct of data from experimental DAQ, processed by
%       preprocessUserDaq()
%   daqOutput - struct of output signals sent by experimental DAQ,
%       processed by preprocessUserDaq() 
%   daqTime - vector of times corresponding to each sample point of daqData
%
% OUTPUTS:
%   leg - output struct with following fields:
%   	frameTimes - start times of each leg video frame, in seconds
%   	trigTimes - start times of each leg video frame trigger, in seconds
%       vidHeight - height in pixels of leg video
%       vidWidth - width in pixels of leg video
%       numVidFrames - number of video frames in leg video file
%
% CREATED: 8/7/20 HHY
%
% UPDATED: 
%   8/7/20 - HHY
%   9/6/20 - HHY - convert output to leg struct 
%   9/14/20 - HHY - include info about video frame size

function leg = preprocessLegVid(legVidPath, daqData, daqOutput, daqTime)

    % frame start indicies, strobe signal - find falling edges
    frameStarts = find(diff(daqData.legCamFrames) < -0.1);
    
    % frame trigger indicies, output sent by experimental DAQ
    frameTrigs = find(diff(daqOutput.legCamFrameStartTrig) > 0.1);
    
%     % check that number of captured frames matches number of triggered
%     %  frames
%     if (length(frameStarts) ~= length(frameTrigs))
%         disp('Warning: Frame count mismatch in leg tracking video');
%     end
    
    % leg vid frame times
    legVidFrameTimes = daqTime(frameStarts + 1);
    
    % leg vid trigger times
    legVidTrigTimes = daqTime(frameTrigs + 1);
    
    % save into leg output struct
    leg.frameTimes = legVidFrameTimes;
    leg.trigTimes = legVidTrigTimes;
    
    % video reader for leg video
    legVidReader = VideoReader(legVidPath);
    
    % save some parameters about the video
    leg.vidHeight = legVidReader.Height;
    leg.vidWidth = legVidReader.Width;
    leg.numVidFrames = legVidReader.NumFrames;
    
end