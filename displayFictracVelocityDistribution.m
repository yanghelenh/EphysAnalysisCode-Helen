% displayFictracVelocityDistribution.m
%
% Function that prompts user for trial.mat file for a fictracOnly type
%  experiment and processes that FicTrac data to display a heatmap
%  representing the distribution of counts across yaw and forward 
%  velocities.
% Used to check FicTrac calibration, as existing data suggests there is a
%  rightward skew. As of 10/9/20
%
% INPUTS:
%   none, but prompts user for trial.mat file
%
% OUTPUTS:
%   none, but plots distribution heatmap
%
% CREATED: 10/9/20 - HHY
%
% UPDATED:
%   10/9/20 - HHY
%
function displayFictracVelocityDistribution()
    % some constants
    fictracParams.dsf = 20; % downsample to 1000 Hz;
    fictracParams.filtParams.padLen = int32(200);
    fictracParams.filtParams.sigmaPos = int32(100); % 100 ms
    fictracParams.filtParams.sigmaVel = int32(50); % 50 ms
    
    BALL_DIAM = 6.46;
    circum = BALL_DIAM * pi; % circumference of ball, in mm
    fictracParams.mmPerDeg = circum / 360; % mm per degree of ball
    fictracParams.degPerMM = 360 / circum; % deg per mm ball
    
    % plotting constants
    xScale = [-250 250 25];
    yScale = [-5 15 30];
    zScale = [0 2000]; % check this
    minNumVals = 20;
    offsets = 0;
    

    disp('Select trial.mat file from FicTrac experiment');
    [ftFileName, ftDir] = uigetfile('*.mat');
    ftFile = [ftDir filesep ftFileName];
    
    % get metadata file in same folder as trial file
    mdFile = [ftDir filesep 'metaDat.mat'];
    
    % load trial.mat file
    load(ftFile, 'inputParams', 'rawData', 'rawOutput');
    % load settings from metadata file
    load(mdFile, 'settings');
    
    % check that this file contains FicTrac data
    if(contains(inputParams.exptCond, 'fictrac','IgnoreCase', true))
        % preprocess data
        [daqData, daqOutput, daqTime] = preprocessUserDaq(inputParams, ...
            rawData, rawOutput, settings);
        fictrac = preprocessFicTrac(daqData, daqTime, ...
            settings.bob.sampRate);
        
        % assumes FicTrac didn't drop 
        fictrac.dropInd = [];
        
        % filter and downsample FicTrac data
        fictracProc = dsFiltFictrac(fictracParams, fictrac);
        
        % plot heatmap
        genHeatmap(fictracProc.yawAngVel', fictracProc.fwdVel', [],...
            'yawAngVel', 'fwdVel', 'counts', xScale, yScale, zScale, ...
            minNumVals, offsets, fictracParams.degPerMM, ftFileName);
        
    else % otherwise, display error message and return
        disp('Selected trial.mat file does not contain fictrac data');
        return;
    end
end