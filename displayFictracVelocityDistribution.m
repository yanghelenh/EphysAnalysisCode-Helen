% displayFictracVelocityDistribution.m
%
% Function that prompts user for trial.mat file for a fictracOnly type
%  experiment and processes that FicTrac data to display a heatmap
%  representing the distribution of counts across yaw and forward 
%  velocities.
% Used to check FicTrac calibration, as existing data suggests there is a
%  rightward skew. As of 10/9/20
% Displays 3D scatterplot showing point cloud around (0,0,0) in yaw,
%  forward and slide axes
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
%   10/12/20 - HHY
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
    xScale = [-50 50 25];
    yScale = [-50 50 25];
    zScale = [0 5000]; % check this
    minNumVals = 20;
    offsets = 0;
    
    xDataName = 'yawAngVel';
    yDataName = 'fwdVel';
    zDataName = 'counts';
    

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
        genHeatmap(fictracProc.(xDataName)', fictracProc.(yDataName)', [],...
            xDataName, yDataName, zDataName, xScale, yScale, zScale, ...
            minNumVals, offsets, fictracParams.degPerMM, ftFileName);
        
        % plot 3D scatterplot, shows point cloud in all 3 dimensions
        figure;
        scatter3(fictracProc.yawAngVel, ...
            fictracProc.fwdVel*fictracParams.degPerMM, ...
            fictracProc.slideVel*fictracParams.degPerMM, ...
            'Marker', '.', ...
            'MarkerFaceAlpha', 0.2, 'MarkerEdgeAlpha', 0.2);
        axLims = [xScale(1) xScale(2)];
        xlim(axLims);
        ylim(axLims);
        zlim(axLims);
        xlabel('yawAngVel (deg/s)');
        ylabel('fwdVel (deg/s)');
        zlabel('slideVel (deg/s)');
        axis square;
        
        
    else % otherwise, display error message and return
        disp('Selected trial.mat file does not contain fictrac data');
        return;
    end
end