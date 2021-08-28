% displayFictracVelocityDistribution_ftOutput.m
%
% Function that prompts user for the .dat output file returned by the
%  FicTrac program and processes that FicTrac data to display a heatmap
%  representing the distribution of counts across yaw and forward 
%  velocities.
% Displays 3D scatterplot showing point cloud around (0,0,0) in yaw,
%  forward and slide axes
% Adaptation of displayFictracVelocityDistribution.m for FicTrac .dat
%  output files instead of experimental file (from DAQ output)
%
% INPUTS:
%   none, but prompts user for outputdatatest.dat file
%
% OUTPUTS:
%   none, but plots distribution heatmap
%
% CREATED: 10/13/20 - HHY
%
% UPDATED:
%   10/13/20 - HHY
%

function displayFictracVelocityDistribution_ftOutput()

    % some constants
    lowPassFilterCutOff = 40;
    
    % plotting constants
    xScale = [-50 50 25];
    yScale = [-50 50 25];
    zScale = [0 500]; % check this
    minNumVals = 20;
    offsets = 0;
    
    xDataName = 'yawAngVel';
    yDataName = 'fwdVel';
    zDataName = 'counts';

    disp('Select outputdatatest.dat file from FicTrac experiment');
    [ftFileName, ftDir] = uigetfile('*.dat');
    ftFile = [ftDir filesep ftFileName];

    % read in file
    fileID = fopen(ftFile, 'r');
    
    % format spec for this file, 23 comma separated values
    formatSpec = '%d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, \n';

    % read in whole file
    sizeOutput = [23 Inf];

    fileOutput = fscanf(fileID, formatSpec, sizeOutput);

    fclose(fileID);

    % position vectors
    yawPos = unwrap(fileOutput(17, :)); % need to unwrap heading
    fwdPos = fileOutput(20,:);
    slidePos = fileOutput(21, :);

    % time
    t = fileOutput(22,:);
    t = t-t(1);
%     sampleRate = 1/median(diff(t));
    sampleRate = 80;
    
    % compute velocities, filter as in standard analysis pipeline
    yawAngVel = computeVelocity(yawPos);
    fwdVel = computeVelocity(fwdPos);
    slideVel = computeVelocity(slidePos);
    
    % plot heatmap
    genHeatmap(eval(xDataName), eval(yDataName), [],...
        xDataName, yDataName, zDataName, xScale, yScale, zScale, ...
        minNumVals, offsets, [], ftFileName);

    % plot 3D scatterplot, shows point cloud in all 3 dimensions
    figure;
    scatter3(yawAngVel, ...
        fwdVel, ...
        slideVel, ...
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
    
    
    % helper function to compute velocity in deg/s from position in radians
    function vel = computeVelocity(pos)
        [b,a] = butter(2 , 0.5, 'low');

        % filter data using butterworth function
        filtPos = filtfilt(b, a, pos);
        filtPos = filtfilt(b, a, filtPos);
        % transform from radians into degrees
        filtPosDeg = (filtPos / (2*pi)) * 360;
        % velocity
        vel = gradient(filtPosDeg) .* sampleRate;
    end
    
    
end
