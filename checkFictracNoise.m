% checkFictracNoise.m
%
% Function to load in FicTrac output files run on same video (testVid.avi
%  from 201104), outputdatatest2-4.dat. To check whether offline run of
%  FicTrac produces the same or different output (especially wrt noise)
%  when run with identical parameters on the same video).

function checkFictracNoise()
    % plotting constants
    xScale = [-25 25 25];
    yScale = [-25 25 25];
    zScale = [0 100]; % check this
    minNumVals = 20;
    offsets = 0;

    xDataName = 'yawAngVel';
    yDataName = 'fwdVel';
    zDataName = 'counts';

    disp('Select outputdatatest.dat file from FicTrac experiment');
    [ftFileNames, ftDir] = uigetfile('*.dat', 'MultiSelect', 'on');

    for i = 1:length(ftFileNames)
        ftFile = [ftDir filesep ftFileNames{i}];

        % read in file
        fileID = fopen(ftFile, 'r');

        % format spec for this file, 23 comma separated values
        formatSpec = '%d, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, \n';

        % read in whole file
        sizeOutput = [23 Inf];

        fileOutput = fscanf(fileID, formatSpec, sizeOutput);

        fclose(fileID);

        % output variables
        % preallocate variables if they don't yet exist
        if (~exist('yawPos'))
            yawPos = zeros(length(ftFileNames), size(fileOutput,2));
            fwdPos = zeros(length(ftFileNames), size(fileOutput,2));
            slidePos = zeros(length(ftFileNames), size(fileOutput,2));
            t = zeros(length(ftFileNames), size(fileOutput,2));
        end

        % position vectors
        yawPos(i,:) = unwrap(fileOutput(17, :)); % need to unwrap heading
        fwdPos(i,:) = fileOutput(20,:);
        slidePos(i,:) = fileOutput(21, :);

        % time
        t(i,:) = fileOutput(22,:);
        t(i,:) = t(i,:)-t(i,1);

    end

    % plot
%     figure;
%     plot(t',yawPos');
% 
%     figure;
%     plot(t',fwdPos');

    sampleRate = 80;

    for i = 1:length(ftFileNames)
        % compute velocities, filter as in standard analysis pipeline
        yawAngVel = computeVelocity(yawPos(i,:));
        fwdVel = computeVelocity(fwdPos(i,:));
        slideVel = computeVelocity(slidePos(i,:));

%         % plot heatmap
        genHeatmap(yawAngVel, fwdVel, [],...
            xDataName, yDataName, zDataName, xScale, yScale, zScale, ...
            minNumVals, offsets, [], sprintf('Run %d',i));

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
        title(sprintf('Run %d',i));

    end

    % compute mean, plot velocity distribution there
    meanYawPos = mean(yawPos,1);
    meanFwdPos = mean(fwdPos,1);
    meanSlidePos = mean(slidePos,1);

    meanYawAngVel = computeVelocity(meanYawPos);
    meanFwdVel = computeVelocity(meanFwdPos);
    meanSlideVel = computeVelocity(meanSlidePos);

    % plot 3D scatterplot, shows point cloud in all 3 dimensions
    figure;
    scatter3(meanYawAngVel, ...
        meanFwdVel, ...
        meanSlideVel, ...
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
    title('Mean');
    
    genHeatmap(meanYawAngVel, meanFwdVel, [],...
        xDataName, yDataName, zDataName, xScale, yScale, zScale, ...
        minNumVals, offsets, [], 'Mean');



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