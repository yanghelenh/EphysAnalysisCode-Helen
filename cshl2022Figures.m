% cshl2022Figures.m
%
% Quick and dirty script for making ephys figures for DNg13 and DNa02 for 
%  2022 CSHL Neuronal Circuits poster
%


%% DNg13
pDataFullPath = '/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/pData_2021cshl/200821_fly01_cell01_trial01.mat';
ttl = 'DNg13';

%% DNa02
pDataFullPath = '/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/pData_220622/211104_fly01_cell01_trial01_pData.mat';
% pDataFullPath = '/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/pData_2021cshl/201228_fly01_cell01_trial03_pData.mat';
ttl = 'DNa02';

%% load data
load(pDataFullPath);

%% t delay = -50 ms
tDelay = -0.05;


%% get FicTrac params during moving bouts only
bufferT = 0.1; % time in sec around bout start/end to remove

% bout start and end times
moveBoutStartT = legTrack.t(moveNotMove.moveBout(:,1)) + bufferT;
moveBoutEndT = legTrack.t(moveNotMove.moveBout(:,2)) - bufferT;

% convert times into FicTrac indices
fictracMoveInd = [];

for i = 1:length(moveBoutStartT)
    validInd = find((fictracProc.t >= moveBoutStartT(i)) & ...
        (fictracProc.t <= moveBoutEndT(i)));
    fictracMoveInd = [fictracMoveInd; validInd];
end

% get Fictrac params at these ind
fwdVel = fictracProc.fwdVel(fictracMoveInd)';
yawAngVel = fictracProc.yawAngVel(fictracMoveInd)';

% spike rate at FicTrac ind, with tDelay offset
spikeRate = interp1(ephysSpikes.t+tDelay, ephysSpikes.spikeRate, ...
    fictracProc.t);
spikeRate = spikeRate(fictracMoveInd);

%%
% get Fictrac params at these ind
fwdVel = fictracProc.fwdVel';
yawAngVel = fictracProc.yawAngVel';

% spike rate at FicTrac ind, with tDelay offset
spikeRate = interp1(ephysSpikes.t+tDelay, ephysSpikes.spikeRate, ...
    fictracProc.t);

%%
% get Fictrac params at these ind
fwdVel = fictrac.fwdVel;
yawAngVel = fictrac.yawAngVel;

% spike rate at FicTrac ind, with tDelay offset
spikeRate = interp1(ephysSpikes.t+tDelay, ephysSpikes.spikeRate, ...
    fictrac.t(1:end-4));




%% heatmaps

xDataName = 'yawAngVel';
yDataName = 'fwdVel';
zDataName = 'spikeRate';

[f, heatmapMat, countsMat] = genHeatmap(yawAngVel, fwdVel, spikeRate,...
    xDataName, yDataName, zDataName, [-500 500 30], [-5 15 30], [0 70], 20,...
    0, [], ttl);


%% kernels

% get Fictrac
fwdVel = fictracProc.fwdVel';
yawAngVel = fictracProc.yawAngVel';

% spike rate at FicTrac ind, with tDelay offset
spikeRate = interp1(ephysSpikes.t+tDelay, ephysSpikes.spikeRate, ...
    fictracProc.t);

kernelParams.winLen = 0.5;
kernelParams.cutFreq = 0;
kernelParams.tauFreq = 0;
kernelParams.fwdKernelBW = 100;
kernelParams.revKernelBW = 100;
kernelParams.sampRate = 1000;

% compute forward kernel
[tempKern, lags, numSeg] = computeWienerKernel(...
    yawAngVel, ...
    spikeRate, ...
    kernelParams.sampRate, ...
    kernelParams.winLen,...
    kernelParams.cutFreq, ...
    kernelParams.tauFreq);

fwdKernel = ...                            
    slepianWinFilter(tempKern, kernelParams.fwdKernelBW, ...
    kernelParams.sampRate);

% compute reverse kernel
[tempKern, ~, numSeg] = computeWienerKernel(...
    spikeRate, ...
    yawAngVel, ...
    kernelParams.sampRate, ...
    kernelParams.winLen,...
    kernelParams.cutFreq, ...
    kernelParams.tauFreq);

revKernel = ...
    slepianWinFilter(tempKern, kernelParams.revKernelBW, ...
    kernelParams.sampRate);


%% stance speed
legIDs.ind = 1:6;

% for each half step

% delays in sec, neg delay is ephys b/f behavior
tDelay = [-0.5 -0.2 -0.1 -0.05 -0.025 0 0.025 0.05 0.1 0.2 0.5]; 

% preallocate for spike rate matrix for (number of steps x 2 (for 2 half
%  steps) x num delays
stepSpikeRate = zeros(size(legSteps.stepInds,1),2, length(tDelay));

% times when all spikes occured
spikeTimes = ephysSpikes.t(ephysSpikes.startInd);


% loop through all time delays
for j = 1:length(tDelay)
    % incorporate time offset b/w ephys and behavior
    spikeTimesDelay = spikeTimes - tDelay(j);

    % loop through all steps
    for i = 1:size(legSteps.stepInds,1)
        % get step start, mid, end times
        stepStartT = legTrack.t(legSteps.stepInds(i,1));
        stepMidT = legTrack.t(legSteps.stepInds(i,2));
        stepEndT = legTrack.t(legSteps.stepInds(i,3));

        % for first half step, figure out how many spikes 
        numSpikesStartMid = sum((spikeTimesDelay >= stepStartT) &...
            (spikeTimesDelay < stepMidT));
        % convert to step rate
        stepSpikeRate(i,1,j) = numSpikesStartMid / (stepMidT-stepStartT);

        % for second half step, figure out how many spikes 
        numSpikesMidEnd = sum((spikeTimesDelay >= stepMidT) &...
            (spikeTimesDelay < stepEndT));
        % convert to step rate
        stepSpikeRate(i,2,j) = numSpikesMidEnd / (stepEndT-stepMidT); 
    end
end

% which parameter, comment out as needed

% thisStepParam = legSteps.stepSpeeds(:,2);
whichParamStr = 'Stance Step Speed (body lengths/s)';
thisStepParam = legSteps.stepLengths(:,2);
whichParamStr = 'Step Lengths (body lengths)';
thisEphysSpikes = stepSpikeRate;
thisWhichLeg = legSteps.stepWhichLeg;

xLimits = [0 0.75];
yLimits = [0 150];
% colormap
colMap = lines(6);

% loop through all tDelay
for j = 4
    
    % spikes just for this delay
    thisDelayEphysSpikes = thisEphysSpikes(:,j);
    
    figure('Position', [10 10 1000 400]);
    
    % ipsi legs (right)
    subplot(1,2,1);
    
    for i = 1:3
        thisLeg = legIDs.ind(i);
        % steps for this leg
        thisLegStepParams = thisStepParam(thisWhichLeg == thisLeg);
        % ephys for steps for this leg
        thisLegStepSpikes = thisDelayEphysSpikes(thisWhichLeg == thisLeg);
        
        % remove outliers
        [~,stepParamOutInd] = rmoutliers(thisLegStepParams);
        stepParamNoOut = thisLegStepParams(~stepParamOutInd);
        spikeNoOut = thisLegStepSpikes(~stepParamOutInd);
        
        % fit line
        p = polyfit(stepParamNoOut, spikeNoOut, 1);
        f = polyval(p,linspace(xLimits(1),xLimits(2)));

        scatter(stepParamNoOut,spikeNoOut, 80, '.', ...
            'MarkerEdgeColor', colMap(i,:), 'MarkerFaceColor', colMap(i,:));
        hold on;
        plot(linspace(xLimits(1),xLimits(2)),f, 'Color', colMap(i,:), ...
            'LineWidth', 2);
        legend;
    end
    ylim(yLimits);
    xlim(xLimits);
    
    % contra legs (left)
    subplot(1,2,2);
    for i = 4:6
        thisLeg = legIDs.ind(i);
        % steps for this leg
        thisLegStepParams = thisStepParam(thisWhichLeg == thisLeg);
        % ephys for steps for this leg
        thisLegStepSpikes = thisDelayEphysSpikes(thisWhichLeg == thisLeg);
        
        % remove outliers
        [~,stepParamOutInd] = rmoutliers(thisLegStepParams);
        stepParamNoOut = thisLegStepParams(~stepParamOutInd);
        spikeNoOut = thisLegStepSpikes(~stepParamOutInd);
        
        % fit line
        p = polyfit(stepParamNoOut, spikeNoOut, 1);
        f = polyval(p,linspace(xLimits(1),xLimits(2)));

        scatter(thisLegStepParams,thisLegStepSpikes, 80, '.',...
            'MarkerEdgeColor', colMap(i,:), 'MarkerFaceColor', colMap(i,:));
        hold on;
        plot(linspace(xLimits(1),xLimits(2)),f, 'Color', colMap(i,:), ...
            'LineWidth', 2);  
        legend;
    end
    ylim(yLimits);
    xlim(xLimits);

end
