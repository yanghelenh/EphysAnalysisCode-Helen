% createEphysLegFictracMovie.m
%
% Quick and dirty script to generate movie showing ephys recording, FicTrac
%  variables, and leg movie. Starts from raw data.
% Made 8/7/20, used for videos for leading edge symposium and U19 talk. Now
%  obsolete. Make new functions if new videos desired.

% full path to data folder
dataPath = '/Users/hyang/Dropbox (HMS)/EphysData_RAW/200826/fly01/cell01';
% trial name
trialName = 'trial01.mat';
% full path to leg video
vidPath = '/Users/hyang/Desktop';
vidName = 'fly01_cell01_trial01_legVid';

tVid = [50 60]; % start and end times, in seconds

cd(dataPath);

% load data
load(trialName);

% load ephys settings
[dataDir, exptFnDir, settings] = ephysSettings();

% preprocess DAQ data
[daqData, daqOutput, daqTime] = preprocessUserDaq(inputParams, ...
    rawData, rawOutput, settings);

% preprocess ephys data
[ephysData, ephysMeta] = preprocessEphysData(daqData, daqOutput, ...
    daqTime, inputParams, settings);
% preprocess leg video trigger data
[legVidFrameTimes, legVidTrigTimes] = preprocessLegVid(...
    daqData, daqOutput, daqTime);
% preprocess FicTrac data
[yawAngVel, yawAngPosWrap, fwdVel, fwdCumPos, slideVel, slideCumPos, ...
    xPos, yPos] = preprocessFicTrac(daqData, settings.bob.sampRate);

ficTracTimes = daqTime;

% FicTrac variables
sampRate = 1/(median(diff(ficTracTimes)));
avgWindow = 0.2;
smoYawAngVel = moveAvgFilt(yawAngVel, sampRate, avgWindow);
smoFwdVel = moveAvgFilt(fwdVel, sampRate, avgWindow);


% legVid variables
legStartInd = find(legVidFrameTimes >= tVid(1), 1, 'first');
legEndInd = find(legVidFrameTimes < tVid(2), 1, 'last');
legVidPath = [vidPath filesep vidName];
legVidImgSize = [448 448];

legVidImg = zeros(legVidImgSize(1),legVidImgSize(2), ...
    legEndInd - legStartInd + 1);
for i = legStartInd:1:legEndInd
    legVidFile = sprintf('legVid-%i.tiff', i-1+inputParams.startLegVidNum);
    legVidImg(:,:,i - legStartInd + 1) = imread(...
        [legVidPath filesep legVidFile]);
end
legVidMinInt = 30;
legVidMaxInt = 160;

legVidFrameRate = 1/(median(diff(legVidFrameTimes)));

% get indicies to display on each frame
ficTracDispInd = zeros(1,legEndInd - legStartInd + 1);
for i = 1:(legEndInd - legStartInd + 1)
    ficTracDispInd(i) = find(ficTracTimes == ...
        legVidFrameTimes(i + legStartInd - 1), 1, 'first');    
end
ephysDispInd = ficTracDispInd;


ficTracStartInd = ficTracDispInd(1);
ficTracFigTimes = ficTracTimes(ficTracDispInd) - ficTracTimes(ficTracStartInd);
smoFwdVelFig = smoFwdVel(ficTracDispInd);
smoAngVelFig = smoYawAngVel(ficTracDispInd);
ephysVFig = ephysData.scaledVoltage;

fwdVelScale = [-5 10];
angVelScale = [-200 200];
ephysVScale = [-55 -25];

timeScale = [0 10];

% make figure
figFrames(legEndInd - legStartInd + 1) = struct('cdata',[],'colormap',[]);
% figFrames(200) = struct('cdata',[],'colormap',[]);
f = figure;
set(f,'Position',[10,10,1000,700]);
for i = 1:(legEndInd - legStartInd + 1)
    ficTracInd = ficTracDispInd(i);
        
    % leg video
    legAx = subplot(3, 3, [1, 4, 7]);
    imagesc(legVidImg(:,:,i),[legVidMinInt legVidMaxInt]);
    axis equal;
    axis tight;
    colormap(legAx,'gray');
    set(gca,'XTick',[],'YTick',[]);
    legPos = get(legAx,'Position');
    newLegPos = [legPos(1)-0.06, legPos(2)-0.2, legPos(3)*1.5, legPos(4)*1.5];
    set(legAx,'Position',newLegPos);
    
    % ephys trace
    ephysTime = daqTime(ephysDispInd(1):ephysDispInd(i)) - daqTime(ephysDispInd(1));
    ephysAx = subplot(3, 3, [2,3]);
    plot(ephysTime, ephysVFig(ephysDispInd(1):ephysDispInd(i)), 'b');
    ylim(ephysVScale);
    xlim(timeScale);
    yticks([-55 -45 -35 -25]);
    xlabel('Time (s)');
    ylabel('Membrane Potential (mV)');
    ephysPos = get(ephysAx,'Position');
    newEphysPos = [ephysPos(1)+0.06, ephysPos(2), ephysPos(3), ephysPos(4)*0.9];
    set(ephysAx,'Position',newEphysPos);
    title('Electrophysiology Recording');
        
    
    % FicTrac forward velocity
    fwdVelAx = subplot(3,3,[5,6]);
    plot(ficTracFigTimes(1:i),...
        smoFwdVelFig(1:i), 'k',...
        'LineWidth',1.5);
    ylim(fwdVelScale);
    xlim(timeScale);
    xlabel('Time (s)');
    ylabel('mm/s');
    fwdVelPos = get(fwdVelAx,'Position');
    newFVPos = [fwdVelPos(1)+0.06, fwdVelPos(2), fwdVelPos(3), fwdVelPos(4)*0.9];
    set(fwdVelAx,'Position',newFVPos);   
    title('Forward Velocity');
    
    % FicTrac yaw velocity
    yawAx = subplot(3,3,[8,9]);
    plot(ficTracFigTimes(1:i),...
        smoAngVelFig(1:i), 'k',...
        'LineWidth',1.5);
    ylim(angVelScale);
    xlim(timeScale);
    xlabel('Time (s)');
    ylabel('deg/s');
    yawPos = get(yawAx,'Position');
    newYawPos = [yawPos(1)+0.06, yawPos(2), yawPos(3), yawPos(4)*0.9];
    set(yawAx,'Position',newYawPos);   
    title('Rotational Velocity');   
    
    drawnow;
    figFrames(i) = getframe(gcf); 
end

% save video
vidName = [pwd filesep 'figureVid.mp4'];

v = VideoWriter(vidName, 'MPEG-4');
v.FrameRate = legVidFrameRate;
open(v);
writeVideo(v, figFrames);
close(v);
