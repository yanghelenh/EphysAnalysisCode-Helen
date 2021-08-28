% temporary script for plotting visual stimulus behavior data

% some constants
fictracParams.dsf = 20; % downsample to 1000 Hz;
fictracParams.filtParams.padLen = int32(200);
fictracParams.filtParams.sigmaPos = int32(100); % 100 ms
fictracParams.filtParams.sigmaVel = int32(50); % 50 ms

BALL_DIAM = 6.46;
circum = BALL_DIAM * pi; % circumference of ball, in mm
fictracParams.mmPerDeg = circum / 360; % mm per degree of ball
fictracParams.degPerMM = 360 / circum; % deg per mm ball

dataPath = '/Users/hyang/Dropbox (HMS)/VisualBehaviorData_RAW/210227/fly02/cell01/trial04.mat';
metadataPath = '/Users/hyang/Dropbox (HMS)/VisualBehaviorData_RAW/210227/fly02/cell01/metaDat.mat';

% load data
load(metadataPath);
load(dataPath);

% preprocess user DAQ
[daqData, daqOutput, daqTime] = preprocessUserDaq(inputParams, ...
    rawData, rawOutput, settings);

% preprocess FicTrac
fictrac = preprocessFicTrac(daqData, daqTime, settings.bob.sampRate);

% drop index
fictrac.dropInd = [];

% filter FicTrac
fictracProc = dsFiltFictrac(fictracParams, fictrac);

% downsample visual stimulus output
dsf = 50;

visstimT = downsample(daqTime, dsf);
visstimX = downsample(daqData.panelsDAC0X, dsf);
visstimY = downsample(daqData.panelsDAC1Y, dsf);

% plot
xRange = [150 200];

figure;

subplot(3,1,1);
plot(fictracProc.t, fictracProc.fwdVel);
xlabel('Time (s)');
ylabel('Forward Velocity (mm/s)');
xlim(xRange);

subplot(3,1,2);
plot(fictracProc.t, fictracProc.yawAngVel);
xlabel('Time (s)');
ylabel('Yaw Velocity (deg/s)');
xlim(xRange);

subplot(3,1,3);
plot(visstimT, visstimX);
xlabel('Time(s)');
ylabel('VisStim X (V)');
xlim(xRange);

% subplot(4,1,4);
% plot(visstimT, visstimY);
% xlabel('Time(s)');
% ylabel('VisStim Y (V)');
% xlim(xRange);



