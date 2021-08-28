% temporary script for plotting ephys data with FicTrac and visual stimuli
% 
% 4/8/21

% some constants
fictracParams.dsf = 20; % downsample to 1000 Hz;
fictracParams.filtParams.padLen = int32(200);
fictracParams.filtParams.sigmaPos = int32(100); % 100 ms
fictracParams.filtParams.sigmaVel = int32(50); % 50 ms

BALL_DIAM = 6.46;
circum = BALL_DIAM * pi; % circumference of ball, in mm
fictracParams.mmPerDeg = circum / 360; % mm per degree of ball
fictracParams.degPerMM = 360 / circum; % deg per mm ball

dataPath = '/Users/hyang/Dropbox (HMS)/EphysDataMain_RAW/210408/fly01/cell01/trial01.mat';
metadataPath = '/Users/hyang/Dropbox (HMS)/EphysDataMain_RAW/210408/fly01/cell01/metaDat.mat';

% load data
load(metadataPath);
load(dataPath);

% preprocess user DAQ
[daqData, daqOutput, daqTime] = preprocessUserDaq(inputParams, ...
    rawData, rawOutput, settings);

% preprocess ephys data
[ephysData, ephysMeta] = preprocessEphysData(daqData, daqOutput, ...
    daqTime, inputParams, settings);

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
xRange = [100 150];

figure;

subplot(4,1,1);
plot(ephysData.t, ephysData.scaledVoltage);
xlabel('Time (s)');
ylabel('mV');
xlim(xRange);

subplot(4,1,2);
plot(fictracProc.t, fictracProc.fwdVel);
xlabel('Time (s)');
ylabel('Forward Velocity (mm/s)');
xlim(xRange);

subplot(4,1,3);
plot(fictracProc.t, fictracProc.yawAngVel);
xlabel('Time (s)');
ylabel('Yaw Velocity (deg/s)');
xlim(xRange);

subplot(4,1,4);
plot(visstimT, visstimY);
xlabel('Time(s)');
ylabel('VisStim Y (V)');
xlim(xRange);

% subplot(4,1,4);
% plot(visstimT, visstimY);
% xlabel('Time(s)');
% ylabel('VisStim Y (V)');
% xlim(xRange);



