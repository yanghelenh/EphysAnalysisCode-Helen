% plotPDataRawTracesEphys.m
%
% Quick script to plot raw traces given user-selected pData file
% Written to plot traces for manuscript. 
% Adaptation from plotPDataRawTraces.m from 2P-Analysis-Code
%
% CREATED: 9/9/23 - HHY

%% plotting parameters
jPot = 13; % junction potential (to be subtracted from Vm)
Vmscale = [-70 -50];
fwdVelScale = [-5 15];
yawVelScale = [-200 200];
legScale = [-1 1];
% xScale = [62.5 64.5];
% xScale = [127.25 129.25];
xScale = [155 157];


%% load data

% curDir = pwd;
% cd(pDataPath()); % go to pData folder
% 
% % ask user to select pData file
pDataPath = '/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/pData_220622';
disp('Select a pData file to display.');
[pDataFName, pDataDirPath] = uigetfile('*pData.mat', 'Select pData file', pDataPath, ...
    'MultiSelect','off');

pDataFullPath = [pDataDirPath filesep pDataFName];
load(pDataFullPath);

%% generate plot
numSubplots = 4;

% downsample ephys data 10X
ephysVmDS = downsample(ephysData.scaledVoltage,10);
ephysTDS = downsample(ephysData.t,10);

% adjust for junction potential
ephysVmDS = ephysVmDS - jPot;

% get plotting indices for Vm
VmStartInd = find(ephysTDS < xScale(1), 1, 'last');
VmEndInd = find(ephysTDS > xScale(2), 1, 'first');

% get plotting indices for FicTrac
ftStartInd = find(fictracProc.t < xScale(1), 1, 'last');
ftEndInd = find(fictracProc.t > xScale(2), 1, 'first');

% get plotting indices for legs
legStartInd = find(legTrack.t < xScale(1), 1, 'last');
legEndInd = find(legTrack.t > xScale(2), 1, 'first');

f = figure('Position', [0 0 1200 900]);

h{1} = subplot(numSubplots, 1, 1);
plot(ephysTDS(VmStartInd:VmEndInd), ephysVmDS(VmStartInd:VmEndInd));
xlim(xScale);
ylim(Vmscale);
ylabel('Membrane potential (mV)');

h{2} = subplot(numSubplots, 1, 2);
plot(fictracProc.t(ftStartInd:ftEndInd), ...
    fictracProc.yawAngVel(ftStartInd:ftEndInd));
xlim(xScale);
ylim(yawVelScale);
ylabel('Yaw velocity (deg/s)');

h{3} = subplot(numSubplots, 1, 3);
plot(fictracProc.t(ftStartInd:ftEndInd), ...
    fictracProc.fwdVel(ftStartInd:ftEndInd));
xlim(xScale);
ylim(fwdVelScale);
ylabel('Forward velocity (mm/s)');

h{4} = subplot(numSubplots, 1, 4);
plot(legTrack.t(legStartInd:legEndInd), ...
    legTrack.srnfLegX(legStartInd:legEndInd, 1:6));
xlim(xScale);
ylim(legScale);
set(gca, 'ydir', 'reverse');
ylabel('Leg A-P Position (body lengths)');
xlabel('Time (s)');


% cd(curDir);