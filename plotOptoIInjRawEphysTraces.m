% plotOptoIInjRawEphysTraces.m
% 
% Simple script for plotting ephys scaled voltage traces during optogenetic
%  stimulation or current injection
% One row for voltage, one row for stimulation
% Option to downsample from acquisition rate for reduced size vectors
% 
% For supplementary figure of 2024 DNa02/DNg13 paper
%
% CREATED: 5/12/24 - HHY
%
% UPDATED:
%   5/12/24 - HHY
%

%% plotting parameters
vmScale = [-70 -35]; % voltage scale
stimScale = [-0.1 1.1]; % scale for stimulation
tScale = [148.8 153.8]; % time scale
tScale = [172 176]; % time scale

%% other parameters
dsRate = 5000; % downsampled aquisition rate, in Hz
junPot = -13; % junction potential correction, subtract from Vm
% pData path
pDataPath = '/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/pData_220622';

%% load data from pData file

% select pData file
[pDataFName, pDataDirPath] = uigetfile('*.mat', ...
    'Select pData file', pDataPath, 'MultiSelect', 'off');

pDataFullPath = [pDataDirPath filesep pDataFName];

% get variables saved in pData file
pDatVars = whos('-file', pDataFullPath);

pDatVarsNames = cell(size(pDatVars));

% convert pDatVars into cell array of just names
for j = 1:length(pDatVars)
    pDatVarsNames{j} = pDatVars(j).name;
end

% check if this pData file has opto or iInj
% load relevant data from pData file
if (any(strcmpi(pDatVarsNames, 'opto')))
    load(pDataFullPath, 'opto', 'ephysData');
elseif (any(strcmpi(pDatVarsNames, 'iInj')))
    load(pDataFullPath, 'iInj', 'ephysData');
end

%% process data for plotting

% remove junction potential from Vm
vm = ephysData.scaledVoltage + junPot;

% convert opto or iInj logical to plotting (0 for off, 1 for on)
if (any(strcmpi(pDatVarsNames, 'opto')))
    stimVal = zeros(size(opto.stimOnLogical));
    stimVal(opto.stimOnLogical) = 1;
elseif (any(strcmpi(pDatVarsNames, 'iInj')))
    load(pDataFullPath, 'iInj', 'ephysData');
    stimVal = zeros(size(iInj.onLogical));
    stimVal(iInj.onLogical) = 1;
end

% if downsampling
dsT = ephysData.t(1):(1/dsRate):ephysData.t(end);
dsT = dsT';

% downsampled signals
vmPlot = interp1(ephysData.t,vm,dsT);
stimPlot = interp1(ephysData.t,stimVal,dsT);

% get start and end indices, for plot times
startInd = find(dsT < tScale(1), 1, 'last');
endInd = find(dsT > tScale(2), 1, 'first');


%% plot figure

numSubplots = 2;

f = figure('Position', [0 0 1500 900]);

% plot voltage trace
subplot(numSubplots, 1, 1);
plot(dsT(startInd:endInd),vmPlot(startInd:endInd));
xlim(tScale);
ylim(vmScale);
ylabel('membrane potential (mV)');

subplot(numSubplots, 1, 2);
plot(dsT(startInd:endInd),stimPlot(startInd:endInd));
xlim(tScale);
ylim(stimScale);
ylabel('stim');

sgtitle(sprintf('%s',pDataFName));




