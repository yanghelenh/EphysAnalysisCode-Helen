% plotSpikePhaseLockExTrace.m
%
% Quick script to plot leg phase and ephys trace

%% plotting parameters
jPot = 13; % junction potential (to be subtracted from Vm)
Vmscale = [-70 -50];
Vmscale = [-65 -45];
phaseScale = [0 360];
xScale = [180 181]; % DNa02  
xScale = [176.58 177.58]; % DNg13
xScale = [176.15 177.15]; % DNg13 2
xScale = [163.2 164.2]; % DNg13 3

% pData file for DNa02
pDataFName = '211115_fly01_cell01_trial12_pData.mat';
% pData file for DNg13
pDataFName = '220620_fly02_cell01_trial01_pData.mat';

% pData path
pDataPath = '/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/pData_220622';

%% prompt user to select pData files
    if isempty(pDataFName)
        [pDataFName, pDataDirPath] = uigetfile('*.mat', ...
            'Select pData file', pDataPath, 'MultiSelect', 'off');
    else
        pDataDirPath = pDataPath;
    end

            pDataName = pDataFName;

        pDataFullPath = [pDataDirPath filesep pDataName];

        % get variables saved in pData file
        pDatVars = whos('-file', pDataFullPath);
    
        pDatVarsNames = cell(size(pDatVars));
        
        % convert pDatVars into cell array of just names
        for j = 1:length(pDatVars)
            pDatVarsNames{j} = pDatVars(j).name;
        end


% load variables from pData
        if (any(strcmpi(pDatVarsNames, 'iInj')))
            load(pDataFullPath, 'legPhase', 'fictracProc', ...
                'legTrack', 'legStepsCont', 'ephysSpikes', 'legSteps', ...
                'iInj', 'ephysData');
        else
            load(pDataFullPath, 'legPhase', 'fictracProc', ...
                'legTrack', 'legStepsCont', 'ephysSpikes', 'legSteps', ...
                'ephysData');
        end
        
%% get phase
        refLegInd = 5; % DNg13
        refLegInd = 2; % DNa02

        % get all mean offsets
        medOffPhase = zeros(size(legPhase));
        for k = 1:6
            thisOffset = deg2rad(legPhase(:,k) - legPhase(:,refLegInd));
            medOff = circ_mean(thisOffset(~isnan(thisOffset)));
            medOffPhase(:,k) = deg2rad(legPhase(:,k)) - medOff;
        end

                phaseVals = nan(size(legPhase,1),1);
        % get mean phase estimate
        for k = 1:size(medOffPhase,1)
            thisRow = medOffPhase(k,:);
            thisRow(isnan(thisRow)) = [];
            if ~isempty(thisRow)
                phaseVals(k) = wrapTo2Pi(circ_mean(thisRow'));
            end
        end

        phaseValsDeg = rad2deg(phaseVals);

%% generate plot
numSubplots = 2;

% downsample ephys data 10X
ephysVmDS = downsample(ephysData.scaledVoltage,10);
ephysTDS = downsample(ephysData.t,10);

% adjust for junction potential
ephysVmDS = ephysVmDS - jPot;

% get plotting indices for Vm
VmStartInd = find(ephysTDS < xScale(1), 1, 'last');
VmEndInd = find(ephysTDS > xScale(2), 1, 'first');

% get plotting indices for legs
legStartInd = find(legTrack.t < xScale(1), 1, 'last');
legEndInd = find(legTrack.t > xScale(2), 1, 'first');

f = figure('Position', [0 0 400 900]);

h{1} = subplot(numSubplots, 1, 1);
plot(ephysTDS(VmStartInd:VmEndInd), ephysVmDS(VmStartInd:VmEndInd));
xlim(xScale);
ylim(Vmscale);
ylabel('Membrane potential (mV)');

h{2} = subplot(numSubplots, 1, 2);
plot(legTrack.t(legStartInd:legEndInd), ...
    phaseValsDeg(legStartInd:legEndInd));
xlim(xScale);
ylim(phaseScale);
ylabel('step phase (deg)');
xlabel('Time (s)');