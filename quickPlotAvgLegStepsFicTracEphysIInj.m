% script for quick plots of average leg steps, FicTrac, iInj data

datPath = '/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/AnalyzedData/220913_iInj';

% 
datNames = {'220912_fly01'};
legInd = [1 2 3 4 5 6; 1 2 3 4 5 6];

durInd = 1; % 0.5 s
ampVal = 100; % amp = 1.3
phaseInd = 1; % stance

varName = 'stepLengths';
yScale = [0.1 0.6];

varName = 'stepSpeeds';
yScale = [0 5];

varName = 'stepDurations';
yScale = [0 0.5];

xScale = [-2.5 2.5];

figure;

% loop through all flies
for i = 1:length(datNames)
    fullDatPath = [datPath filesep datNames{i} '_avgLegStepsFicTracEphysIinj.mat'];

    load(fullDatPath);

    ampInd = find(amps == ampVal);

    for j = 1:6 % all legs
        subplot(2,3,j);
        thisLegInd = legInd(i,j);

        plot(legStepsIinj.durTs{durInd},...
            legStepsIinj.(varName).means{ampInd,durInd,thisLegInd,phaseInd});

        hold on;

        line([0 0],yScale,'Color','k');
        line([durs(durInd) durs(durInd)],yScale,'Color','k');
        
        ylim(yScale);
        xlim(xScale);
    end

end

%%

yScale = [0 160];
figure;

for i = 1:length(datNames)
    fullDatPath = [datPath filesep datNames{i} '_avgLegStepsFicTracEphysIinj.mat'];

    load(fullDatPath);

    ampInd = find(amps == ampVal);

    plot(ephysIinj.durTs{durInd}, ephysIinj.spikeRate.means{ampInd,durInd});

    hold on;

    line([0 0],yScale,'Color','k');
    line([durs(durInd) durs(durInd)],yScale,'Color','k');  

    ylim(yScale);
    xlim(xScale);

end





