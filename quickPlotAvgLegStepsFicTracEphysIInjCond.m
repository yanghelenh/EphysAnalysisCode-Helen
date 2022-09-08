% script for quick plots of average leg steps, FicTrac, iInj data

datPath = '/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/AnalyzedData/220825_iInj_cond_200msStart_200msWin_g20yawAngVel';

% 
datNames = {'220614_fly01','220621_fly01'};
legInd = [1 2 3 4 5 6; 1 2 3 4 5 6];

durInd = 2; % 0.5 s
ampVal = -50; % amp = 1.3
phaseInd = 1; % stance

varName = 'stepLengths';
yScale = [0.1 0.6];

varName = 'stepSpeeds';
yScale = [0 5];

varName = 'stepDurations';
yScale = [0 0.5];

xScale = [-5 5.5];

figure;

% loop through all flies
for i = 1:length(datNames)
    fullDatPath = [datPath filesep datNames{i} '_avgLegStepsFicTracEphysIinj_cond.mat'];

    load(fullDatPath);

    ampInd = find(amps == ampVal);

    for j = 1:6 % all legs
        subplot(2,3,j);
        thisLegInd = legInd(i,j);

        plot(legStepsCond.durTs{durInd},...
            legStepsCond.(varName).means{ampInd,durInd,thisLegInd,phaseInd});

        hold on;

        line([0 0],yScale,'Color','k');
        line([durs(durInd) durs(durInd)],yScale,'Color','k');
        
        ylim(yScale);
        xlim(xScale);
    end

end

%%

varName = 'fwdVel';
yScale = [-5,15];

varName = 'yawAngVel';
yScale = [-200, 200];

figure;

% loop through all flies
for i = 1:length(datNames)
    fullDatPath = [datPath filesep datNames{i} '_avgLegStepsFicTracEphysIinj_cond.mat'];

    load(fullDatPath);

    ampInd = find(amps == ampVal);

        plot(fictracCond.durTs{durInd},...
            fictracCond.(varName).means{ampInd,durInd});

        hold on;

        line([0 0],yScale,'Color','k');
        line([durs(durInd) durs(durInd)],yScale,'Color','k');
        
        ylim(yScale);
        xlim(xScale);

end

%%

yScale = [0 140];
figure;

for i = 1:length(datNames)
    fullDatPath = [datPath filesep datNames{i} '_avgLegStepsFicTracEphysIinj_cond.mat'];

    load(fullDatPath);

    ampInd = find(amps == ampVal);

    plot(ephysCond.durTs{durInd}, ephysCond.spikeRate.means{ampInd,durInd});

    hold on;

    line([0 0],yScale,'Color','k');
    line([durs(durInd) durs(durInd)],yScale,'Color','k');  

    ylim(yScale);
    xlim(xScale);

end





