% script for quick plots of average leg steps, FicTrac, opto data

datPath = '/Users/hyang/Dropbox (HMS)/OptoBehaviorAnalysis-Helen/AnalyzedData/220802_LegStepsFictracOpto_1Fly';

% unilateral
datNames = {'220415_fly01','220417_fly01','220419_fly01','220428_fly01',...
    '220517_fly02'};
legInd = [4 5 6 1 2 3; 4 5 6 1 2 3; 1 2 3 4 5 6; 1 2 3 4 5 6; 1 2 3 4 5 6];

% both
datNames = {'220331_fly03', '220426_fly01'};
legInd = [1 2 3 4 5 6; 1 2 3 4 5 6];

% neither
datNames = {'220428_fly02','220526_fly01'};
legInd = [1 2 3 4 5 6; 1 2 3 4 5 6];

durInd = 1; % dur = 0.2 s
NDInd = 3; % ND = 1.3
phaseInd = 1; % stance

varName = 'stepLengths';
yScale = [0.2 0.7];

varName = 'stepSpeeds';
yScale = [1 6];

varName = 'stepDurations';
yScale = [0 0.3];

xScale = [-5 5];

figure;

% loop through all flies
for i = 1:length(datNames)
    fullDatPath = [datPath filesep datNames{i} '_avgLegStepsFicTracOpto.mat'];

    load(fullDatPath);

    for j = 1:6 % all legs
        subplot(2,3,j);
        thisLegInd = legInd(i,j);

        plot(legStepsOpto.durTs{durInd},...
            legStepsOpto.(varName).means{NDInd,durInd,thisLegInd,phaseInd});

        hold on;

        line([0 0],yScale,'Color','k');
        line([durs(durInd) durs(durInd)],yScale,'Color','k');
        
        ylim(yScale);
        xlim(xScale);
    end

end






