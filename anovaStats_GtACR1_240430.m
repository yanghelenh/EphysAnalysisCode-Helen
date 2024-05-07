%% ANOVA stats for DNa02 SPARC2D GtACR1

datDir = '/Users/hyang/Dropbox (HMS)/OptoBehaviorAnalysis-Helen/AnalyzedData';
whichParam = 'stepXLengths';
whichPhase = 'stance';
vels = 0;
NDs = 0;
yScale = [-0.2 0.2];
plotAvg = 'mean';
plotIndiv = true;
plotDiff = false;

uni = plotOptomotorOptoDiffLegStepParam_allFlies(datDir, whichParam, whichPhase, vels, NDs, yScale, plotAvg, plotIndiv, plotDiff);
neither = plotOptomotorOptoDiffLegStepParam_allFlies(datDir, whichParam, whichPhase, vels, NDs, yScale, plotAvg, plotIndiv, plotDiff);
% both = plotOptomotorOptoDiffLegStepParam_allFlies(datDir, whichParam, whichPhase, vels, NDs, yScale, plotAvg, plotIndiv, plotDiff);

legInds = 1:6;
uni = squeeze(uni(5,legInds,:));
neither = squeeze(neither(5,legInds,:));
% both = squeeze(both(5,legInds,:));

uni = uni';
neither = neither';
% both = both';

allDat = [uni(:); neither(:)];
% allDat = [both(:); neither(:)];

legG1 = [repmat("R1",17,1); repmat("R2",17,1); repmat("R3",17,1); ...
    repmat("L1",17,1); repmat("L2",17,1); repmat("L3",17,1)];
% legG1 = [repmat("R1",4,1); repmat("R2",4,1); repmat("R3",4,1); ...
%     repmat("L1",4,1); repmat("L2",4,1); repmat("L3",4,1)];

legG2 = [repmat("R1",7,1); repmat("R2",7,1); repmat("R3",7,1);...
    repmat("L1",7,1); repmat("L2",7,1); repmat("L3",7,1)];
legG = [legG1; legG2];

stimG = [repmat("uni",102,1); repmat("neither",42,1)];
% stimG = [repmat("both",24,1); repmat("neither",42,1)];

[p,~,stats] = anovan(allDat, {legG, stimG}, 'model', 'interaction');

figure;
[c,m,h,gnames] = multcompare(stats,'Dimension', [1 2], 'CType','tukey-kramer');

%% label legs based only on side

datDir = '/Users/hyang/Dropbox (HMS)/OptoBehaviorAnalysis-Helen/AnalyzedData';
whichParam = 'stepPEPX';
whichPhase = 'stance';
vels = 0;
NDs = 0;
yScale = [-0.1 0.1];
plotAvg = 'mean';
plotIndiv = true;
plotDiff = false;

uni = plotOptomotorOptoDiffLegStepParam_allFlies(datDir, whichParam, ...
    whichPhase, vels, NDs, yScale, plotAvg, plotIndiv, plotDiff);
neither = plotOptomotorOptoDiffLegStepParam_allFlies(datDir, whichParam,...
    whichPhase, vels, NDs, yScale, plotAvg, plotIndiv, plotDiff);

legInds = 1:6;
uni = squeeze(uni(5,legInds,:));
neither = squeeze(neither(5,legInds,:));

uni = uni';
neither = neither';

allDat = [uni(:); neither(:)];

legG1 = [repmat("I",17,1); repmat("I",17,1); repmat("I",17,1); ...
    repmat("C",17,1); repmat("C",17,1); repmat("C",17,1)];
legG2 = [repmat("I",7,1); repmat("I",7,1); repmat("I",7,1);...
    repmat("C",7,1); repmat("C",7,1); repmat("C",7,1)];
legG = [legG1; legG2];
stimG = [repmat("uni",102,1); repmat("neither",42,1)];

[p,~,stats] = anovan(allDat, {legG, stimG}, 'model', 'interaction');

figure;
[c,m,h,gnames] = multcompare(stats,'Dimension', [1 2], 'CType','tukey-kramer');

%% ANOVA tests for DNa02 opto - circular

datDir = '/Users/hyang/Dropbox (HMS)/OptoBehaviorAnalysis-Helen/AnalyzedData';
whichParam = 'stepDirections';
whichPhase = 'stance';
vels = 0;
NDs = 0;
yScale = [-30 30];
plotAvg = 'mean';
plotIndiv = true;
plotDiff = false;

% uni = plotOptomotorOptoDiffLegStepParam_allFlies(datDir, whichParam, whichPhase, vels, NDs, yScale, plotAvg, plotIndiv, plotDiff);
neither = plotOptomotorOptoDiffLegStepParam_allFlies(datDir, whichParam, whichPhase, vels, NDs, yScale, plotAvg, plotIndiv, plotDiff);
both = plotOptomotorOptoDiffLegStepParam_allFlies(datDir, whichParam, whichPhase, vels, NDs, yScale, plotAvg, plotIndiv, plotDiff);

legInds = 1:6;
% uni = squeeze(uni(5,legInds,:));
neither = squeeze(neither(5,legInds,:));
both = squeeze(both(5,legInds,:));

% uni = uni';
neither = neither';
both = both';

% allDat = [uni(:); neither(:)];
allDat = [both(:); neither(:)];

% legG1 = [repmat("R1",17,1); repmat("R2",17,1); repmat("R3",17,1); ...
%     repmat("L1",17,1); repmat("L2",17,1); repmat("L3",17,1)];
legG1 = [repmat("R1",4,1); repmat("R2",4,1); repmat("R3",4,1); ...
    repmat("L1",4,1); repmat("L2",4,1); repmat("L3",4,1)];

legG2 = [repmat("R1",7,1); repmat("R2",7,1); repmat("R3",7,1);...
    repmat("L1",7,1); repmat("L2",7,1); repmat("L3",7,1)];
legG = [legG1; legG2];

% stimG = [repmat("uni",102,1); repmat("neither",42,1)];
stimG = [repmat("both",24,1); repmat("neither",42,1)];


[p, stats] = circ_hktest(deg2rad(allDat), legG, stimG, 1, {'Leg','Stim'});

% post-hoc tests
whichLeg = 6;
[p, stats] = circ_wwtest(deg2rad(uni(:,whichLeg)), deg2rad(neither(:,whichLeg)));