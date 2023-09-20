%% ANOVA tests for IInj

datDir = '/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/AnalyzedData';
whichParam = 'stepDirections';
whichPhase = 'stance';
durs = 1;
amps = [-75 100];
yScale = [-20 20];
plotAvg = 'mean';
plotIndiv = true;
plotDiff = true;

% legInds = 4:6;
% legInds = 1:3;
legInds = 1:6;

allFliesMeans = plotIInjLegStepParam_allFlies(datDir, whichParam, whichPhase,...
    durs, amps, yScale, plotAvg, plotIndiv, plotDiff);


hyperpol = squeeze(allFliesMeans(2,legInds,:))';
depol = squeeze(allFliesMeans(3,legInds,:))';
% 
diff = hyperpol - depol;
diff = diff(:);

allDat = [hyperpol(:); depol(:)];


% legG = [repmat("L1",6,1); repmat("L2",6,1); repmat("L3",6,1)];
% legG = [legG; legG];

legG = [repmat("R1",6,1); repmat("R2",6,1); repmat("R3",6,1);...
    repmat("L1",6,1); repmat("L2",6,1); repmat("L3",6,1)];
legG = [legG; legG];

% stimG = [repmat("hyperpol",18,1); repmat("depol",18,1)];
stimG = [repmat("hyperpol",36,1); repmat("depol",36,1)];

% blockG = (1:6)';
% blockG = repmat(blockG,6,1);

blockG = (1:6)';
blockG = repmat(blockG,12,1);

% p = anovan(allDat, {blockG, legG, stimG});
[p, ~, stats] = anovan(allDat, {blockG, legG, stimG}, 'model', 'interaction');
% [p, ~, stats] = anovan(allDat, {legG, stimG}, 'model', 'interaction');
% [p, ~, stats] = anova1(diff, legG);

% [c,m,h,gnames] = multcompare(stats,'Dimension', [1 2], 'CType','tukey-kramer');
% 
% [c,m,h,gnames] = multcompare(stats);


%% ANOVA tests for DNa02 opto

datDir = '/Users/hyang/Dropbox (HMS)/OptoBehaviorAnalysis-Helen/AnalyzedData';
whichParam = 'stepDirections';
whichPhase = 'stance';
durs = 0.2;
NDs = 1.3;
yScale = [-0.2 0.2];
plotAvg = 'mean';
plotIndiv = true;
plotDiff = true;


allFliesMeans = plotOptoLegStepParam_allFlies(datDir, whichParam, whichPhase,...
    durs, NDs, yScale, plotAvg, plotIndiv, plotDiff);

legInds = 1:6;

uni = squeeze(allFliesMeans(7,legInds,:));

neither = squeeze(allFliesMeans(7,legInds,:));

uni = uni';
neither = neither';

allDat = [uni(:); neither(:)];


legG1 = [repmat("R1",8,1); repmat("R2",8,1); repmat("R3",8,1); ...
    repmat("L1",8,1); repmat("L2",8,1); repmat("L3",8,1)];
legG2 = [repmat("R1",7,1); repmat("R2",7,1); repmat("R3",7,1);...
    repmat("L1",7,1); repmat("L2",7,1); repmat("L3",7,1)];
legG = [legG1; legG2];
stimG = [repmat("uni",48,1); repmat("neither",42,1)];

% legG1 = [repmat("R1",8,1); repmat("R2",8,1); repmat("R3",8,1)];
% legG2 = [repmat("R1",7,1); repmat("R2",7,1); repmat("R3",7,1)];
% legG = [legG1; legG2];
% 
% stimG = [repmat("uni",24,1); repmat("neither",21,1)];


[p,~,stats] = anovan(allDat, {legG, stimG}, 'model', 'interaction');
% [p,~,stats] = anovan(allDat, {legG, stimG});

[c,m,h,gnames] = multcompare(stats,'Dimension', [1 2], 'CType','hsd');
