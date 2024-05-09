% anovaStats_optoIInj_240509.m
%
% Script to perform ANOVA statistical tests and post-hoc tests for
%  manipulation experiments (CsCh and GtACR1 opto and IInj)
%
% CREATED: 5/9/24 - HHY
%
% UPDATED: 5/9/24 - HHY
%

%% DNa02>>SPARC2D-CsCh, non-circular
% as difference from no stim condition

% which parameter
% whichParam = 'stepXLengthsAbs';
% whichParam = 'stepAEPX';
whichParam = 'stepPEPX';

% input parameters
datDir = '/Users/hyang/Dropbox (HMS)/OptoBehaviorAnalysis-Helen/AnalyzedData';
whichPhase = 'stance';
durs = 0.2;
NDs = 1.3;
yScale = [-0.2 0.2];
plotAvg = 'mean';
plotIndiv = true;
plotDiff = true;

% leg info
legInds = 1:6;

% which index of output matrices
ndDurInd = 7;

% get data points
% select unilateral
uni = plotOptoLegStepParam_allFlies(datDir, whichParam, ...
    whichPhase, durs, NDs, yScale, plotAvg, plotIndiv, plotDiff);
% select neither
neither = plotOptoLegStepParam_allFlies(datDir, whichParam, ...
    whichPhase, durs, NDs, yScale, plotAvg, plotIndiv, plotDiff);
% select both
both = plotOptoLegStepParam_allFlies(datDir, whichParam, ...
    whichPhase, durs, NDs, yScale, plotAvg, plotIndiv, plotDiff);

% close all figures
close all

% get only relevant ND, duration combo
uni = squeeze(uni(ndDurInd,legInds,:));
neither = squeeze(neither(ndDurInd,legInds,:));
both = squeeze(both(ndDurInd,legInds,:));

% build tables for ranova function
% measurements - uni and neither
uniNeiCond = ...
    [repmat("uni",size(uni,2),1); repmat("neither", size(neither,2),1)];
uniNeiData = cat(2,uni,neither)';

uniNeiTbl = table(uniNeiCond,uniNeiData(:,1),uniNeiData(:,2),...
    uniNeiData(:,3), uniNeiData(:,4), uniNeiData(:,5), uniNeiData(:,6), ...
    'VariableNames',{'Opto','y1', 'y2', 'y3', 'y4', 'y5', 'y6'});
% within
uniNeiWithin = table([repmat("ipsi",3,1); repmat("contra",3,1)], ...
    ["front"; "mid"; "hind"; "front"; "mid"; "hind"], 'VariableNames', ...
    {'side', 'leg'});
uniNeiWithin2 = table(["R1";"R2";"R3";"L1";"L2";"L3"], 'VariableNames',...
    {'leg'});

% ranova
rm = fitrm(uniNeiTbl,'y1-y6~Opto','WithinDesign',uniNeiWithin);
uniNeiRanova = ranova(rm,'WithinModel','side+leg');

rm2 = fitrm(uniNeiTbl,'y1-y6~Opto','WithinDesign',uniNeiWithin2);
uniNeiRanova2 = ranova(rm2, "WithinModel",'leg');

% post-hoc tests
uniNeiPostHoc = multcompare(rm,'Opto','By','side',...
    'ComparisonType','tukey-kramer');
uniNeiPostHoc2 = multcompare(rm2,'Opto','By','leg',...
    'ComparisonType','tukey-kramer');



% build tables for ranova function
% measurements - both and neither
% condition labels
bothNeiCond = ...
    [repmat("both",size(both,2),1); repmat("neither", size(neither,2),1)];
% combined data
bothNeiData = cat(2,both,neither)';

% table
bothNeiTbl = table(bothNeiCond,bothNeiData(:,1),bothNeiData(:,2),...
    bothNeiData(:,3), bothNeiData(:,4), bothNeiData(:,5), bothNeiData(:,6), ...
    'VariableNames',{'Opto','y1', 'y2', 'y3', 'y4', 'y5', 'y6'});

% within
bothNeiWithin = table([repmat("ipsi",3,1); repmat("contra",3,1)], ...
    ["front"; "mid"; "hind"; "front"; "mid"; "hind"], 'VariableNames', ...
    {'side', 'leg'});
bothNeiWithin2 = table(["R1";"R2";"R3";"L1";"L2";"L3"], 'VariableNames',...
    {'leg'});


% ranova
rm = fitrm(bothNeiTbl,'y1-y6~Opto','WithinDesign',bothNeiWithin);
bothNeiRanova = ranova(rm,'WithinModel','side+leg');

rm2 = fitrm(bothNeiTbl,'y1-y6~Opto','WithinDesign',bothNeiWithin2);
bothNeiRanova2 = ranova(rm2, "WithinModel",'leg');


% post-hoc tests
bothNeiPostHoc = multcompare(rm,'Opto','By','side',...
    'ComparisonType','tukey-kramer');
bothNeiPostHoc2 = multcompare(rm2,'Opto','By','leg',...
    'ComparisonType','tukey-kramer');


%% DNa02>>SPARC2D-GtACR1, non-circular
% static - gray

% which parameter
whichParam = 'stepXLengthsAbs';
whichParam = 'stepAEPX';
whichParam = 'stepPEPX';

% input parameters
datDir = '/Users/hyang/Dropbox (HMS)/OptoBehaviorAnalysis-Helen/AnalyzedData';
whichPhase = 'stance';
vels = 0;
NDs = 0;
yScale = [-0.2 0.2];
plotAvg = 'mean';
plotIndiv = true;
plotDiff = false;

% leg info
legInds = 1:6;


% index for particular ND
thisNDInd = 5;

% get data points
% select unilateral, 1. gray, 2. static
uni = plotOptomotorOptoDiffLegStepParam_allFlies(...
    datDir, whichParam, whichPhase, vels, NDs, yScale, plotAvg, ...
    plotIndiv, plotDiff);
% select neither, 1. gray, 2. static
neither = plotOptomotorOptoDiffLegStepParam_allFlies(...
    datDir, whichParam, whichPhase, vels, NDs, yScale, plotAvg, ...
    plotIndiv, plotDiff);
% select both, 1. gray, 2. static
both = plotOptomotorOptoDiffLegStepParam_allFlies(...
    datDir, whichParam, whichPhase, vels, NDs, yScale, plotAvg, ...
    plotIndiv, plotDiff);

% close all figures
close all

% get only relevant ND
uni = squeeze(uni(thisNDInd,legInds,:));
neither = squeeze(neither(thisNDInd,legInds,:));
both = squeeze(both(thisNDInd,legInds,:));


% build tables for ranova function
% measurements - uni and neither
% condition labels
uniNeiCond = ...
    [repmat("uni",size(uni,2),1); repmat("neither", size(neither,2),1)];
% combined data
uniNeiData = cat(2,uni,neither)';

% table
uniNeiTbl = table(uniNeiCond,uniNeiData(:,1),uniNeiData(:,2),...
    uniNeiData(:,3), uniNeiData(:,4), uniNeiData(:,5), uniNeiData(:,6), ...
    'VariableNames',{'Opto','y1', 'y2', 'y3', 'y4', 'y5', 'y6'});

% within
uniNeiWithin = table([repmat("ipsi",3,1); repmat("contra",3,1)], ...
    ["front"; "mid"; "hind"; "front"; "mid"; "hind"], 'VariableNames', ...
    {'side', 'leg'});
uniNeiWithin2 = table(["R1";"R2";"R3";"L1";"L2";"L3"], 'VariableNames',...
    {'leg'});


% ranova
rm = fitrm(uniNeiTbl,'y1-y6~Opto','WithinDesign',uniNeiWithin);
uniNeiRanova = ranova(rm,'WithinModel','side+leg');

rm2 = fitrm(uniNeiTbl,'y1-y6~Opto','WithinDesign',uniNeiWithin2);
uniNeiRanova2 = ranova(rm2, "WithinModel",'leg');


% post-hoc tests
uniNeiPostHoc = multcompare(rm,'Opto','By','side',...
    'ComparisonType','tukey-kramer');
uniNeiPostHoc2 = multcompare(rm2,'Opto','By','leg',...
    'ComparisonType','tukey-kramer');



% build tables for ranova function
% measurements - both and neither
% condition labels
bothNeiCond = ...
    [repmat("both",size(both,2),1); repmat("neither", size(neither,2),1)];
% combined data
bothNeiData = cat(2,both,neither)';

% table
bothNeiTbl = table(bothNeiCond,bothNeiData(:,1),bothNeiData(:,2),...
    bothNeiData(:,3), bothNeiData(:,4), bothNeiData(:,5), bothNeiData(:,6), ...
    'VariableNames',{'Opto','y1', 'y2', 'y3', 'y4', 'y5', 'y6'});

% within
bothNeiWithin = table([repmat("ipsi",3,1); repmat("contra",3,1)], ...
    ["front"; "mid"; "hind"; "front"; "mid"; "hind"], 'VariableNames', ...
    {'side', 'leg'});
bothNeiWithin2 = table(["R1";"R2";"R3";"L1";"L2";"L3"], 'VariableNames',...
    {'leg'});


% ranova
rm = fitrm(bothNeiTbl,'y1-y6~Opto','WithinDesign',bothNeiWithin);
bothNeiRanova = ranova(rm,'WithinModel','side+leg');

rm2 = fitrm(bothNeiTbl,'y1-y6~Opto','WithinDesign',bothNeiWithin2);
bothNeiRanova2 = ranova(rm2, "WithinModel",'leg');


% post-hoc tests
bothNeiPostHoc = multcompare(rm,'Opto','By','side',...
    'ComparisonType','tukey-kramer');
bothNeiPostHoc2 = multcompare(rm2,'Opto','By','leg',...
    'ComparisonType','tukey-kramer');



%% DNg13 IInj, non-circular

% which parameter
% whichParam = 'stepXLengthsAbs';
% whichParam = 'stepAEPX';
whichParam = 'stepPEPX';

% input parameters
datDir = '/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/AnalyzedData';
whichPhase = 'stance';
durs = 1;
amps = [-75 100];
yScale = [-0.1 0.1];
plotAvg = 'mean';
plotIndiv = true;
plotDiff = true;

% leg info
legInds = 1:6;

% get data points
allFliesMeans = plotIInjLegStepParam_allFlies(datDir, whichParam, whichPhase,...
    durs, amps, yScale, plotAvg, plotIndiv, plotDiff);

% close figures
close all

% get hyperpol and depol separately
hyperpol = squeeze(allFliesMeans(2,legInds,:))';
depol = squeeze(allFliesMeans(3,legInds,:))';

% label on conditions
condLbl = [repmat("hyperpol",size(hyperpol,2),1); ...
    repmat("depol", size(depol,2),1)];
% combined data
cmbData = cat(1,hyperpol,depol);
% table
dng13Tbl = table(condLbl,cmbData(:,1),cmbData(:,2),cmbData(:,3),...
    cmbData(:,4), cmbData(:,5), cmbData(:,6), 'VariableNames', {'iInj', ...
    'y1','y2','y3','y4','y5','y6'});

% within
dng13Within1 = table([repmat("ipsi",3,1); repmat("contra",3,1)], ...
    ["front"; "mid"; "hind"; "front"; "mid"; "hind"], 'VariableNames', ...
    {'side', 'leg'});
dng13Within2 = table(["R1";"R2";"R3";"L1";"L2";"L3"], 'VariableNames',...
    {'leg'});


% ranova
rm = fitrm(dng13Tbl,'y1-y6~iInj','WithinDesign',dng13Within1);
dng13Ranova1 = ranova(rm,'WithinModel','side+leg');

rm2 = fitrm(dng13Tbl,'y1-y6~iInj','WithinDesign',dng13Within2);
dng13Ranova2 = ranova(rm2, "WithinModel",'leg');


% post-hoc tests
dng13PostHoc1 = multcompare(rm,'iInj','By','side',...
    'ComparisonType','tukey-kramer');
dng13PostHoc2 = multcompare(rm2,'iInj','By','leg',...
    'ComparisonType','tukey-kramer');



%% DNa02>>SPARC2D-CsCh, circular

datDir = '/Users/hyang/Dropbox (HMS)/OptoBehaviorAnalysis-Helen/AnalyzedData';
whichParam = 'stepDirections';
whichPhase = 'stance';
durs = 0.2;
NDs = 1.3;
yScale = [-25 25];
plotAvg = 'mean';
plotIndiv = true;
plotDiff = true;

% leg indices
legInds = 1:6;
numLegs = length(legInds);

% which index of output matrices
ndDurInd = 7;

% get all fly means
uni = plotOptoLegStepParam_allFlies(datDir, whichParam, whichPhase,...
    durs, NDs, yScale, plotAvg, plotIndiv, plotDiff);

neither = plotOptoLegStepParam_allFlies(datDir, whichParam, whichPhase,...
    durs, NDs, yScale, plotAvg, plotIndiv, plotDiff);

both = plotOptoLegStepParam_allFlies(datDir, whichParam, whichPhase,...
    durs, NDs, yScale, plotAvg, plotIndiv, plotDiff);

numFliesUni = size(uni,3);
numFliesNei = size(neither,3);
numFliesBoth = size(both,3);

uni = squeeze(uni(ndDurInd,legInds,:));

both = squeeze(both(ndDurInd,legInds,:));

neither = squeeze(neither(ndDurInd,legInds,:));

uni = uni';
both = both';
neither = neither';

uniNeiDat = [uni(:); neither(:)];
bothNeiDat = [both(:); neither(:)];


legBoth = [repmat("R1",numFliesBoth,1); repmat("R2",numFliesBoth,1); ...
    repmat("R3",numFliesBoth,1); repmat("L1",numFliesBoth,1); ...
    repmat("L2",numFliesBoth,1); repmat("L3",numFliesBoth,1)];
legUni = [repmat("R1",numFliesUni,1); repmat("R2",numFliesUni,1); ...
    repmat("R3",numFliesUni,1); repmat("L1",numFliesUni,1); ...
    repmat("L2",numFliesUni,1); repmat("L3",numFliesUni,1)];
legNei = [repmat("R1",numFliesNei,1); repmat("R2",numFliesNei,1); ...
    repmat("R3",numFliesNei,1); repmat("L1",numFliesNei,1); ...
    repmat("L2",numFliesNei,1); repmat("L3",numFliesNei,1)];
legUniNei = [legUni; legNei];
legBothNei = [legBoth; legNei];

stimBothNei = [repmat("both",numFliesBoth * numLegs,1); ...
    repmat("neither",numFliesNei * numLegs,1)];
stimUniNei = [repmat("uni",numFliesUni * numLegs,1); ...
    repmat("neither",42,1)];

% ANOVA: uni-neither
[pUniNei, statsUniNei] = circ_hktest(deg2rad(uniNeiDat), legUniNei,...
    stimUniNei, 1, {'Leg','Stim'});

% post-hoc tests
pUniNeiPH = [];
statsUniNeiPH = [];
for i = 1:numLegs
    whichLeg = legInds(i);
    [thisP, thisStats] = circ_wwtest(deg2rad(uni(:,whichLeg)), ...
        deg2rad(neither(:,whichLeg)));
    pUniNeiPH = [pUniNeiPH; thisP];
    statsUniNeiPH = [statsUniNeiPH; thisStats];
end

% ANOVA: both-neither
[pBothNei, statsBothNei] = circ_hktest(deg2rad(bothNeiDat), legBothNei,...
    stimBothNei, 1, {'Leg','Stim'});

% post-hoc tests
pBothNeiPH = [];
statsBothNeiPH = [];
for i = 1:numLegs
    whichLeg = legInds(i);
    [thisP, thisStats] = circ_wwtest(deg2rad(both(:,whichLeg)), ...
        deg2rad(neither(:,whichLeg)));
    pBothNeiPH = [pBothNeiPH; thisP];
    statsBothgNeiPH = [statsBothNeiPH; thisStats];
end

%% DNa02>>SPARC2D-GtACR1, circular

datDir = '/Users/hyang/Dropbox (HMS)/OptoBehaviorAnalysis-Helen/AnalyzedData';
whichParam = 'stepDirections';
whichPhase = 'stance';
vels = 0;
NDs = 1.3;
yScale = [-25 25];
plotAvg = 'mean';
plotIndiv = true;
plotDiff = true;

% leg indices
legInds = 1:6;
numLegs = length(legInds);

% index for particular ND
thisNDInd = 5;

% get all fly means
% get data points
% select unilateral, 1. gray, 2. static
uni = plotOptomotorOptoDiffLegStepParam_allFlies(...
    datDir, whichParam, whichPhase, vels, NDs, yScale, plotAvg, ...
    plotIndiv, plotDiff);
% select neither, 1. gray, 2. static
neither = plotOptomotorOptoDiffLegStepParam_allFlies(...
    datDir, whichParam, whichPhase, vels, NDs, yScale, plotAvg, ...
    plotIndiv, plotDiff);
% select both, 1. gray, 2. static
both = plotOptomotorOptoDiffLegStepParam_allFlies(...
    datDir, whichParam, whichPhase, vels, NDs, yScale, plotAvg, ...
    plotIndiv, plotDiff);

% close all figures
close all

% num flies
numFliesUni = size(uni,3);
numFliesNei = size(neither,3);
numFliesBoth = size(both,3);

% get only relevant ND
uni = squeeze(uni(thisNDInd,legInds,:));
neither = squeeze(neither(thisNDInd,legInds,:));
both = squeeze(both(thisNDInd,legInds,:));

uni = uni';
both = both';
neither = neither';

uniNeiDat = [uni(:); neither(:)];
bothNeiDat = [both(:); neither(:)];


legBoth = [repmat("R1",numFliesBoth,1); repmat("R2",numFliesBoth,1); ...
    repmat("R3",numFliesBoth,1); repmat("L1",numFliesBoth,1); ...
    repmat("L2",numFliesBoth,1); repmat("L3",numFliesBoth,1)];
legUni = [repmat("R1",numFliesUni,1); repmat("R2",numFliesUni,1); ...
    repmat("R3",numFliesUni,1); repmat("L1",numFliesUni,1); ...
    repmat("L2",numFliesUni,1); repmat("L3",numFliesUni,1)];
legNei = [repmat("R1",numFliesNei,1); repmat("R2",numFliesNei,1); ...
    repmat("R3",numFliesNei,1); repmat("L1",numFliesNei,1); ...
    repmat("L2",numFliesNei,1); repmat("L3",numFliesNei,1)];
legUniNei = [legUni; legNei];
legBothNei = [legBoth; legNei];

stimBothNei = [repmat("both",numFliesBoth * numLegs,1); ...
    repmat("neither",numFliesNei * numLegs,1)];
stimUniNei = [repmat("uni",numFliesUni * numLegs,1); ...
    repmat("neither",42,1)];

% ANOVA: uni-neither
[pUniNei, statsUniNei] = circ_hktest(deg2rad(uniNeiDat), legUniNei,...
    stimUniNei, 1, {'Leg','Stim'});

% post-hoc tests
pUniNeiPH = [];
statsUniNeiPH = [];
for i = 1:numLegs
    whichLeg = legInds(i);
    [thisP, thisStats] = circ_wwtest(deg2rad(uni(:,whichLeg)), ...
        deg2rad(neither(:,whichLeg)));
    pUniNeiPH = [pUniNeiPH; thisP];
    statsUniNeiPH = [statsUniNeiPH; thisStats];
end

% ANOVA: both-neither
[pBothNei, statsBothNei] = circ_hktest(deg2rad(bothNeiDat), legBothNei,...
    stimBothNei, 1, {'Leg','Stim'});

% post-hoc tests
pBothNeiPH = [];
statsBothNeiPH = [];
for i = 1:numLegs
    whichLeg = legInds(i);
    [thisP, thisStats] = circ_wwtest(deg2rad(both(:,whichLeg)), ...
        deg2rad(neither(:,whichLeg)));
    pBothNeiPH = [pBothNeiPH; thisP];
    statsBothgNeiPH = [statsBothNeiPH; thisStats];
end

%% DNg13 IInj, circular

datDir = '/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/AnalyzedData';
whichParam = 'stepDirections';
whichPhase = 'stance';
durs = 1;
amps = [-75 100];
yScale = [-20 20];
plotAvg = 'mean';
plotIndiv = true;
plotDiff = true;

% leg indices
legInds = 1:6;
numLegs = length(legInds);

% get all fly means
allFliesMeans = plotIInjLegStepParam_allFlies(datDir, whichParam, whichPhase,...
    durs, amps, yScale, plotAvg, plotIndiv, plotDiff);

close all

% get number of flies
numFlies = size(allFliesMeans,3);

% get specific data for test
hyperpol = squeeze(allFliesMeans(2,legInds,:))';
depol = squeeze(allFliesMeans(3,legInds,:))';

allDat = [hyperpol(:); depol(:)];


% generate factor matricies
legG = [repmat("R1",numFlies,1); repmat("R2",numFlies,1); ...
    repmat("R3",numFlies,1); repmat("L1",numFlies,1); ...
    repmat("L2",numFlies,1); repmat("L3",numFlies,1)];
legG = [legG; legG];

stimG = [repmat("hyperpol",numLegs * numFlies,1); ...
    repmat("depol",numLegs * numFlies,1)];

[p, stats] = circ_hktest(deg2rad(allDat), legG, stimG, 1, {'Leg','Stim'});

% post-hoc tests
p = [];
stats = [];
for i = 1:numLegs
    whichLeg = legInds(i);
    [thisP, thisStats] = circ_wwtest(deg2rad(hyperpol(:,whichLeg)), ...
        deg2rad(depol(:,whichLeg)));
    p = [p; thisP];
    stats = [stats; thisStats];
end
