% analyzeKernelsPooled.m
%
% Script to compute kernels averaged over multiple trials. Saves the
% kernels as well as the corresponding figures
%
% run extractKernels.m on all selected trials, save data, save figures
%
% UPDATED: 
%   9/27/20 - HHY

datPath = '/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/AnalyzedData/200928_kernels-autoCorr';
figPath = '/Users/hyang/Dropbox (HMS)/EphysAnalysis-Helen/Figures/200928_kernels-autoCorr';

vars = {'Exclude', 'FlyID'};
conds = {'%s~=1', 'ismember(%s,[12,15])'};

kernelParams2.winLen = 1;
kernelParams2.revTauFreq = 4;
kernelParams2.revCutFreq = 20;
kernelParams2.fwdTauFreq = 2;
kernelParams2.fwdCutFreq = 10;
kernelParams2.fwdKernelBW = 10;
kernelParams2.revKernelBW = 10;
kernelParams2.sampRate = 250;

autoCorrParams.maxLag = 1.5;

% degPerMM = 17.738631428198860;
degPerMM = 0;

whichEphys = 'spikeRate';

% load y axis scale and labels; variables: yScale, yLabels, yScaleAllDeg,
% yLabelsAllDeg
[yScale, yLabels] = loadKernelYScaleYLabels(whichEphys, degPerMM);

% autocorrelation scale bars
acYScale = {[-0.5 1], [-0.5 1], [-0.5 1], [-0.5 1], [-0.5 1], [-0.5 1],...
    [-0.5 1], [-0.5 1], [-0.5 1]};

%%
% load metadata spreadsheet
metaDat = loadMetadataSpreadsheet();

    % select appropriate pData
[~, selMetaDat] = returnSelectMetaDat(metaDat, vars, conds);

% compute all kernels, this takes a while
[kernels, autoCorr, exptNames, kernelParams, autoCorrParams] = ...
    extractKernels(...
    selMetaDat, pDataDir(), kernelParams2, autoCorrParams);
    
% save data
save([datPath filesep 'g13_kernels-autoCorr.mat'], 'selMetaDat', ...
    'kernels', 'autoCorr', 'exptNames', 'kernelParams', ...
    'autoCorrParams', 'vars', 'conds', '-v7.3');

% generate kernel figures

kernelSEMFig = plotMeanKernels(kernels, kernelParams, 'spikeRate', ...
    yScale, 'g13', yLabels, degPerMM, 1, 0, 0);
kernelIndivFig = plotMeanKernels(kernels, kernelParams, 'spikeRate',...
    yScale, 'g13', yLabels, degPerMM, 0, 1, 0);
kernelMeanOnlyFig = plotMeanKernels(kernels, kernelParams, 'spikeRate', ...
    yScale, 'g13', yLabels, degPerMM, 0, 0, 0);

% generate autocorrelation figures
autoCorrSEMFig = plotMeanAutoCorr(autoCorr, autoCorrParams, acYScale, ...
    'g13', 1, 0);
autoCorrIndivFig = plotMeanAutoCorr(autoCorr, autoCorrParams, acYScale, ...
    'g13', 0, 1);
autoCorrMeanFig = plotMeanAutoCorr(autoCorr, autoCorrParams, acYScale, ...
    'g13', 0, 0);



% save figures
saveas(kernelSEMFig, ...
    [figPath filesep 'g13' '_kernels_semFig'], 'fig');
saveas(kernelIndivFig, ...
    [figPath filesep 'g13' '_kernels_indivFig'], 'fig');
saveas(kernelMeanOnlyFig, ...
    [figPath filesep 'g13' '_kernels_meanOnlyFig'], 'fig');

saveas(autoCorrSEMFig, ...
    [figPath filesep 'g13' '_autoCorr_semFig'], 'fig');
saveas(autoCorrIndivFig, ...
    [figPath filesep 'g13' '_autoCorr_indivFig'], 'fig');
saveas(autoCorrMeanFig, ...
    [figPath filesep 'g13' '_autoCorr_meanOnlyFig'], 'fig');