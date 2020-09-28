% computeHeatmap_pData.m
%
% Function that takes in user-specified pData file, which has been run
%  through filtering of FicTrac and spike calling for ephys data, and
%  returns heatmap(s) relating spike rate/voltage to forward
%  velocity and yaw velocity or the counts 
%
% Note: written 9/16/20 as probably temporary function just for data
%  visualization. Will add functions that pool across trials and cells
%
% INPUTS:
%
% 
% OUTPUTS:
%   kernels - struct with all the kernel data
%       fFwdVel, fYawVel, fYawSpd, fTotSpd, fSlideVel, rFwdVel, rYawVel,
%           rYawSpd, rTotSpd, rSlideVel
%       for each of the above fields: fields spikeRate, medFiltV
%   kernelParams - struct of kernel parameters
%       winLen - length of window, in seconds that averages are computed 
%           over
%       cutFreq - cutoff frequency (f_cut) in attenuation applied to
%           frequency domain filter, a la Nagel and Wilson 2011. Set to 0 
%           with tauFreq if no attenuation desired.
%       tauFreq - f_tau in attenuation applied to frequency domain filter, 
%           a la Nagel and Wilson 2011. Set to 0 with cutFreq if no 
%           attenuation desired.
%       sampRate - sampling rate to convert dF/F and FicTrac data to, and
%           to calculate kernel at
%       fwdKernelBW - full bandwidth of Slepian window filter applied to
%           forward kernels post-hoc
%       revKernelBW - full bandwidth of Slepian window filter applied to
%           reverse kernels post-hoc
%       t - kernel times in sec
%   heatmapMat
%       spikeRate - 3D matrix: fwd vel, rot vel, t delay
%       medFiltV
%       countsMat - 3D matrix number of values per bin
%
% CREATED: 9/16/20 - HHY
%
% UPDATED: 9/16/20 - HHY
%
function [heatmapMat] = computeHeatmap_pData()

    heatmapMat = [];

    kernelParams.winLen = 0.5;
    kernelParams.cutFreq = 0;
    kernelParams.tauFreq = 0;
    kernelParams.fwdKernelBW = 100;
    kernelParams.revKernelBW = 100;
    kernelParams.sampRate = 400;

    % prompt user to select pData file
    [pDataFName, pDataPath] = uigetfile('*.mat', 'Select pData file', ...
        pDataDir());
    
    pDataFullPath = [pDataPath pDataFName];
    
    % check that selected pData contains relevant processed data
    pDatVars = struct2cell(whos('-file', pDataFullPath));
    pDatVarsNames = pDatVars(1,:);
    contProcFictrac = sum(contains(pDatVarsNames, 'fictracProc'));
    contEphysSpikes = sum(contains(pDatVarsNames, 'ephysSpikes'));
    
    if (contProcFictrac && contEphysSpikes)
        % load processed data
        load(pDataFullPath, 'fictracProc', 'ephysSpikes');
        
        % compute kernels
        % kernel length
        kernelLen = (kernelParams.sampRate * 2*kernelParams.winLen) - 1;
        
        eachSub.kernel = zeros(1, kernelLen);
%         eachSub.varExp = [];
        
        % pre-allocate
        oneKernels.fFwdVel = eachSub;
        oneKernels.fYawVel = eachSub;
        oneKernels.fSlideVel = eachSub;
        oneKernels.fYawSpd = eachSub;
        oneKernels.fTotSpd = eachSub;
        oneKernels.rFwdVel = eachSub;
        oneKernels.rYawVel = eachSub;
        oneKernels.rSlideVel = eachSub;
        oneKernels.rYawSpd = eachSub;
        oneKernels.rTotSpd = eachSub;
        
        kernels.spikeRate = oneKernels;
        kernels.medFiltV = oneKernels;
        
        % turn FicTrac values to NaNs where it dropped

        % hack to delete zeros from dropInd, bug in when it was
        %  computed -- FIX
        fictracProc.dropInd(fictracProc.dropInd < 1) = [];

        fictracProc.fwdVel(fictracProc.dropInd) = nan;
        fictracProc.slideVel(fictracProc.dropInd) = nan;
        fictracProc.yawAngVel(fictracProc.dropInd) = nan;
        fictracProc.yawAngSpd(fictracProc.dropInd) = nan;
        fictracProc.totAngSpd(fictracProc.dropInd) = nan;
        
        FwdVel = fictracProc.fwdVel';
        SlideVel = fictracProc.slideVel';
        YawVel = fictracProc.yawAngVel';
        YawSpd = fictracProc.yawAngSpd';
        TotSpd = fictracProc.totAngSpd';
        
        % convert to same time scale
        spikeRate = interp1(ephysSpikes.t, ephysSpikes.spikeRate, ...
            fictracProc.t);
        medFiltV = interp1(ephysSpikes.t, ephysSpikes.medFiltV, ...
            fictracProc.t);
        
        % names of behavioral variables
        behVars = {'FwdVel', 'SlideVel', 'YawVel', 'YawSpd',...
            'TotSpd'};
        actVars = {'spikeRate', 'medFiltV'};
        
%         for i = 1:length(actVars)
%             for k = 1:length(behVars)
%                 fKernName = ['f' behVars{k}];
%                 rKernName = ['r' behVars{k}];
%                 
%                 % compute forward kernel
%                 [tempKern, lags, numSeg] = computeWienerKernel(...
%                     eval(behVars{k}), ...
%                     eval(actVars{i}), ...
%                     kernelParams.sampRate, ...
%                     kernelParams.winLen,...
%                     kernelParams.cutFreq, ...
%                     kernelParams.tauFreq);
%                 
%                 kernels.(actVars{i}).(fKernName).kernel = ...                             slepianWinFilter(tempKern, ...
%                     slepianWinFilter(tempKern, kernelParams.fwdKernelBW, ...
%                     kernelParams.sampRate);
%                 
% %                 predResp = computeLinearPrediction(...
% %                     kernels.(actVars{i}).(fKernName).kernel,...
% %                     eval(behVars{k}));
% %                 [kernels.(actVars{i}).(fKernName).varExp,~] = ...
% %                     computeVarianceExplained(ephysSpikes(actVars{i}), ...
% %                     predResp, 'poly1');
%                 
%                 % compute reverse kernel
%                 [tempKern, ~, numSeg] = computeWienerKernel(...
%                     eval(actVars{i}), ...
%                     eval(behVars{k}), ...
%                     kernelParams.sampRate, ...
%                     kernelParams.winLen,...
%                     kernelParams.cutFreq, ...
%                     kernelParams.tauFreq);
%                 
%                 kernels.(actVars{i}).(rKernName).kernel = ...
%                     slepianWinFilter(tempKern, kernelParams.revKernelBW, ...
%                     kernelParams.sampRate);
%                 
% %                 predResp = computeLinearPrediction(...
% %                     kernels.(actVars{i}).(rKernName).kernel,...
% %                     ephysSpikes.(actVars{i}));
% %                 
% %                 [kernels.(actVars{i}).(rKernName).varExp,~] = ...
% %                     computeVarianceExplained(eval(behVars{k}), ...
% %                     predResp, 'poly1'); 
%             end
%         end
        
        dateName = pDataFName(1:6);
        flyName = pDataFName(8:12);
        cellName = pDataFName(14:19);
        trialName = pDataFName(21:27);
        
        
        % generate heat maps
        xDataName = 'yawAngVel';
        yDataName = 'fwdVel';
%         zDataName = 'spikeRate';
        zDataName = 'counts';
        xScale = [-300 300 30];
        yScale = [-5 15 30];
%         zScale = [0 150];
        zScale = [0 1000];
        minNumVals = 20;
        offsets = 0;
%         offsets = [-4000 -2000 -1500 -1000 -800 -600 -500 -400 -300 -200 -100 -50 0 100 200 500 1000 1500 2000 4000];
%         offsets = [-2000 0 2000];
        degPerMM = [];
        ttl = [dateName ' ' flyName ' ' cellName ' ' trialName];
        
        [f, heatmapMat, countsMat] = genHeatmap(YawVel, FwdVel, spikeRate,...
            xDataName, yDataName, zDataName, xScale, yScale, zScale, minNumVals,...
            offsets, degPerMM, ttl);

    else
        fprintf('%s does not contain fictracProc and/or ephysSpikes\n', ...
            pDataFName);
    end
    
%     kernelParams.t = lags;
end