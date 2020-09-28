% computeKernels_pData.m
%
% Function that takes in user-specified pData file, which has been run
%  through filtering of FicTrac and spike calling for ephys data, and
%  returns kernels relating spike rate and median-filtered voltage to
%  FicTrac parameters.
%
% Note: written 9/16/20 as probably temporary function just for data
%  visualization. Will add functions that pool across trials and cells
%
% INPUTS:
%   kernelParams - struct of kernel parameters
%       winLen - length of window, in seconds that averages are computed 
%           over
%       fwdCutFreq - cutoff frequency (f_cut) in attenuation applied to
%           frequency domain filter, a la Nagel and Wilson 2011. Set to 0 
%           with tauFreq if no attenuation desired. For forward kernels
%       revCutFreq - cutoff frequency (f_cut) in attenuation applied to
%           frequency domain filter, a la Nagel and Wilson 2011. Set to 0 
%           with tauFreq if no attenuation desired. For reverse kernels
%       fwdTauFreq - f_tau in attenuation applied to frequency domain filter, 
%           a la Nagel and Wilson 2011. Set to 0 with cutFreq if no 
%           attenuation desired. For forward kernels
%       revTauFreq - f_tau in attenuation applied to frequency domain filter, 
%           a la Nagel and Wilson 2011. Set to 0 with cutFreq if no 
%           attenuation desired. For reverse kernels
%       sampRate - sampling rate to convert dF/F and FicTrac data to, and
%           to calculate kernel at
%       fwdKernelBW - full bandwidth of Slepian window filter applied to
%           forward kernels post-hoc
%       revKernelBW - full bandwidth of Slepian window filter applied to
%           reverse kernels post-hoc
%
% 
% OUTPUTS:
%   kernels - struct with all the kernel data
%       fFwdVel, fYawVel, fYawSpd, fTotSpd, fSlideVel, rFwdVel, rYawVel,
%           rYawSpd, rTotSpd, rSlideVel
%       for each of the above fields: fields spikeRate, medFiltV
%   kernelParams - struct of kernel parameters, from inputs, adds:
%       t - kernel times in sec
%
% CREATED: 9/16/20 - HHY
%
% UPDATED: 9/23/20 - HHY
%
function [kernels, kernelParams] = computeKernels_pData(...
    kernelParams)

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
        eachSub.varExp = [];
        
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
        
        % convert to same time scale
        newT = ...
            fictracProc.t(1):(1/kernelParams.sampRate):fictracProc.t(end);
        % FicTrac
        FwdVel = interp1(fictracProc.t, fictracProc.fwdVel', newT);
        SlideVel = interp1(fictracProc.t, fictracProc.slideVel', newT);
        YawVel = interp1(fictracProc.t, fictracProc.yawAngVel', newT);
        YawSpd = interp1(fictracProc.t, fictracProc.yawAngSpd' ,newT);
        TotSpd = interp1(fictracProc.t, fictracProc.totAngSpd', newT);
        
        % Ephys
        spikeRate = interp1(ephysSpikes.t, ephysSpikes.spikeRate, newT);
        medFiltV = interp1(ephysSpikes.t, ephysSpikes.medFiltV, newT);
        
        % names of behavioral variables
        behVars = {'FwdVel', 'SlideVel', 'YawVel', 'YawSpd',...
            'TotSpd'};
        actVars = {'spikeRate', 'medFiltV'};
        
        for i = 1:length(actVars)
            for k = 1:length(behVars)
                fKernName = ['f' behVars{k}];
                rKernName = ['r' behVars{k}];
                
                % compute forward kernel
                [tempKern, lags, numSeg] = computeWienerKernel(...
                    eval(behVars{k}), ...
                    eval(actVars{i}), ...
                    kernelParams.sampRate, ...
                    kernelParams.winLen,...
                    kernelParams.fwdCutFreq, ...
                    kernelParams.fwdTauFreq);
                
                kernels.(actVars{i}).(fKernName).kernel = ...                             slepianWinFilter(tempKern, ...
                    slepianWinFilter(tempKern, kernelParams.fwdKernelBW, ...
                    kernelParams.sampRate);
                
                predResp = computeLinearPrediction(...
                    kernels.(actVars{i}).(fKernName).kernel,...
                    eval(behVars{k}));
                [kernels.(actVars{i}).(fKernName).varExp,~] = ...
                    computeVarianceExplained(eval(actVars{i}), ...
                    predResp, 'poly1');
                
                % compute reverse kernel
                [tempKern, ~, numSeg] = computeWienerKernel(...
                    eval(actVars{i}), ...
                    eval(behVars{k}), ...
                    kernelParams.sampRate, ...
                    kernelParams.winLen,...
                    kernelParams.revCutFreq, ...
                    kernelParams.revTauFreq);
                
                kernels.(actVars{i}).(rKernName).kernel = ...
                    slepianWinFilter(tempKern, kernelParams.revKernelBW, ...
                    kernelParams.sampRate);
                
                predResp = computeLinearPrediction(...
                    kernels.(actVars{i}).(rKernName).kernel,...
                    eval(actVars{i}));
                
                [kernels.(actVars{i}).(rKernName).varExp,~] = ...
                    computeVarianceExplained(eval(behVars{k}), ...
                    predResp, 'poly1'); 
            end
        end
        kernelParams.t = lags;
    end
    
    
end