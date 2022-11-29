% extractEphysKernelsSelectPData.m
%
% Function to compute kernels and autocorrelation for multiple pData files
%  from a single fly. Relates spike rate and median filtered Vm and FicTrac 
%  parameters: forward velocity, yaw velocity, yaw speed, total speed, and 
%  slide velocity.
% Operates only when fly is moving
% Works only on trials without current injection
%
% 11/22/22 - temporary function. Rewrite later to deal with multiple flies
%
% INPUTS:
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
%   autoCorrParams - struct of autocorrelation parameters
%       maxLag - maximum lag for autocorrelation, in seconds
%   
% OUTPUT:
%   kernels - struct with all the kernel data
%       fFwdVel, fYawVel, fYawSpd, fTotSpd, fSlideVel, rFwdVel, rYawVel,
%           rYawSpd, rTotSpd, rSlideVel
%       for each of the above fields: fields spikeRate and medFiltV
%       for each of the above fields: kernel (kernelLen vector),
%           varExpl (variance explained, scalar)
%   kernelParams - struct of kernel parameters (as input, but adds)
%       t - kernel times in sec
%   autoCorrParams - struct of kernel parameters (as input, but adds)
%       sampRate - sampling rate autocorrelation was computed at, is equal
%           to kernelParams.sampRate
%       lags - autocorrelation times in sec
%   autoCorr - struct with all of the autocorrelation data for
%       spikeRate, medFiltV, FwdVel, YawVel, SlideVel, YawSpd, TotSpd
%    
% CREATED: 11/22/22 - HHY
%
% UPDATED:
%   11/22/22 - HHY
%
function [kernels, kernelParams, autoCorrParams, autoCorr] = ...
    extractEphysKernelsSelectPData(kernelParams, autoCorrParams)

    % preallocate struct for kernels 
    kernelLen = (kernelParams.sampRate * 2*kernelParams.winLen) - 1;
    % max lag of autocorrelation, in samples
    autoCorrMaxLagSamp = (autoCorrParams.maxLag * autoCorrParams.sampRate);
    % length of autocorrelation is 1 more than max lag
    autoCorrLen = autoCorrMaxLagSamp + 1;

    % struct structure for spikeRate and medFiltV fields
    oneCellStrct.kernel = zeros(1, kernelLen);
    oneCellStrct.varExpl = [];
    
    % use oneCellStrct to compose struct for 1 kernel condition
    oneKernelCondStrct.spikeRate = oneCellStrct;
    oneKernelCondStrct.medFiltV = oneCellStrct;
    
    % use oneKernelCondStrct to compose full kernels struct
    kernels.fFwdVel = oneKernelCondStrct;
    kernels.fYawVel = oneKernelCondStrct;
    kernels.fSlideVel = oneKernelCondStrct;
    kernels.fYawSpd = oneKernelCondStrct;
    kernels.fTotSpd = oneKernelCondStrct;
    kernels.rFwdVel = oneKernelCondStrct;
    kernels.rYawVel = oneKernelCondStrct;
    kernels.rSlideVel = oneKernelCondStrct;
    kernels.rYawSpd = oneKernelCondStrct;
    kernels.rTotSpd = oneKernelCondStrct;
    
    % preallocate for autocorrelation
    autoCorr.spikeRate = zeros(1, autoCorrLen);
    autoCorr.medFiltV = zeros(1, autoCorrLen);
    autoCorr.FwdVel = zeros(1, autoCorrLen);
    autoCorr.YawVel = zeros(1, autoCorrLen);
    autoCorr.SlideVel = zeros(1, autoCorrLen);
    autoCorr.YawSpd = zeros(1, autoCorrLen);
    autoCorr.TotSpd = zeros(1, autoCorrLen);


    % prompt user to select pData files
    [pDataFNames, pDataPath] = uigetfile('*.mat', 'Select pData files', ...
        pDataDir(), 'MultiSelect', 'on');
    
    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(pDataFNames))
        numPDataFiles = length(pDataFNames);
    else
        numPDataFiles = 1;
    end

    % preallocate 
    % valid time vector, used to compute contribution of each trial to 
    %  fly's total
    oneCellTrialStrct.validTime = zeros(numPDataFiles,1);
    % matrix for all trial kernels
    oneCellTrialStrct.allTrialKernels = zeros(numPDataFiles, kernelLen);
    % vector for variance explained
    oneCellTrialStrct.allTrialVarExpl = zeros(numPDataFiles,1);
    
    oneCellCondStrct.spikeRate = oneCellTrialStrct;
    oneCellCondStrct.medFiltV = oneCellTrialStrct;
    
    oneFlyKernel.fFwdVel = oneCellCondStrct;
    oneFlyKernel.fYawVel = oneCellCondStrct;
    oneFlyKernel.fSlideVel = oneCellCondStrct;
    oneFlyKernel.fYawSpd = oneCellCondStrct;
    oneFlyKernel.fTotSpd = oneCellCondStrct;
    oneFlyKernel.rFwdVel = oneCellCondStrct;
    oneFlyKernel.rYawVel = oneCellCondStrct;
    oneFlyKernel.rSlideVel = oneCellCondStrct;
    oneFlyKernel.rYawSpd = oneCellCondStrct;
    oneFlyKernel.rTotSpd = oneCellCondStrct;
    
    % preallocate for autocorrelation
    oneTrialACStrct.validTime = zeros(numPDataFiles,1);
    oneTrialACStrct.allTrialAutoCorr = zeros(numPDataFiles, autoCorrLen);
    
    oneFlyAutoCorr.spikeRate = oneTrialACStrct;
    oneFlyAutoCorr.medFiltV = oneTrialACStrct;
    oneFlyAutoCorr.FwdVel = oneTrialACStrct;
    oneFlyAutoCorr.YawVel = oneTrialACStrct;
    oneFlyAutoCorr.SlideVel = oneTrialACStrct;
    oneFlyAutoCorr.YawSpd = oneTrialACStrct;
    oneFlyAutoCorr.TotSpd = oneTrialACStrct;
    
    % loop through all pData files
    for j = 1:numPDataFiles
    
        % handle whether it's a cell array or not
        if (iscell(pDataFNames))
            pDataName = pDataFNames{j};
        else
            pDataName = pDataFNames;
        end

        % save fly name as first pDataName's date and fly (12 characters)
        if (j == 1)
            flyName = pDataName(1:12);
        end
        
        pDataFullPath = [pDataPath pDataName];
        
        % get variables in pData file
        pDataMatObj = matfile(pDataFullPath);
        pDataVarsStrct = whos(pDataMatObj);
        pDataVars = struct2cell(pDataVarsStrct);
        pDataVarsNames = pDataVars(1,:); % cell array of names
    
        % check that pData has fictracProc, ephysData, and moveNotMove 
        % structs, otherwise, skip
        if (any(contains(pDataVarsNames, 'fictracProc')) && ...
            any(contains(pDataVarsNames,'ephysData')) && ...
            any(contains(pDataVarsNames,'moveNotMove')))

            % load pData
            load(pDataFullPath, 'fictracProc', 'ephysSpikes', 'moveNotMove');
            
            % turn FicTrac values to NaNs where it dropped
            
            % hack to delete zeros from dropInd, bug in when it was
            %  computed -- FIX
            fictracProc.dropInd(fictracProc.dropInd < 1) = [];

            fictracProc.fwdVel(fictracProc.dropInd) = nan;
            fictracProc.slideVel(fictracProc.dropInd) = nan;
            fictracProc.yawAngVel(fictracProc.dropInd) = nan;
            fictracProc.yawAngSpd(fictracProc.dropInd) = nan;
            fictracProc.totAngSpd(fictracProc.dropInd) = nan;

            % turn FicTrac values to NaNs where fly isn't walking
            fictracProc.fwdVel(moveNotMove.ftNotMoveInd) = nan;
            fictracProc.slideVel(moveNotMove.ftNotMoveInd) = nan;
            fictracProc.yawAngVel(moveNotMove.ftNotMoveInd) = nan;
            fictracProc.yawAngSpd(moveNotMove.ftNotMoveInd) = nan;
            fictracProc.totAngSpd(moveNotMove.ftNotMoveInd) = nan;

            % turn spikeRate and medFiltV values to NaNs where fly isn't
            % walking
            % get not moving ind
            bouts = fictracProc.t(moveNotMove.ftNotMoveBout);
            ephysNotMoveInd = convertBoutsToInd(bouts, ephysSpikes.t);

            ephysSpikes.spikeRate(ephysNotMoveInd) = nan;
            ephysSpikes.medFiltV(ephysNotMoveInd) = nan;

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
            
            
            % loop through and compute all kernels for this trial
            % loop through all behavioral variables
            for k = 1:length(behVars)

                fKernName = ['f' behVars{k}];
                rKernName = ['r' behVars{k}];
                
                % loop through all ephys (spikeRate, medFiltV)
                for l = 1:length(actVars)
                        
                    % compute forward kernel
                    [tempKern, lags, numSeg] = computeWienerKernel(...
                        eval(behVars{k}), ...
                        eval(actVars{l}), ...
                        kernelParams.sampRate, ...
                        kernelParams.winLen,...
                        kernelParams.fwdCutFreq, ...
                        kernelParams.fwdTauFreq);
                    oneFlyKernel.(fKernName).(actVars{l}).validTime(j) = ...
                        numSeg;
                    oneFlyKernel.(fKernName).(actVars{l}).allTrialKernels(j,:) = ...
                        slepianWinFilter(tempKern, ...
                        kernelParams.fwdKernelBW, ...
                        kernelParams.sampRate);
                    predResp = computeLinearPrediction(...
                        oneFlyKernel.(fKernName).(actVars{l}).allTrialKernels(j,:),...
                        eval(behVars{k}));
                    [oneFlyKernel.(fKernName).(actVars{l}).allTrialVarExpl(j),...
                        ~] = computeVarianceExplained(eval(actVars{l}), ...
                        predResp, 'poly1');

                    % compute reverse kernel
                    [tempKern, ~, numSeg] = computeWienerKernel(...
                        eval(actVars{l}), ...
                        eval(behVars{k}), ...
                        kernelParams.sampRate, ...
                        kernelParams.winLen,...
                        kernelParams.revCutFreq, ...
                        kernelParams.revTauFreq);
                    oneFlyKernel.(rKernName).(actVars{l}).validTime(j) = ...
                        numSeg;
                    oneFlyKernel.(rKernName).(actVars{l}).allTrialKernels(j,:) = ...
                        slepianWinFilter(tempKern, ...
                        kernelParams.revKernelBW, ...
                        kernelParams.sampRate);
                    predResp = computeLinearPrediction(...
                        oneFlyKernel.(rKernName).(actVars{l}).allTrialKernels(j,:),...
                        eval(actVars{l}));
                    [oneFlyKernel.(rKernName).(actVars{l}).allTrialVarExpl(j),...
                        ~] = computeVarianceExplained(...
                        eval(behVars{k}), predResp, 'poly1');  
                end  
            end
            
            % loop through and compute behavioral autocorrelations for
            %  this trial
            for k = 1:length(behVars)
                [oneFlyAutoCorr.(behVars{k}).allTrialAutoCorr(j,:),...
                    autoCorrLags] = normAutoCorrWNan(eval(behVars{k}),...
                    autoCorrMaxLagSamp);
                % valid time, number of non NaN elements
                oneFlyAutoCorr.(behVars{k}).validTime(j) = ...
                    sum(~isnan(eval(behVars{k})));
            end
            
            % loop through and compute imaging autocorrelations for
            %  this trial
            for k = 1:length(actVars)
                % if trial has this imaging type, compute autocorr
                [oneFlyAutoCorr.(actVars{k}).allTrialAutoCorr(j,:),...
                    autoCorrLags] = normAutoCorrWNan(...
                    eval(actVars{k}), autoCorrMaxLagSamp);
                oneFlyAutoCorr.(actVars{k}).validTime(j) = ...
                    length(eval(actVars{k}));
                    
            end

        end
    end

    % get average kernel for fly
    ofkFN = fieldnames(oneFlyKernel);
    % loop through all kernel types
    for m = 1:length(ofkFN)
        occFN = fieldnames(oneFlyKernel.fFwdVel);
        % loop through all img cond/cells
        for n = 1:length(occFN)
            validTime = oneFlyKernel.(ofkFN{m}).(occFN{n}).validTime(...
                ~isnan(oneFlyKernel.(ofkFN{m}).(occFN{n}).validTime));
            
            % check that this condition has any trials
            if ~isempty(validTime)
                weighting = validTime ./ (sum(validTime));

                % remove NaNs
                oneFlyKernel.(ofkFN{m}).(occFN{n}).allTrialKernels(...
                    isnan(oneFlyKernel.(ofkFN{m}).(occFN{n}).allTrialKernels)) = []; 
                oneFlyKernel.(ofkFN{m}).(occFN{n}).allTrialVarExpl(...
                    isnan(oneFlyKernel.(ofkFN{m}).(occFN{n}).allTrialVarExpl)) = [];

                % weight kernel
                weightTrialKernels = ...
                    oneFlyKernel.(ofkFN{m}).(occFN{n}).allTrialKernels .* ...
                    weighting;
                % sum to get weighted average kernel
                thisKernel = sum(weightTrialKernels, 1);
                
                % weight variance explained
                weightVarExpl = ...
                    oneFlyKernel.(ofkFN{m}).(occFN{n}).allTrialVarExpl .* ...
                    weighting;
                % sum to get weighted average variance explained
                thisVarExpl = sum(weightVarExpl);
                
            else
                thisKernel = NaN(1,kernelLen);
                thisVarExpl = nan;
            end
            
            % save into main struct
            kernels.(ofkFN{m}).(occFN{n}).kernel = thisKernel;
            kernels.(ofkFN{m}).(occFN{n}).varExpl = thisVarExpl;    
        end
    end
    
    % get average autocorrelation for fly
    ofacFN = fieldnames(oneFlyAutoCorr);
    % loop through all fields, separate autocorrelations
    for m = 1:length(ofacFN)
        validTime = oneFlyAutoCorr.(ofacFN{m}).validTime(...
            ~isnan(oneFlyAutoCorr.(ofacFN{m}).validTime));
        
        % check that this condition has any trials
        if ~isempty(validTime)
            weighting = validTime ./ (sum(validTime));
            
            % remove NaNs
            oneFlyAutoCorr.(ofacFN{m}).allTrialAutoCorr(...
                isnan(...
                oneFlyAutoCorr.(ofacFN{m}).allTrialAutoCorr)) = []; 
            
            % weight autocorrelation
            weightTrialAutoCorr = ...
                oneFlyAutoCorr.(ofacFN{m}).allTrialAutoCorr .* ...
                weighting;
            % sum to get weighted average autocorrelation
            thisAutoCorr = sum(weightTrialAutoCorr,1);
            
        else % if no trials
            thisAutoCorr = NaN(1,autoCorrLen);
        end
        
        % save into main struct
        autoCorr.(ofacFN{m}) = thisAutoCorr;
    end  

    kernelParams.t = lags;
    autoCorrParams.lags = autoCorrLags;
    autoCorrParams.t = autoCorrLags / autoCorrParams.sampRate;

end
