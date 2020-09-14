% detectSpikes.m
%
% Function that takes in voltage trace and returns indicies of when action
%  potentials occurred. Action potentials defined by threshold crossings in
%  first derivative of voltage (dv/dt). Adacent threshold crossings must be
%  at least refractoryPeriod length of time from each other.
% 
% Returns indicies of peak, start, and end. Start and end defined by
%  specified thresholds. Peak is max between them.
% 
% If start and end vectors aren't the same length (error in calling
%  spikes), peak indicies are not returned.
%
% When two threshold crossings are too close together, assumes the latter
%  is wrong.
%
% INPUTS:
%   voltage - voltage trace
%   t - time trace
%   dvdtThresh - struct for setting thresholds on first derivative for
%           spike detection
%       start - threshold to cross from below for start of spike
%       end - threshold to cross from below for at end of spike
%   refractoryPeriod - length of time in seconds where another spike can't
%       occur
%
% OUTPUTS:
%   peakInd - indicies for peak/max of each action potential
%   startInd - indicies for start of each action potential
%
% CREATED: 9/13/20 - HHY
%
% UPDATED:
%   9/13/20 - HHY
%
function [peakInd, startInd, endInd] = detectSpikes(voltage, t, ...
    dvdtThresh, refractoryPeriod)

    % find inter-sample interval (i.e. inverse of sample rate), in seconds
    isi = median(diff(t));
    
    % find refractory period in samples
    refrPrdSamp = refractoryPeriod / isi;
    
    % find first derivative of voltage
    dvdt = diff(voltage) ./ diff(t);
    
    % find when first derivative of voltage exceeds threshold
    abvThresh = dvdt > dvdtThresh.start;
    blwThresh = dvdt < dvdtThresh.end;
    
	% find indicies for when these threshold crossings occur
    spikeStartInd = find(diff(abvThresh) == 1); % 1st to exceed thresh
    spikeEndInd = find(diff(blwThresh) == -1); % last to be below thresh
    
    % correct for taking difference
    spikeStartInd = spikeStartInd + 1;
    spikeEndInd = spikeEndInd + 1;
    
    % correct for if there are threshold crossings too close together
    %  (within refractory period)
    spikeStartIndRefrCorr = rmRefrPrdThreshCross(spikeStartInd, ...
        refrPrdSamp);
    spikeEndIndRefrCorr = rmRefrPrdThreshCross(spikeEndInd, refrPrdSamp);
    
    % first spike end index must not be prior to first spike start index
    %  delete first index if that's the case
    if (spikeEndIndRefrCorr(1) < spikeStartIndRefrCorr(1))
        spikeEndIndRefrCorr(1) = [];
    end
    
    % last spike start ind must not be later than last spike end ind; if it
    %  is, delete last element in spike start ind array
    if (spikeStartIndRefrCorr(end) > spikeEndIndRefrCorr(end))
        spikeStartIndRefrCorr(end) = [];
    end
    
    % 
    
    % length of spikeStart and spikeEnd vectors should be same at this
    %  point; if not, bigger error than can be corrected here; peakInd = []
    if (length(spikeStartIndRefrCorr) ~= length(spikeEndIndRefrCorr))
        peakInd = [];
    else % otherwise, find peaks
        % find peak of spike as max between the start and end indicies for
        %  each spike
        peakInd = zeros(size(spikeStartIndRefrCorr)); % preallocate
        for i = 1:length(spikeStartIndRefrCorr)
            % voltage snippet
            spikeVoltage = voltage(...
                spikeStartIndRefrCorr(i):spikeEndIndRefrCorr(i));
            [~, localPeakInd] = max(spikeVoltage);
            peakInd(i) = localPeakInd + spikeStartIndRefrCorr(i) - 1;
        end
    end
    
    % outputs
    startInd = spikeStartIndRefrCorr;
    endInd = spikeEndIndRefrCorr;
    
    
%     % test plotting: dots at spike start, peak, end
%     spikeStartTimes = t(startInd);
%     spikeEndTimes = t(endInd);
%     spikeStartVals = voltage(startInd);
%     spikeEndVals = voltage(endInd);
%     spikePeakTimes = t(peakInd);
%     spikePeakVals = voltage(peakInd);
%     
%     figure; 
%     plot(t, voltage);
%     hold on;
%     plot(spikeStartTimes, spikeStartVals, 'LineStyle', 'none', ...
%         'Marker', '.');
%     plot(spikeEndTimes, spikeEndVals, 'LineStyle', 'none', 'Marker', ...
%         '.', 'MarkerEdgeColor', 'magenta');
%     plot(spikePeakTimes, spikePeakVals, 'LineStyle', 'none', 'Marker', ...
%         '.', 'MarkerEdgeColor', 'green');
end
