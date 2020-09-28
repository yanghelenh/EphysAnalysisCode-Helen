% computeWienerKernel.m
%
% Function to compute first-order Wiener kernel, as estimate of linear
%  filter, using the method described in French 1976, as in Clark et al.
%  2011 and Nagel and Wilson 2011.
% Attenuate high-frequency content following filtering method of Nagel and
%  Wilson 2011. c(w) = e^(-|w-f_cut|/f_tau) for |w|>=f_cut
% Follows method of Matlab function xcorr to get both positive and negative
%  lags.
% Deals with NaNs in input or output by excluding any segment that contains
%  any.
% Note: input and output must be same length and have same sampling rate.
%
% INPUT:
%   input - stimulus signal (i.e. signal whose autocorrelation appears in
%       the denominator of the kernel computation)
%   output - response signal (i.e. signal that is being cross-correlated in
%       the numerator of the kernel computation).
%   sampRate - sampling rate of input and output signals, in Hz
%   winLen - length of window, in seconds, that filter is computed over; 
%       returns filter of length 2*winLen-1 for negative and positive lags
%   cutFreq - cutoff frequency (f_cut) in attenuation applied to frequency
%       domain filter, a la Nagel and Wilson 2011. Set to 0 with tauFreq if
%       no attenuation desired.
%   tauFreq - f_tau in attenuation applied to frequency domain filter, a la
%       Nagel and Wilson 2011. Set to 0 with cutFreq if no attenuation
%       desired.
%
% OUTPUT:
%   kernel - time-domain estimate of kernel
%   lags - time points for each kernel value, in seconds
%   numSeg - number of segments used to compute kernel
%
% CREATED: 4/1/19
% UPDATED: 5/21/19 - HHY
%   8/28/19 - HHY - also return numSeg
%

function [kernel, lags, numSeg] = computeWienerKernel(input, output, ...
    sampRate, winLen, cutFreq, tauFreq)

    winSamps = floor(winLen * sampRate);

    % get overlapping segments of input and output of length window
    % total number of segments
    numSeg = length(input) - winSamps + 1;
    % preallocate
    inputSeg = zeros(numSeg, winSamps);
    outputSeg = inputSeg;
    
    % extract segments
    for i = 1:numSeg
        endInd = i + winSamps - 1;
        inputSeg(i,:) = input(i:endInd);
        outputSeg(i,:) = output(i:endInd);
    end
    
    % remove any segments with NaN in either input or output
    % gets indicies of segments with at least 1 NaN
    inputNaNs = find(sum(isnan(inputSeg), 2)); 
    outputNaNs = find(sum(isnan(outputSeg), 2));
    % segments that have NaNs in input or output or both
    nanInd = union(inputNaNs, outputNaNs);
    
    % remove segments with NaNs
    inputSeg(nanInd,:) = [];
    outputSeg(nanInd,:) = [];
    
    % compute Wiener Kernel
    numSeg = size(inputSeg, 1); % update numSeg after NaN segments removed
    
    % get fft of each segment, input and output
%     n = (window-1) * 2; % new input length for zero-padding
    mxl = winSamps - 1; % max lags to keep, following xcorr code
    n = 2^nextpow2(2*winSamps-1); % number of samples in fft
    inFFT = fft(inputSeg, n, 2);
    outFFT = fft(outputSeg, n, 2);
    inConj = conj(inFFT); % complex conjugate of input
    
    % first order Wiener Kernel is <Y(w)X*(w)>/<X(w)X*(w)>
    numerator = mean(outFFT .* inConj, 1);
    denominator = mean(inFFT .* inConj, 1);
    
    kernelFFT_noAtt = numerator ./ denominator;
    
    % attenuate high frequency content
    if ((cutFreq > 0) && (tauFreq > 0))
        % actual frequencies for first half
        halfFreq = sampRate * (0:(n/2))/n;    
        % frequencies for whole fft, absolute value of
        freq = [halfFreq(1:(end-1)) fliplr(halfFreq(1:(end-1)))];
        
        % weights at each frequency, following e^(-|w-f_cut|/f_tau) for 
        %  |w|>=f_cut
        filtWeights = exp(-1 * abs(freq - cutFreq) / tauFreq);
        filtWeights(abs(freq) < cutFreq) = 1;
        
        % apply attenuation to kernel
        kernelFFT = kernelFFT_noAtt .* filtWeights;
    else
        kernelFFT = kernelFFT_noAtt;
    end
    
    % low pass filter numerator
%     if (lowPassCutoff > 0) % apply low pass filtering to numerator
%         % actual frequencies for first half
%         halfFreq = sampRate * (0:(n/2))/n; 
%         % actual frequencies for whole fft
%         freq = [halfFreq(1:(end-1)) fliplr(halfFreq(1:(end-1)))];
%         % weights at each frequency, for first order filter with specified
%         %  low pass cutoff
%         filtWeights = 1 ./ (1 + (freq./(2 * pi * lowPassCutoff)));
%         
%         % low pass filter the numerator
%         lpfNumerator = numerator .* filtWeights;
%         
%         % compute 1st order wiener kernel
%         kernelFFT = lpfNumerator ./ denominator;
%           
%     else % no filtering
%         % compute 1st order wiener kernel 
%         kernelFFT = numerator ./ denominator;
%     end
    
    % transform frequency domain wiener kernel into time domain, following
    %  xcorr code
    kernelTime = real(ifft(kernelFFT,[],2));
    
    % keep only lags we want and move negative lags before positive lags
    kernel = [kernelTime(n - mxl + (1:mxl)), kernelTime(1:mxl+1)];
    
    % calculate lags
    lags = ([fliplr(1:mxl) * -1, 0:mxl]) /sampRate;
    
%     % same kernel with no filtering
%     % transform frequency domain wiener kernel into time domain
%     kernelTime_noAtt = real(ifft(kernelFFT_noAtt,[],2));
%     
%     % keep only lags we want and move negative lags before positive lags
%     kernel_noAtt = [kernelTime_noAtt(n - mxl + (1:mxl)), ...
%         kernelTime_noAtt(1:mxl+1)];
        
    
end