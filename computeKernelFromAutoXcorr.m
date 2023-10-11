% computeKernelFromAutoXcorr.m
%
% Function to take pre-computed autocorrelations and cross-correlations
%  (outputs of getSpikerateAutocorr_cell() and
%  getXCorrEphysContParam_cell()) and convert it to a first-order Wiener
%  kernel by transforming the autocorrelation and cross-correlation into
%  frequency space and then dividing the cross-correlation by the
%  autocorrelation.
% Select through GUI files for autocorrelations and for cross-correlations.
%  Matches them based on date_fly_cell name.
% Attenuate high-frequency content following filtering method of Nagel and
%  Wilson 2011. c(w) = e^(-|w-f_cut|/f_tau) for |w|>=f_cut
%
% INPUTS:
%   datDir - full path to location of autocorr and xcorr files
%   sampRate - sampling rate to interpolate autocorr and xcorr to, in Hz
%   winLen - length of window, in seconds, that filter is computed over; 
%       returns filter of length 2*winLen+1 samp for negative and positive lags
%   cutFreq - cutoff frequency (f_cut) in attenuation applied to frequency
%       domain filter, a la Nagel and Wilson 2011. Set to 0 with tauFreq if
%       no attenuation desired.
%   tauFreq - f_tau in attenuation applied to frequency domain filter, a la
%       Nagel and Wilson 2011. Set to 0 with cutFreq if no attenuation
%       desired.
%   saveFileDir - full path to folder in which to save output file
% 
% OUTPUTS:
%   none, but saves output file
%
% CREATED: 10/5/23 - HHY
%
% UPDATED:
%   10/5/23 - HHY
%
function computeKernelFromAutoXcorr(datDir, sampRate, winLen, ...
    cutFreq, tauFreq, saveFileDir)

    % generate parameters for fft and kernel computation
    winSamps = floor(winLen * sampRate);
    mxl = winSamps - 1; % max lags to keep, following xcorr code
    n = 2^nextpow2(2*winSamps-1); % number of samples in fft
    lagsK = ([fliplr(1:mxl) * -1, 0:mxl]) /sampRate;


    % prompt user to select cross-correlation files
    disp('Select cross-correlation files');
    [xCorrFNames, xCorrDirPath] = uigetfile('*.mat', ...
        'Select all cross-correlation files', datDir, 'MultiSelect', 'on');

    % prompt user to select autocorrelation files
    disp('Select autocorrelation files');
    [autoCorrFNames, autoCorrDirPath] = uigetfile('*.mat', ...
        'Select all autocorrelation files', datDir, 'MultiSelect', 'on');

    % if only 1 pData file selected, not cell array; make sure loop still
    %  works 
    if (iscell(xCorrFNames))
        numFiles = length(xCorrFNames);
    else
        numFiles = 1;
    end

    % loop through all files, assumes that they match across flies
    for i = 1:numFiles
        % handle whether it's a cell array or not
        if (iscell(xCorrFNames))
            xCorrName = xCorrFNames{i};
        else
            xCorrName = xCorrFNames;
        end

        cellName = xCorrName(1:19);

        xCorrFullPath = [xCorrDirPath filesep xCorrName];

        % load cross-correlation
        load(xCorrFullPath, 'lagsT','xCorr');
        xCorrT = lagsT;

        % find corresponding autocorrelation file
        if (iscell(autoCorrFNames))
            autoCorrName = autoCorrFNames{contains(autoCorrFNames, cellName)};

        else
            if (contains(autoCorrFNames, cellName))
                autoCorrName = autoCorrFNames;
            else
                return;
            end
        end

        autoCorrFullPath = [autoCorrDirPath filesep autoCorrName];

        % load autocorrelation
        load(autoCorrFullPath, 'lagsT', 'autoCorr');
        autoCorrT = lagsT;

        
        % interpolate autocorrelation and cross-correlation sampling time
        xCorrInt = interp1(xCorrT, xCorr, lagsK, 'linear','extrap');
        autoCorrInt = interp1(autoCorrT, autoCorr, lagsK, 'linear','extrap');

        % get FFT
        xCorrFFT = fft(xCorrInt,n);
        autoCorrFFT = fft(autoCorrInt,n);

        % get kernel, in freq 
        kernelFFT_noAtt = xCorrFFT ./ autoCorrFFT;

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

        % transform kernel into time domain
        kernelTime = real(ifft(kernelFFT));

        % keep only lags we want and move negative lags before positive lags
        kernel = [kernelTime(n - mxl + (1:mxl)), ...
            kernelTime(1:mxl+1)];

        kernel = kernel';

        % get time of peak in kernel
        [~, peakInd] = max(kernel);
        peakT = lagsK(peakInd);

        % save output file for this cell, with kernel, lags, filtering
        %  params
        outFileFullPath = [saveFileDir filesep cellName '_kernel.mat'];

        save(outFileFullPath, 'kernel', 'lagsK', 'sampRate', 'cutFreq', ...
            'tauFreq', 'winLen', 'peakT');
    end

end
