% gaussSmooth.m
%
% Function to smooth data by bidirectional convolution with 3.5 sigma
%  truncated gaussian kernal.
% Adaptation of Stephen's python function safe_interp_conv_smooth_ball_data
%  in a2lib/proc_utils
%
% INPUTS:
%   vec - vector to smooth
%   padLen - length of padding, in samples
%   kernLen - length of Gaussian kernal for smoothing, in samples
%
% OUTPUTS:
%   smoVec - smoothed vector
%
% CREATED: 8/30/21 - HHY
%
% UPDATED:
%   8/31/21 - HHY - fixed bug with values going outside of data range with
%       interpolation and resampling
%

function smoVec = gaussSmooth(vec, padLen, kernLen)

    % get padding
    prePad = getInterpPad(vec, padLen, 'start');
    postPad = getInterpPad(vec, padLen, 'end');
    
    % add padding to vector
    if iscolumn(vec) % column vector
        padVec = [prePad; vec; postPad];
    else % row vector
        padVec = [prePad vec postPad];
    end
    
    % upsample by interpolation
    upSampFac = 10; % factor to upsample by
    ind = 1:length(padVec); % original sample points
    % upsampled points
    upSampInd = linspace(ind(1),ind(end),length(padVec) * upSampFac);
    % upsampled vector
    padVecUs = interp1(ind, padVec, upSampInd, 'spline');
    
    % gaussian kernel for kernLen samples (take upsampling into account)
    kernLenUs = kernLen * upSampFac;
    gaussCutoff = 3.5; % number of sigma to extend gaussian
    % sample points for gaussian, spanning +/- 3.5 sigma
    gaussPts = linspace(-1 * gaussCutoff, gaussCutoff, kernLenUs);
    
    % kernal
    kernUs = normpdf(gaussPts);
    kernUs = kernUs / sum(kernUs); % area under curve = 1
    
    % convolve in both directions and average
    vecConvUsFwd = conv(padVecUs, kernUs, 'same');
    vecConvUsRev = fliplr(conv(fliplr(padVecUs), fliplr(kernUs), 'same'));
    % average
    vecConvUs = (vecConvUsFwd + vecConvUsRev) / 2;
    
    % downsample, by interpolation
    vecConv = interp1(upSampInd, vecConvUs, ind, 'spline');
    
    % remove padding
    smoVec = vecConv((padLen + 1):(length(vecConv)-padLen));
    
    % if vec was originally column vector, change smoVec to column vector
    if iscolumn(vec)
        smoVec = smoVec';
    end
end