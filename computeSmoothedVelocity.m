% computeSmoothedVelocity.m
%
% Helper function called by dsFiltFictrac.m that takes in downsampled
%  FicTrac position and smoothing parameters and returns smoothed position
%  and velocity. Performs Gaussian kernel smoothing on position, computes
%  velocity, and performs Gaussian kernel smoothing on velocity.
%
% INPUT:
%   pos - input position values
%   padLen - padding length for Gaussian kernel smoothing, shared for
%       smoothing on position and on velocity
%   sigmaPos - standard deviation of Gaussian kernel for position
%   sigmaVel - standard deviation of Gaussian kernel for velocity
%   sampRate - sampling rate of pos, in Hz
%
% OUTPUT:
%   smoPos - smoothed position values
%   smoVel - smoothed velocity values
%
% CREATED: 8/28/19 - HHY
%
% UPDATED: 8/28/19 - HHY
%

function [smoPos, smoVel] = computeSmoothedVelocity(pos, padLen, ...
    sigmaPos, sigmaVel, sampRate)

    % check that pos is row vector
    if ~isrow(pos)
        pos = pos';
    end
    
    % first, smooth position by convolving with gaussian kernel
    smoPosPy = py.proc_utils.safe_interp_conv_smooth_ball_data(...
        pos, padLen, sigmaPos);
    % convert from python to matlab data format
    smoPos = cell2mat(cell(smoPosPy.tolist()));

    % compute velocity
    vel = gradient(smoPos) .* sampRate;
    
    % smooth again in velocity
    smoVelPy = py.proc_utils.safe_interp_conv_smooth_ball_data(...
        vel, padLen, sigmaVel);
    % convert from python to matlab data format
    smoVel = cell2mat(cell(smoVelPy.tolist()));
end