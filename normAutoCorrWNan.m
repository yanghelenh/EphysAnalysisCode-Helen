% normAutoCorrWNan.m
%
% Function to compute normalized autocorrelation (Pearson's linear
%  correlation coefficient) while ignoring NaNs (using 'Rows', 'pairwise'
%  option of corr() - rho(i,j) is only computed when (i,j) lacks NaN).
% Output is 1-sided under assumption that input data is real and therefore
%  autocorrelation is real and even (symmetric about 0).
% Specify lags (in samples).
%
% INPUTS:
%   data - vector of data to compute autocorrelation on
%   mLag - max lag, in samples
%
% OUTPUTS:
%   out - normalized autocorrelation (1 sided); length is 1+nLags
%   lags - actual lags, in samples; length is 1+nLags
%
% CREATED: 9/4/19
%
% UPDATED:
%   9/4/19
%

function [out, lags] = normAutoCorrWNan(data, mLag)

    % ensure that data is a column vector
    wasCol = 1; % logical for whether data was column vector or not
    if ~isvector(data)
        disp('Input data must be a vector');
        out = [];
        lags = [];
        return;
    elseif ~iscolumn(data)
        % make row vector a column
        data = data';
        wasCol = 0;
    end
    
    % preallocate
    out = zeros(mLag+1, 1);
    
    out(1) = 1; % correlation with 0 lag is 1
    
    for i = 2:(mLag+1)
        out(i) = corr(data(i:end), data(1:(end-i+1)), 'rows', 'complete');
    end
    
    lags = (0:mLag)';
    
    % if input data was a row vector, return a row vector; otherwise, keep
    %  as column vector
    if ~wasCol
        out = out';
        lags = lags';
    end 
end