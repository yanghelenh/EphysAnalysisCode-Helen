% computeVarianceExplained.m
%
% Function that takes in predicted response and actual response and returns
%  the variance explained given specified form of nonlinearity.
% NOTE: The predicted response and actual response should be on the same
%  time basis.
%
% INPUTS:
%   actResp - actual response
%   predResp - predicted response
%   modelType - type of curve to fit to relate actResp and predResp, must
%       be one of the model types that can be fed in as fittype to fit
%       function; as string of name
% 
% OUTPUTS:
%   varExpl - variance explained (r-squared of fit), as value from 0 to 1
%   fitObj - fit object returned by the fit function for relating actual
%       and predicted responses
%
% CREATED: 5/23/19
% UPDATED: 5/23/19
%   9/26/19 - deal with what happens if no valid data
%

function [varExpl, fitObj] = computeVarianceExplained(actResp, predResp,...
    modelType)

    % make column vectors if not already
    if (~iscolumn(actResp))
        actResp = actResp';
    end
    if (~iscolumn(predResp))
        predResp = predResp';
    end
    
    % remove NaN values from actResp and predResp
    % gets indicies of segments with at least 1 NaN
    actNaNs = find(isnan(actResp)); 
    predNaNs = find(isnan(predResp));
    % segments that have NaNs in input or output or both
    nanInd = union(actNaNs, predNaNs);
    
    % remove segments with NaNs
    actResp(nanInd) = [];
    predResp(nanInd) = [];
    
    if (~isempty(actResp) && ~isempty(predResp))
    
        [fitObj, gof] = fit(predResp, actResp, fittype(modelType));

        varExpl = gof.rsquare;
    else
        fitObj = [];
        varExpl = nan;
    end
end