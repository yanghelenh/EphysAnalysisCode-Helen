% fit2PiPerSine.m
%
% Helper function to fit sine of period 2 pi to the input data. Returns
%  fittype model.
% Equation: a*sin(x+c) + d
% 
% INPUTS:
%   x - x values for fit, in radians
%   y - y values for fit (independent variable)
%   startVals - start values for a, c, and d, as vector
%
% OUTPUTS:
%   fittedMdl - fittype model as output
%
% CREATED: 9/18/23 - HHY
%
% UPDATED:
%   9/18/23 - HHY
%
function fittedMdl = fit2PiPerSine(x, y, startVals)
    % make sure x and y inputs are columns
    if ~iscolumn(x)
        x = x';
    end
    if ~iscolumn(y)
        y = y';
    end

    % define model: sine of period 2 pi, with amplitude (a), shift (c), and
    %  offset (d)
    mdl = fittype('a*sin(x+c) + d','indep','x');

    % get fitted model
    fittedMdl = fit(x, y, mdl, 'start', startVals);
end