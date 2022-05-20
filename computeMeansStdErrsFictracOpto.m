% computeMeansStdErrsFictracOpto.m
%
% Helper function for computeAvgFictracOpto_1Fly() that takes in cell array
%  of stimulation reps and returns 2 cell arrays of the same size, one with
%  the mean and one with the standard error across reps
%
% INPUTS:
%   reps - cell array of size # NDs x # durations, where each element is
%       matrix of reps where each row is 1 rep and each column is a time 
%       point
%
% OUTPUTS:
%   means - cell array of size # NDs x # durations, where each element is
%       vector corresponding to the mean 
%   stdErrs - cell array of size # NDs x # durations, where each element is
%       vector corresponding to the standard error
%
% CREATED: 3/29/22 - HHY
%
% UPDATED:
%   3/29/22 - HHY
%   5/18/22 - HHY - modify to remove all reps that contain NaNs from being
%       considered in calculating mean and standard error
%
function [means, stdErrs] = computeMeansStdErrsFictracOpto(reps)
    for i = 1:size(reps,1)
        for j = 1:size(reps,2)
            thisRep = reps{i,j};

            % get indices of all reps that contain NaNs in data (from
            %  FicTrac dropping, usually)
            rowsWNan = [];
            for k = 1:size(thisRep,1)
                thisRow = thisRep(k,:);
                if (any(isnan(thisRow)))
                    rowsWNan = [rowsWNan k];
                end
            end
            % remove any rows that contain NaNs
            thisRep(rowsWNan,:) = [];

            % get mean
            means{i,j} = mean(thisRep,1);

            % get std error
            stdErrs{i,j} = std(thisRep,0,1) / sqrt(size(thisRep,1));
        end
    end
end