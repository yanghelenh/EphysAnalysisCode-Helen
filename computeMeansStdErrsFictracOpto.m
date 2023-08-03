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
%   isCirc - boolean for whether to compute circular stats
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
%   8/2/22 - HHY - modify to ignore NaNs in calculating means 
%   7/19/23 - HHY - allow for circular stats
%
function [means, stdErrs] = computeMeansStdErrsFictracOpto(reps, isCirc)
    for i = 1:size(reps,1)
        for j = 1:size(reps,2)
            thisRep = reps{i,j};

            % 5/18/22 version
%             % get indices of all reps that contain NaNs in data (from
%             %  FicTrac dropping, usually)
%             rowsWNan = [];
%             for k = 1:size(thisRep,1)
%                 thisRow = thisRep(k,:);
%                 if (any(isnan(thisRow)))
%                     rowsWNan = [rowsWNan k];
%                 end
%             end
%             % remove any rows that contain NaNs
%             thisRep(rowsWNan,:) = [];
% 
%             % get mean
%             means{i,j} = mean(thisRep,1);
% 
%             % get std error
%             stdErrs{i,j} = std(thisRep,0,1) / sqrt(size(thisRep,1));

            thisRepMean = zeros(1,size(thisRep,2));
            thisRepSEM = zeros(1,size(thisRep,2));
            % loop over each column
            for k = 1:size(thisRep,2)
                thisCol = thisRep(:,k);

                % number of not NaN elements
                numValEle = sum(~isnan(thisCol));

                % remove NaNs from this column
                thisCol(isnan(thisCol)) = [];

                % compute mean and std error for this time point
                % NaN when no valid elements
                if (numValEle ~=0)
                    if (~isCirc) % not circular stats
                        thisRepMean(k) = mean(thisCol);
                        thisRepSEM(k) = std(thisCol,0,1) / sqrt(numValEle);
                    else % circular stats
                        thisColRad = deg2rad(thisCol);
                        thisRepMean(k) = rad2deg(circ_mean(thisColRad));
                        thisRepSEM(k) = rad2deg(circ_std(thisColRad)) / ...
                            sqrt(numValEle);
                    end
                else
                    thisRepMean(k) = nan;
                    thisRepSEM(k) = nan;
                end

            end

            means{i,j} = thisRepMean;
            stdErrs{i,j} = thisRepSEM;
        end
    end
end