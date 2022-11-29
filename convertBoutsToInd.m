% convertBoutsToInd.m
%
% Helper function that takes in start and end times of bouts and a time
%  vector and returns all indices of times that occur during that bout
%
% INPUTS:
%   bouts - n by 2 matrix of start and end times of bouts (n bouts), 1st
%       column is start time, 2nd is end time
%   t - time vector for which to return indices
%
% OUTPUTS:
%   ind - indices into t within bouts
%
% CREATED: 11/22/22 - HHY
%
% UPDATED:
%   11/22/22 - HHY
%
function ind = convertBoutsToInd(bouts, t)
    % initialize ind
    ind = [];

    for i = 1:size(bouts,1)
        theseInd = find((t > bouts(i,1)) & (t <= bouts(i,2)));

        ind = [ind; theseInd];
    end
end