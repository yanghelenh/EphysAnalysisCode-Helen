% replaceValuesOutsideThresholdBound.m
% 
% REPLACEVALUESOUTSIDETHRESOLDBOUNDS replaces all values that are above or 
%  below a threshold with the last preceding value that was within the 
%  allow threshold bounds
%
% INPUT:
%   data - 1 dimentional array to be processed
%   threshold - number that defines the +/- threshold to use
%
% OUTPUT:
%   out - processed verion of data with these values replaced
%
%
% NOTE: If the first value in the array is above threshold this value will
%   be set to zero
% 
% CREATED: Yvette Fisher 1/2018
% UPDATED: 12/3/18 HHY

function [out] = replaceValuesOutsideThresholdBound(data, THRESHOLD_BOUND)

    % TODO: add in buffered region around this value that should be removed 

    % if abs value of data above threshold set to NaN
    data ( abs (data) > THRESHOLD_BOUND) = NaN;

    % check that first value in array in not Nan
    if( isnan (data(1)) )
        data(1) = 0; % set first value to zero if it was Nan
    end

    % replace Nan with proceding value
    nanIDX = find( isnan(data) );

    while(~isempty(nanIDX))
        data(nanIDX) = data(nanIDX - 1);

        % find any remaining NaN
        nanIDX  = find( isnan(data) );
    end


    out = data;
end

