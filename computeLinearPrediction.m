% computeLinearPrediction.m
%
% Function that uses linear filter and input signal to produce an predicted
%  output by convolution.
% Filter and input must be on same time basis.
% Run after running computeWienerKernel()
% NOTE: NaNs okay in input, just means that predicted output will also have
%   NaNs
%
% INPUTS:
%   kernel - linear filter/kernel, to use to generate predicted output
%   input - the input that the kernel is acting on to produce the predicted
%       response
%
% OUTPUTS:
%   predResp - predicted response
%
% CREATED: 5/23/19
% UPDATED: 5/23/19
%

function predResp = computeLinearPrediction(kernel, input)
    
    % ensure that both kernel and input are row vectors
    if (~isrow(input))
        input = input';
    end
    if (~isrow(kernel))
        kernel = kernel';
    end
    
    predResp = conv(input, kernel,'same');
end
    