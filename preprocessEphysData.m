% preprocessEphysData.m
%
% Function to take output from preprocessUserDaq.m and extract
%  appropriately scaled and named ephys data (voltage, current, scaled
%  current/voltage, voltage command or current injection command, gain, 
%  mode, freq)
%
%
% INPUTS:
%   daqData - data collected on DAQ, with fields labeled
%   daqOutput - signal output on DAQ during experiment, with fields labeled
%   daqTime - timing vector for daqData and daqOutput
%   inputParams - input parameters from trial function (e.g. ephysRecording)
%   settings - settings struct from ephysSettings
%
% OUTPUTS:
%   ephysData - struct of appropriately scaled ephys data, named fields
%   ephysMeta - struct of ephys metadatam from decoding telegraphed output,
%       trial parameters
%
% CREATED: 2/26/20
%
% UPDATED:
%   2/26/20 - HHY
%   9/6/20 - HHY - updated to preprocess external commands (voltage or
%       current injection)
%   9/13/20 - HHY - fixed bug in processing external commands
%

function [ephysData, ephysMeta] = preprocessEphysData(daqData, ...
    daqOutput, daqTime, inputParams, settings)
    
    % decode telegraphed output
    ephysMeta.gain = decodeTelegraphedOutput(daqData.ampGain, 'gain');
    ephysMeta.freq = decodeTelegraphedOutput(daqData.ampFreq, 'freq');
    ephysMeta.mode = decodeTelegraphedOutput(daqData.ampMode, 'mode');
    
    % process non-scaled output
    % voltage in mV
    ephysData.voltage = settings.Vm.softGain .* daqData.amp10Vm;
    % current in pA
    ephysData.current = settings.I.softGain .* daqData.ampI;
    
    % process scaled output, external commands
    switch ephysMeta.mode
        case {'Track','V-Clamp'} % voltage clamp, scaled out is current
            % I = alpha * beta mV/pA (1000 is to convert from V to mV)
            ephysMeta.softGain = 1000 / ...
                (ephysMeta.gain * settings.amp.beta);
            % scaled out is current, in pA
            ephysData.scaledCurrent = ephysMeta.softGain .* ...
                daqData.ampScaledOut;
            % if there's a external voltage command
            if(isfield(daqOutput, 'ampExtCmdIn'))
                % external command is voltage, in mV
                ephysData.vCmd = settings.VOut.zeroV + ...
                    (daqOutput.ampExtCmdIn ./ settings.VOut.VConvFactor);
            end
        % current clamp, scaled out is voltage
        case {'I=0','I-Clamp Normal','I-Clamp Fast'}
            % V = alpha mV / mV (1000 for V to mV)
        	ephysMeta.softGain = 1000 / ephysMeta.gain;
            % scaled out is voltage, in mV
            ephysData.scaledVoltage = ephysMeta.softGain .* ...
                daqData.ampScaledOut;     
            % if there's a external voltage command
            if(isfield(daqOutput, 'ampExtCmdIn'))
                % external command is current, in pA
                ephysData.iInj = settings.VOut.zeroI + ...
                    (daqOutput.ampExtCmdIn ./ settings.VOut.IConvFactor);
            end
    end
    
    % copy over some metadata
    ephysMeta.startTimeStamp = inputParams.startTimeStamp;
    ephysData.t = daqTime'; % time points of recording

end