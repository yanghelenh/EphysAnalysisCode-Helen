% decodeTelegraphedOutput.m
%
% Function that decodes gain, frequency, and mode telegraphed from Axopatch
%  200B amplifier
%
% Adapted from Yvette's function of the same name
%
% INPUT:
%   telOut - raw output telegraphed from amplifier
%   chType - which channel type we're decoding ('gain', 'freq', 'mode')
%
% OUTPUT:
%   decVal - decoded value, appropriate to chType specified; gain and freq 
%       as double, mode as string 
%
% CREATED: 1/22/20
% UPDATED: 1/22/20 - HHY
%

function decVal = decodeTelegraphedOutput(telOut, chType)

    % value to decode
    meanTelOut = mean(telOut);

    switch chType
        case 'gain' % gain
            if meanTelOut> 0        && meanTelOut< 0.75
                decVal = 0.05;
            elseif meanTelOut> 0.75 && meanTelOut< 1.25
                decVal = 0.1;
            elseif meanTelOut> 1.25 && meanTelOut< 1.75
                decVal = 0.2;
            elseif meanTelOut> 1.75 && meanTelOut< 2.34
                decVal = 0.5;
            elseif meanTelOut> 2.34 && meanTelOut< 2.85
                decVal = 1;
            elseif meanTelOut> 2.85 && meanTelOut< 3.34
                decVal = 2;
            elseif meanTelOut> 3.34 && meanTelOut< 3.85
                decVal = 5;
            elseif meanTelOut> 3.85 && meanTelOut< 4.37
                decVal = 10;
            elseif meanTelOut> 4.37 && meanTelOut< 4.85
                decVal = 20;
            elseif meanTelOut> 4.85 && meanTelOut< 5.34
                decVal = 50;
            elseif meanTelOut> 5.34 && meanTelOut< 5.85
                decVal = 100;
            elseif meanTelOut> 5.85 && meanTelOut< 6.37
                decVal = 200;
            elseif meanTelOut> 6.37 && meanTelOut< 6.85
                decVal = 500;
            end
        case 'freq' % freq
            if meanTelOut > 0 && meanTelOut < 3
                decVal = 1;
            elseif meanTelOut > 3 && meanTelOut < 5
                decVal = 2;
            elseif meanTelOut > 5 && meanTelOut < 7
                decVal = 5;
            elseif meanTelOut > 7 && meanTelOut < 9
                decVal = 10;
            elseif meanTelOut > 9
                decVal = 100;
            end
        case 'mode' % mode
            if meanTelOut> 0 && meanTelOut< 1.5
                decVal = 'I-Clamp Fast';
            elseif meanTelOut> 1.5 && meanTelOut< 2.5
                decVal = 'I-Clamp Normal';
            elseif meanTelOut> 2.5 && meanTelOut< 3.5
                decVal = 'I=0';
            elseif meanTelOut> 3.5 && meanTelOut< 5
                decVal = 'Track';
            elseif meanTelOut> 5
                decVal = 'V-Clamp';
            end
    end

end
