% Calculate peak to trough time (ms) of spike wave forms
function [peak2tr,M] = peak2trough(ch1,ch2,ch3,ch4,Fs)

% No spikes
if isempty(ch1)
    peak2tr = NaN;
    return;
end

% Mean waveforms
M1 = mean(ch1,1);
M2 = mean(ch2,1);
M3 = mean(ch3,1);
M4 = mean(ch4,1);

% Choose largest amp channel
[~,Mind] = max([max(M1)-min(M1) max(M2)-min(M2) max(M3)-min(M3) max(M4)-min(M4)]);
switch Mind
    case 1
        M = M1;
    case 2
        M = M2;
    case 3
        M = M3;
    case 4
        M = M4;
end

% Peak
[~,PeakInd] = max(M);

% Trough after peak
[~,TrInd] = min(M(PeakInd:end));
TrInd = TrInd + PeakInd -1;

% Peak to trough (ms)
peak2tr = (TrInd-PeakInd)*(1/Fs)*1000;
end