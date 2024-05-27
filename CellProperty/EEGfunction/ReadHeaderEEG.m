function [Frequency, ADBitVolts] = ReadHeaderEEG(Header)
% extract ADBvolt for that EEG
j = 1;
while isempty(strfind(Header{j}, '-ADBitVolts'))
    j=j+1;
end

if j <= length(Header)
    bvInd = strfind(Header{j}, '-ADBitVolts');
    ADBitVolts = str2double(Header{j}(11+bvInd:end));
else
    sprintf('No ADBvolts detected')
    ADBitVolts = 6.1e-8;
end

% Extrac sampling frequency for EEG
k = 1;
while isempty(strfind(Header{k}, '-SamplingFrequency'))
    k=k+1;
end

if k <= length(Header)
    bvInd = strfind(Header{k}, '-SamplingFrequency');
    Frequency = str2double(Header{k}(18+bvInd:end));
else
    sprintf('No Frequency detected')
    Frequency = 2000;
end
end
