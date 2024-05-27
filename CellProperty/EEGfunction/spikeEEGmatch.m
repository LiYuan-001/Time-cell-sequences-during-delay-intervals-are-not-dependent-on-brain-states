% assign spike to most close EEG time
% it could be used for either sampling frequency is same or different
% between spikes and EEG
% Li Yuan, 01-Mar-2020, UCSD
function [spikeEEGidx,outIdx] = spikeEEGmatch(ts,EEGts,maxDiff)

% Check spike number
if isempty(ts)
    spikeEEGidx = NaN;
    outIdx = [];
    return;
end

outIdx = [];
spikeEEGidx = nan(length(ts),1);

% get start and end idx for EEG
timeDiff = EEGts-ts(1);
timeDiff = abs(timeDiff);
[minDiff,minIdx] = min(timeDiff);
if minDiff > maxDiff
    fprintf('%s\n', 'Start spike is outside range')
    startIdx = 1;
else
    startIdx = minIdx;
end

% get start and end idx for EEG
timeDiff = EEGts-ts(end);
timeDiff = abs(timeDiff);
[minDiff,minIdx] = min(timeDiff);
if minDiff > maxDiff
    fprintf('%s\n', 'End spike is outside range')
    endIdx = length(EEGts);
else
    endIdx = minIdx;
end

% use matrix substract to get difference between every EEGts and spike
if size(ts,1) == 1
    ts = ts';
end

if size(EEGts) ~= 1
    EEGts = EEGts';
end

EEGts = EEGts(startIdx:endIdx);
%--------------------------------------------------------------------------
% % MATRIX method is fast but easily exceed memory limite
% % So I will change to for loop
% % Li Yuan
% timeDiff = EEGts(startIdx:endIdx) - ts;
% % get the absolute value of the difference to find out the smallest
% % differences between EEGts and spike time
% timeDiff = abs(timeDiff);
% 
% [timeDiff3, spikeEEGidx] = min(timeDiff,[],2);
% 
% outIdx = find(timeDiff3 > maxDiff);
% if any (timeDiff3 > maxDiff)  
%     fprintf('%s%d%s\n', 'There are ',sum(outIdx),' spikes outside EEG range')
% end
% spikeEEGidx(outIdx) = [];

for i = 1:length(ts)
    timeDiff = EEGts-ts(i);
    timeDiff2 = abs(timeDiff);
    [minDiff,minIdx] = min(timeDiff2);
    if minDiff > maxDiff
        outIdx = [outIdx,i];
    else
        spikeEEGidx(i) = minIdx+startIdx-1;
    end
end
end
