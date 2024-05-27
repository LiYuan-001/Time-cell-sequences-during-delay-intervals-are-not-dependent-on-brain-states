% This code detects oscillation events based on evelop graph of filtered EEG
% if one cycle of the envelop is longer than defined Time and
% peak if over threshold, this cycle is treated as an event episode
% Li, 07-Feb-2020
function [ThetaStartInd,ThetaEndInd,tsStart,tsEnd] = ThetaFinder3(LFP,ts,Fs,ThresSD,TimeLim,Plot)
if nargin == 3
    ThresSD = 1;
    TimeLim = 1*0.5;
    Plot = 0;
elseif nargin == 4
    TimeLim = 1*0.5;
    Plot = 0;
elseif nargin == 5
    Plot = 0;
end

% Direction of vector
if size(LFP,1)>1
    LFP=LFP';
end

[EnvUpper,EnvLower] = envelope(LFP,200,'peak');

[Env_peaks,Env_peakIdx] = findpeaks(EnvUpper);

Env_mean = nanmean(EnvUpper);
Env_SD = std(EnvUpper,'omitnan');
Env_thres = Env_mean+ThresSD*Env_SD;
PEAKidx = find(EnvUpper(Env_peakIdx) >= Env_thres); % PEAK which over thre1 idx on peakidx
peakLabel = zeros(1,length(Env_peaks));
peakLabel(PEAKidx) = 1;

maxEnv = max(EnvUpper);
InverseEnv = maxEnv-EnvUpper;
[EnvTrough,Env_TroughIdx] = findpeaks(InverseEnv);
troughLabel = zeros(1,length(Env_TroughIdx))-1;

peakTroughIdx = [Env_peakIdx,Env_TroughIdx];
peakTroughLabel = [peakLabel,troughLabel];

[peakTroughIdx,sortIdx] = sort(peakTroughIdx);
peakTroughLabel = peakTroughLabel(sortIdx);

StartIdx = [];
EndIdx = [];
for i = 1:length(peakTroughLabel)-2
    if peakTroughLabel(i) == 0&& peakTroughLabel(i+1)==-1 && peakTroughLabel(i+2)==1
        StartIdx = [StartIdx,i+1];
    elseif peakTroughLabel(i) == 1 && peakTroughLabel(i+1)==-1 && peakTroughLabel(i+2)==0
        EndIdx = [EndIdx,i+1];
    end
end

if EndIdx(1) < StartIdx(1)
    StartIdx = [1,StartIdx];
end

if StartIdx(end) > EndIdx(end)
    EndIdx = [EndIdx,length(peakTroughLabel)];
end

ThetaStartIdx = peakTroughIdx(StartIdx);
ThetaEndIdx = peakTroughIdx(EndIdx);

Label1 = zeros(1,length(EnvUpper));
Label1(ThetaStartIdx) = 1;

% tsStart = ts(Env_TroughIdx(StartIdx));
% tsEnd = ts(Env_TroughIdx(EndIdx));

% if there are next start is within 0.5 sec, join two theta as one theta

nTheta = 0;
breakLength = 0.5*Fs;
i = 1;
while i <=length(ThetaStartIdx)
    Start = ThetaStartIdx(i);
    Ending = ThetaEndIdx(i);
    nTheta = nTheta+1;
    
    if Ending+breakLength <= length(LFP)
        step = breakLength;
    else
        step = length(LFP)-Ending;
    end
    
    
    while any(Label1(Ending+1:Ending+step)) == 1        
        if i+1 <= length(ThetaStartIdx) 
            i = i+1;
        end
        Ending = ThetaEndIdx(i);
        
        if Ending+breakLength <=length(LFP)
            step = breakLength;
        else
            step = length(LFP)-Ending;
        end        
    end       
    i=i+1;
    ThetaStartInd(nTheta) = Start;
    ThetaEndInd(nTheta) = Ending;
    
end

k=zeros(1,length(ThetaStartInd));

for j = 1:length(ThetaStartInd)
    if ts(ThetaEndInd(j))-ts(ThetaStartInd(j)) >= TimeLim
        k(j) = 1;
    end
end
ThetaStartInd = ThetaStartInd(k==1);
ThetaEndInd = ThetaEndInd(k==1);

tsStart = ts(ThetaStartInd);
tsEnd = ts(ThetaEndInd);
nTheta = sum(k==1);

% plot figures to see ripples
if Plot
    figure
    plot(LFP,'b')
    hold on
    plot(EnvUpper,'k')
    plot(Env_peakIdx,EnvUpper(Env_peakIdx),'r^')
    plot(Env_TroughIdx,EnvUpper(Env_TroughIdx),'g.')
end

if Plot
    for i = 1:nTheta
    plot(ThetaStartInd(i):ThetaEndInd(i),LFP(ThetaStartInd(i):ThetaEndInd(i))-500)
    end
end
hold off
end
