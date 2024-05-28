% This code detects oscillation events based on evelop graph of filtered EEG
% if one cycle of the envelop is longer than defined Time and
% peak if over threshold, this cycle is treated as an event episode
% Input: LFP, filtered LFP in the interested frequency range
% Fs: Sampling frequency
% ThresSD: threshold of signal envelop being included (mean+ThresSD*SD) 
% TJoingtLim: maximum time between two events to be concatenated as one.(unit: sec)
% TLthLim: Minimum time length for an event episode.(unit: sec)
% PLOT: plot raw signal and chosen events episodes for visual check
% Output: evtOnIdx: Event onset inidx in filtered signal
% evtOffIdx: Event offset inidx in filtered signal
% eegLabel: filtered signal index in the events are labeled as 1, else 0
%
% running example: [evtOnIdx,evtOffIdx,eegLabel] = eventFinder(EEGt,2000,0.5,0.2,1,1);
% Li, 07-Feb-2020
function [evtOnIdx,evtOffIdx,eegLabel] = eventFinder(LFP_Raw,LFP_Fil,Fs,ThresSD,TJoingtLim,TLthLim,PLOT)
if nargin == 2
    ThresSD = 1;
    TJoingtLim = 0.2;
    PLOT = 0;
elseif nargin == 2
    TJoingtLim = 0.2;
    PLOT = 0;
elseif nargin == 4
    PLOT = 0;
end

% Direction of vector
if size(LFP_Fil,1)>1
    LFP_Fil=LFP_Fil';
end

eegLabel = zeros(length(LFP_Fil),1);

% get the upper boundry of filtered LFP signal envlop
% Hilbert transform
Z = hilbert(LFP_Fil);
% Wave amplitude
EnvUpper = abs(Z);
[Env_peaks,Env_peakIdx] = findpeaks(EnvUpper);
% get mean and SD of signal envlop
Env_mean = nanmean(EnvUpper);
Env_SD = std(EnvUpper,'omitnan');

% calculate threshold value for defining the signal episode
Env_thres = Env_mean+ThresSD*Env_SD;

% get peak and trough index of envelop (where later start and end of an episode
% is located on)
% envelop peak over threshold is labeled as 1
% envelop peak below threshold is labeled as 0
% envelop trough is labeled as -1

% the peak of envelop
PEAKidx = find(EnvUpper(Env_peakIdx) >= Env_thres); % PEAK which over thre1 idx on peakidx
peakLabel = zeros(1,length(Env_peaks));
peakLabel(PEAKidx) = 1;
% to calculate envelop trough index
maxEnv = max(EnvUpper);
InverseEnv = maxEnv-EnvUpper;
[EnvTrough,Env_TroughIdx] = findpeaks(InverseEnv);
troughLabel = zeros(1,length(Env_TroughIdx))-1;

% combine peak and trough index into one vector
peakTroughIdx = [Env_peakIdx,Env_TroughIdx];
peakTroughLabel = [peakLabel,troughLabel];
[peakTroughIdx,sortIdx] = sort(peakTroughIdx);
peakTroughLabel = peakTroughLabel(sortIdx);

% start detecting event
StartIdx = [];
EndIdx = [];
% find start of the over threshold peak and end of the overthreshold
% peak
for i = 1:length(peakTroughLabel)-2
    % event start is from an envelop trough
    if peakTroughLabel(i) == 0&& peakTroughLabel(i+1)==-1 && peakTroughLabel(i+2)==1
        StartIdx = [StartIdx,i+1];
        % event end is an envelop trough
    elseif peakTroughLabel(i) == 1 && peakTroughLabel(i+1)==-1 && peakTroughLabel(i+2)==0
        EndIdx = [EndIdx,i+1];
    end
end

% exclude special cases at boundries
if EndIdx(1) < StartIdx(1)
    StartIdx = [1,StartIdx];
end
if StartIdx(end) > EndIdx(end)
    EndIdx = [EndIdx,length(peakTroughLabel)];
end

eventStartIdx = peakTroughIdx(StartIdx);
eventEndIdx = peakTroughIdx(EndIdx);

Label1 = zeros(1,length(EnvUpper));
Label1(eventStartIdx) = 1;


% if there are next start is within TimeJointLim sec, join two event as one
% event
nEvent = 0;
breakLength = TJoingtLim*Fs;
i = 1;
while i <=length(eventStartIdx)
    Start = eventStartIdx(i);
    Ending = eventEndIdx(i);
    nEvent = nEvent+1;
     
    if Ending+breakLength <= length(LFP_Fil)
        step = breakLength;
    else
        step = length(LFP_Fil)-Ending;
    end
    
    % check whithin the timelimit, if there is new start
    while any(Label1(Ending+1:Ending+step)) == 1        
        if i+1 <= length(eventStartIdx) 
            i = i+1;
        end
        Ending = eventEndIdx(i);
        
        % redeine step again incase i is close to end of sequence
        if Ending+breakLength <=length(LFP_Fil)
            step = breakLength;
        else
            step = length(LFP_Fil)-Ending;
        end        
    end       
    i=i+1;
    evtOnIdx(nEvent) = Start;
    evtOffIdx(nEvent) = Ending;
    
end

% if from an event start to end is < timelength set, exclude this event
k=zeros(1,length(evtOnIdx));
for j = 1:length(evtOnIdx)
    if evtOffIdx(j)-evtOnIdx(j) >= TLthLim*Fs
        k(j) = 1;
    end
end
evtOnIdx = evtOnIdx(k==1);
evtOffIdx = evtOffIdx(k==1);
nEvent = sum(k==1);

% label eeg index inside event as 1
for m = 1:nEvent
    eegLabel(evtOnIdx(m):evtOffIdx(m)) = 1;
end

% plot figures to see ripples
if PLOT
     % plot original EEG to see singal quality
    figure
    plot(LFP_Raw+1000)   
    hold on
    plot(LFP_Fil,'b')
    plot(EnvUpper,'k','LineWidth',1)
    plot(Env_peakIdx,EnvUpper(Env_peakIdx),'r^')
    plot(Env_TroughIdx,EnvUpper(Env_TroughIdx),'g.')
    
    for i = 1:nEvent
        plot(evtOnIdx(i):evtOffIdx(i),LFP_Fil(evtOnIdx(i):evtOffIdx(i))-500)
    end
end
ylim('auto')
hold off
end