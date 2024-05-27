
function EEG = cycleLabel(phaseT)

%Find the trough of each theta cycle
[EEG.TroughPhases,EEG.TroughLocs]= findpeaks(phaseT);

%     InverseTheta = 360-phaseT;
%     [ThetaPeaks,ThetaTroughLocs] = findpeaks(InverseTheta);

EEG.Cycle = ones(1,length(phaseT));
if EEG.TroughLocs(end)~=length(phaseT)
    EEG.Cycle(EEG.TroughLocs(end)+1:length(phaseT)) = length(EEG.TroughLocs)+1;
end

for i=2:length(EEG.TroughLocs)
    EEG.Cycle(EEG.TroughLocs(i-1)+1:EEG.TroughLocs(i)) = i;
end

% double check if the cycle assignment is correct or not
recheck=find(diff(EEG.Cycle)==1);
if length(recheck)~=length(EEG.TroughLocs) || any(recheck~=EEG.TroughLocs)
    error('Theta cycle label went wrong.')
end