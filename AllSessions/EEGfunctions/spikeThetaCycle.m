function phaseTCycleIdx = spikeThetaCycle(spikeEEGidx,cycleStartIdx,cycleEndIdx)

phaseTCycleIdx = nan(length(spikeEEGidx),1);

for i = 1:length(spikeEEGidx)
    
    ind = cycleStartIdx<=spikeEEGidx(i) & cycleEndIdx>=spikeEEGidx(i);
    phaseTCycleIdx(i) = find(ind==1);

end