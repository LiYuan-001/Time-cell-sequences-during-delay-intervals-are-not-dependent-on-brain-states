function lock = phaseLock2(spikePhaseT,p)

if isempty(spikePhaseT)
    lock.numT = 0;
    lock.spikePhaseT = NaN;
    lock.spikeRadPhaseT = NaN;
    lock.ppcT=NaN;
    lock.rT=NaN;
    lock.D_rT=NaN;
    lock.angT=NaN;
    lock.pT=NaN;
else
    

    % Number of spikes
    lock.numT = length(spikePhaseT);
    % Degree --> Radian spike phase
    spikeRadPhaseT = circ_ang2rad(spikePhaseT);
    
    lock.spikePhaseT = spikePhaseT;
    lock.spikeRadPhaseT = spikeRadPhaseT;
    % Pairwise phase consistency
    if lock.numT<2; lock.ppcT=NaN; else lock.ppcT = ppc(spikeRadPhaseT); end
    % Resultant vector length
    if lock.numT==0; lock.rT=NaN; else lock.rT = circ_r(spikeRadPhaseT); end
    if p.subSmaple
        if lock.numT<p.spkMin; lock.D_rT=NaN; else [lock.D_rT,~,~] = subSmaple(spikeRadPhaseT,p.spkMin); end
    end
    % Resultant vector phase [deg]
    if lock.numT==0; lock.angT=NaN; else lock.angT = circ_rad2ang(circ_mean(spikeRadPhaseT)); end
    % Rayleigh test (circ_rtest)
    if lock.numT==0; lock.pT=NaN; else lock.pT = circ_rtest(spikeRadPhaseT); end
end
end
