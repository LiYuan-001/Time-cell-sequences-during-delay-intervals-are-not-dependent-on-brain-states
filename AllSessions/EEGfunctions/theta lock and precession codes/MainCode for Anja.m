%% settings for analysis
% change settings based on the recording seetings
p.spikeFs = 30000; % unit recording sampling rate
p.Fs = 2000;       % LFP recording sampling rate

% theta lock settings
% Output mode, normalized spike counts (=1), spike counts (=2)
p.mode = 1;
% Degree bin 30 or 45 (default 30)
p.degBin = 30;
% von Mises fitting (=1) or not (=0)
p.fit = 1;
% plot spike raster (=1) or not (=0)
p.raster = 0;
% Sub-sampling of spikes (=1) or not (=0)
% NeuroImage2010. The pairwise phase consistency: A bias-free measure of
% rhythmic neuronal synchronization.
p.subSmaple = 1;
% Sub-sampling spike number
p.spkMin = 60;

% for theta lock figure plot
MinH = 0;
MaxH = 0.3;
xDetail = 0:360;
MaxSpike = 3000;    % max spike number for raster
MaxT = 0.3;

%% step 1
% ignore this step if you already have a method to assign theta phase
% assign phases to theta oscillation
% EEGtheta is the theta fileter EEG
% EEGt unit: usec
[phaseT,AmpT] = thetaPhase2(EEGtheta);
%% step 2
% assign spike time to theta phase based on time match
maxDiff = 10^6/p.spikeFs*p.Fs;
% EEGTs: timestamps of theta oscillation
% ts: timestamps of spikes
% ts and EEGTs units are usec here
[spikeEEGidx,outIdx] = spikeEEGmatch(ts,EEGTs,maxDiff);
spikePhaseT = phaseT(spikeEEGidx);
%% step 3
% calculate theta lock 
lock = phaseLock2(spikePhaseT,p);
% plot figure
% normalized spike count with 30deg (or 45deg) bins
if p.degBin==30
    x=15:30:345;
else
    x=22.5:45:360;
end
nT = hist(lock.spikePhaseT,x);

if p.mode==1
    nT = nT./lock.numT;
elseif p.mode==2
    disp('Output: Original spike counts')
end

figure
subplot(1,2,1)
% rose plot
CircularPlot(lock.numT,lock.spikeRadPhaseT,lock.rT,lock.angT,p,MaxT)
TITLE1 = 'Theta Lock';
title(TITLE1,'Interpreter','None');
subplot(1,2,2)
% Histgram
PlotThetaHist(nT,lock.spikePhaseT,lock.spikeRadPhaseT,xDetail,p,MinH,MaxH,lock.pT,MaxSpike,x)
xlabel('Theta phase (deg)'); ylabel('Normalized spike counts')


%% calculate theta precession
% depends on what spike you want to analyse for precession
% all spikes from a field, or each run thru a field, or all spikes
% here assuming phases from spike of intrest are spikePhase2
% spike positions are normalized by the length of the run/field/whole 
% e.g, on linear track, spkX_inField is X coord of spikes in field, field
% pos is coord for place field
% spkX_inFieldNorm = (spkX_inField-min(fieldPos))/(max(fieldPos)-min(fieldPos));

[Circ,lin]=thetaPrecess(spkPhase_inField,spkX_inFieldNorm);

% this plot part is specific to my linear track. not necessary useful 
figure
precessionPlotOrigin(spkPhase_inField,spkX_inField,fieldPos,Circ.pValue,Circ.Alpha,Circ.Phi0)