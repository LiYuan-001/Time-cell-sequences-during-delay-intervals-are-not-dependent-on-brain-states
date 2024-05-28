function Fig8TreadmillDelayThetaCycle_Pairwise_examplePlot(inFile,AnalyzeSes)

% specify which cells to plot out for visualization
TTListTemp = {'TT1_03','TT3_01','TT7_01','TT9_02','TT11_01','TT11_04','TT12_02','TT12_03'};

p.savePlot = 0;
p.writeToFile = 0;

% set parameter
p.delta = [1 4];
p.theta = [6 10];


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
p.MinH = 0;
p.MaxH = 0.3;
xDetail = 0:360;
MaxSpike = 3000;    % max spike number for raster
p.MaxT = 0.3;

% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'
p.spikeFs = 30000;

% spike time differences to be considered using simlar theta
p.cycleThres = 2;

% parameters for import EEG from neuralynx / matlab utilities
FieldSelectionFlags = [1 0 0 0 1];
HeaderExtractionFlag = 1;
ExtractionMode = 1;
ExtractionModeVector = [];

% -------------------------------------------------------------------------

% Read in input information
sessInfo = SessInfoImport(inFile);
close all

for i = AnalyzeSes(1:end)
    
    
    if p.writeToFile
        savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
        if ~exist(savedir2, 'dir')
            mkdir(savedir2);
        end
    end
    
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % get valid cell ind
    % time cell / non-time cell
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5;
    
    % initiate the data
    ThetaLock_LR_Pairwise.rat = sessInfo(i).animal;
    ThetaLock_LR_Pairwise.day = sessInfo(i).day;
    
 
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);

    for k = 1:clusterNum
        ThetaLock_LR_Pairwise.tList{k} = TList{k}(1:end-2);
    end
    
    TTListTempInd = ismember(ThetaLock_LR_Pairwise.tList,TTListTemp);
    TTListTempId = find(TTListTempInd == 1);
    
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    % get each phase names (no delay etc)
    sesDirs = sessInfo(i).sessDirs;
    % load eeg file
    eegInd = sprintf('%s%d%s', 'CSC',sessInfo(i).EEGch,'.ncs');
    cscFile = fullfile(sessInfo(i).mainDir,eegInd);
    % extract EEG by Nlx import
    [Timestamps,Samples,Header] = Nlx2MatCSC(cscFile, FieldSelectionFlags,HeaderExtractionFlag, ExtractionMode, ExtractionModeVector);
    p.t1 = Timestamps(1);
    
    
    % reshape EEG samples
    EEG=reshape(Samples,length(Samples(:)),1);
    % get Frequency and ADBvolts for this channel
    [p.Fs, p.ADBVolts] = ReadHeaderEEG(Header);
    % get EEG to uV
    EEG = EEG*p.ADBVolts*10^6;
    EEGt = fftbandpass(EEG,p.Fs,p.theta(1)-1,p.theta(1),p.theta(2),p.theta(2)+1);
    [phaseT,AmpT] = thetaPhase2(EEGt);  
    
    % theta phase accumulate from each cycle
    % detect theta cycles
    ThetaCycle = cycleLabel(phaseT);
    % cycle in each region and spike count
    cycleCount = length(unique(ThetaCycle.Cycle));
    [~,cycleStartIdx] = unique(ThetaCycle.Cycle);
    cycleEndIdx = cycleStartIdx-1;
    cycleEndIdx = cycleEndIdx(2:end);
    cycleEndIdx(cycleCount) = length(ThetaCycle.Cycle);
        
    % unit: sec
    EEGTs = timeStampResample(Samples,Timestamps,p.Fs)./10^6;          
    
    [ThetaStartInd,ThetaEndInd,tsStart,tsEnd] = ThetaFinder3(EEGt,EEGTs,p.Fs,0.5,0.5);   

    figure
    hold on
    
    for nn = 1:length(ThetaStartInd)
        plot(EEGTs(ThetaStartInd(nn):10:ThetaEndInd(nn)),EEGt(ThetaStartInd(nn):10:ThetaEndInd(nn))-200)
    end
    
    for j = 1:length(sesDirs)
         % load analyzed positions
        delayFile = fullfile(mainDir,sesDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        % load path zone time
        pathZoneFile = fullfile(mainDir,sesDirs{j}, 'PathZone.mat');
        load(pathZoneFile);
        % load session x,y,t
        pathDataFile = fullfile(mainDir,sesDirs{j},'pathData.mat');
        pathData = load(pathDataFile);
        % load turning directions and correctness
        tInfoFile = fullfile(mainDir,sesDirs{j}, 'trialInfo.mat');
        tInfo = load(tInfoFile);
        
        EEGIdx = EEGTs >= pathData.t(1) & EEGTs <= pathData.t(end);
        EEGIdx = find(EEGIdx>0);
        EEGIdx2 = EEGIdx(1):10:EEGIdx(end);
        if contains(sesDirs{j},'on')
            plot(EEGTs(EEGIdx2),EEG(EEGIdx2)+600,'b');
        else
            plot(EEGTs(EEGIdx2),EEG(EEGIdx2)+600,'k');
        end
        
        armLable = [];
        turnLable = [];
        % getting turning directions
        % center arm: left turn 0, right turn 1
        % return arm: left arm 0, right arm 1
        map_1D = load(fullfile(sessInfo(i).mainDir,sesDirs{j},'ratesByECLR.mat'));        
        for m = 1:sum(map_1D.ratesByECLR.valid)
            if map_1D.ratesByECLR.ECLR(m) == 1
                armLable(m) = 1;
                turnLable(m) = 0;
            elseif map_1D.ratesByECLR.ECLR(m) == 2
                armLable(m) = 0;
                turnLable(m) = 0;
            elseif map_1D.ratesByECLR.ECLR(m) == 3
                armLable(m) = 0;
                turnLable(m) = 1;
            elseif map_1D.ratesByECLR.ECLR(m) == 4
                armLable(m) = 1;
                turnLable(m) = 1;
            else
                error('Turning label ERROR');
            end
        end
%         % find Left / Right turn trials
%         leftIdx = (armLable(:,1)==1); % idx on all trials
%         rightIdx = (armLable(:,1)==2); % idx on all trials
        
        
        % -----------------------------------------------------------------
        % region start and end matrix
        % [region, trial start, trial end]
        % 1 return, 2 delay, 3 stem, 4 choice, 5 reward
        regTimeBorder = [PathZone.posStartT.Return,PathZone.posEndT.Base]; % return+base
%         % def1: delay starts at barrier
%         % def2: delay starts at entrance
%         regTimeBorder(2,:,:) = [Fig8DelayZonePos.delayPos1.startT',Fig8DelayZonePos.delayPos1.endT']; % delay
%         regTimeBorder(3,:,:) = [Fig8DelayZonePos.delayPos1.endT'+0.02,PathZone.posEndT.Center]; % stem
%         regTimeBorder(4,:,:) = [PathZone.posStartT.Choice,PathZone.posEndT.Choice]; % choice
%         regTimeBorder(5,:,:) = [PathZone.posStartT.Reward,PathZone.posEndT.Reward]; % reward
%         regionTimeSum = sum(regTimeBorder(:,:,2)-regTimeBorder(:,:,1),2);       

        % def1: delay starts at barrier
        % def2: delay starts at entrance
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
        delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
        delayTend2 = Fig8DelayZonePos.delayPos2.endT;
        
        trialNum = size(delayTstart1,2);
        if contains(sesDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
            delayTend2_2 = delayTstart2+maxT;
        elseif contains(sesDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
            delayTend2_2 = delayTstart2+maxT;
        else
            error('Delay time is wrong')
        end
        
        
        % for each trial, take down spike theta phase in delay &
        % return+base part
        % only in significant theta period
        ts.spkPhaseT_Cycle_First = cell(1,clusterNum);
        ts.spkPhaseT_Cycle_Avg = cell(1,clusterNum);
        ts.tSp = cell(1,clusterNum);
        ts.tSp2 = cell(1,clusterNum);
        for k = 1:clusterNum
            % get each spike time, change unit to msec from sec
            % ts unit: sec
            tSp = Spike_Session.(sesDirs{j}){k};
            
            if ~isempty(tSp)
                % spike phase
                [spikeEEGidx,outIdx] = spikeEEGmatch(tSp,EEGTs,1./p.Fs);
                if ~isempty(outIdx)
                    fprintf('%d%s\n',length(outIdx),' spikes outside EEG range');
                end
                
                spikePhaseT = phaseT(spikeEEGidx);
                if ~isempty(outIdx)
                    spikeEEGidx(outIdx) = NaN;
                    spikePhaseT(outIdx) = NaN;
                    tSp(outIdx) = NaN;
                end
                
                phaseTCycleIdx = spikeThetaCycle(spikeEEGidx,cycleStartIdx,cycleEndIdx);                
                spikePhaseT_Accum = spikePhaseT+360*(phaseTCycleIdx'-1);     
                            
                % for spikes in same cycles, keep first spike phase or
                % average phases from all spikes
                spkPhaseT_Cycle_First = [];
                spkPhaseT_Cycle_Avg = [];
                tSp2 = [];
                for indTemp = 1:max(phaseTCycleIdx)
                    ind = phaseTCycleIdx == indTemp;
                    spikeInd = find(ind==1);
                    if sum(ind) > 0
                        spkPhaseT_Cycle_First = [spkPhaseT_Cycle_First,spikePhaseT_Accum(spikeInd(1))];
                        spkPhaseT_Cycle_Avg = [spkPhaseT_Cycle_Avg,mean(spikePhaseT_Accum(spikeInd))];
                        tSp2 = [tSp2;tSp(spikeInd(1))];
                    end
                end
       
                ts.spkPhaseT_Cycle_First{k} = spkPhaseT_Cycle_First;
                ts.spkPhaseT_Cycle_Avg{k} = spkPhaseT_Cycle_Avg;
                ts.tSp{k} = tSp;
                ts.tSp2{k} = tSp2;
                
                if sum(ismember(TTListTempId,k))>0
                    kkk = find(ismember(TTListTempId,k)==1);
                    plot(tSp2,-1000-kkk*50,'k*');
                end  
            
            end
        end
        
        return_L_FirstPhase = cell(1,clusterNum);
        return_R_FirstPhase = cell(1,clusterNum);
        delay_L_FirstPhase = cell(1,clusterNum);
        delay_R_FirstPhase = cell(1,clusterNum);
        
        return_L_AvgPhase = cell(1,clusterNum);
        return_R_AvgPhase = cell(1,clusterNum);
        delay_L_AvgPhase = cell(1,clusterNum);
        delay_R_AvgPhase = cell(1,clusterNum);
        
        return_L_FirstTs = cell(1,clusterNum);
        return_R_FirstTs = cell(1,clusterNum);
        delay_L_FirstTs = cell(1,clusterNum);
        delay_R_FirstTs = cell(1,clusterNum);
        
        
        % each trial
        
        for m = 1:trialNum
            
            plot(delayTstart2(m):delayTend2_2(m),-800,'m.');
            
            return_spikePhaseTemp_First = cell(1,clusterNum);
            return_spikePhaseTemp_Avg = cell(1,clusterNum);
            
            delay_spikePhaseTemp_First = cell(1,clusterNum);
            delay_spikePhaseTemp_Avg = cell(1,clusterNum);
            
            return_TsTemp_First = cell(1,clusterNum);         
            delay_TsTemp_First = cell(1,clusterNum);
            
            %--------------------------------------------------------------
            % retuern + base
            thetaStartTemp1 = [];
            thetaStartTemp2 = [];
            thetaEndTemp1 = [];
            thetaEndTemp2 = [];
            % return + base spikes under theta
            % condition 1: theta start after location start
            thetaInd = tsStart >= regTimeBorder(m,1) & tsStart <= regTimeBorder(m,2);
            if sum(thetaInd) > 0
                thetaStartTemp1 = tsStart(thetaInd);
                thetaEndTemp1 = tsEnd(thetaInd);
                if thetaEndTemp1(end) > regTimeBorder(m,2)
                    thetaEndTemp1(end) = regTimeBorder(m,2);
                end
            end
            % condition2: theta start before location start
            if any(tsStart < regTimeBorder(m,1) & tsEnd > regTimeBorder(m,1))
                thetaInd = (tsStart < regTimeBorder(m,1) & tsEnd > regTimeBorder(m,1));
                thetaStartTemp2 = regTimeBorder(m,1);
                thetaEndTemp2 = tsEnd(thetaInd);                
                if thetaEndTemp2(end) > regTimeBorder(m,2)
                    thetaEndTemp2(end) = regTimeBorder(m,2);
                end
            end
            
            thetaStart = sort([thetaStartTemp1;thetaStartTemp2]);
            thetaEnd = sort([thetaEndTemp1;thetaEndTemp2]);
            
            if ~isempty(thetaStart)
                for kk = 1:length(thetaStart)
                    for k = 1:clusterNum
                        spikeInd = ts.tSp2{k}>=thetaStart(kk) & ts.tSp2{k}<=thetaEnd(kk);
                        
                        % take down phase
                        return_spikePhaseTemp_First{k} = [return_spikePhaseTemp_First{k},ts.spkPhaseT_Cycle_First{k}(spikeInd)];
                        return_spikePhaseTemp_Avg{k} = [return_spikePhaseTemp_Avg{k},ts.spkPhaseT_Cycle_Avg{k}(spikeInd)];
                        % take down ts
                        return_TsTemp_First{k} = [return_TsTemp_First{k},ts.tSp2{k}(spikeInd)'];
                    end
                end
            end
            
            %--------------------------------------------------------------
            % delay area spike under theta
            thetaStartTemp1 = [];
            thetaStartTemp2 = [];
            thetaEndTemp1 = [];
            thetaEndTemp2 = [];
            % return + base spikes under theta
            % condition 1: theta start after location start
            thetaInd = tsStart >= delayTstart2(m) & tsStart <= delayTend2_2(m);
            if sum(thetaInd) > 0
                thetaStartTemp1 = tsStart(thetaInd);
                thetaEndTemp1 = tsEnd(thetaInd);
                if thetaEndTemp1(end) > delayTend2_2(m)
                    thetaEndTemp1(end) = delayTend2_2(m);
                end
            end
            % condition2: theta start before location start
            if any(tsStart < delayTstart2(m) & tsEnd > delayTstart2(m))
                thetaInd = (tsStart < delayTstart2(m) & tsEnd > delayTstart2(m));
                thetaStartTemp2 = delayTstart2(m);
                thetaEndTemp2 = tsEnd(thetaInd);                
                if thetaEndTemp2(end) > delayTend2_2(m)
                    thetaEndTemp2(end) = delayTend2_2(m);
                end
            end
            
            thetaStart = sort([thetaStartTemp1;thetaStartTemp2]);
            thetaEnd = sort([thetaEndTemp1;thetaEndTemp2]);
            
            if ~isempty(thetaStart)
                for kk = 1:length(thetaStart)
                    for k = 1:clusterNum
                        spikeInd = ts.tSp2{k}>=thetaStart(kk) & ts.tSp2{k}<=thetaEnd(kk);
                        delay_spikePhaseTemp_First{k} = [delay_spikePhaseTemp_First{k},ts.spkPhaseT_Cycle_First{k}(spikeInd)];
                        delay_spikePhaseTemp_Avg{k} = [delay_spikePhaseTemp_Avg{k},ts.spkPhaseT_Cycle_Avg{k}(spikeInd)];
                        delay_TsTemp_First{k} = [delay_TsTemp_First{k},ts.tSp2{k}(spikeInd)'];
                    end
                end
            end   

        
            % center arm: left turn 0, right turn 1
            % return arm: left arm 0, right arm 1
            if armLable(m) == 0
                for k = 1:clusterNum
                    return_L_FirstPhase{k} = [return_L_FirstPhase{k},return_spikePhaseTemp_First{k}];
                    return_L_AvgPhase{k} = [return_L_AvgPhase{k},return_spikePhaseTemp_Avg{k}];
                    return_L_FirstTs{k} = [return_L_FirstTs{k},return_TsTemp_First{k}];
                end
            else
                for k = 1:clusterNum
                    return_R_FirstPhase{k} = [return_R_FirstPhase{k},return_spikePhaseTemp_First{k}];
                    return_R_AvgPhase{k} = [return_R_AvgPhase{k},return_spikePhaseTemp_Avg{k}];
                    return_R_FirstTs{k} = [return_R_FirstTs{k},return_TsTemp_First{k}];
                end
            end
            
            if turnLable(m) == 0
                for k = 1:clusterNum
                    delay_L_FirstPhase{k} = [delay_L_FirstPhase{k},delay_spikePhaseTemp_First{k}];
                    delay_L_AvgPhase{k} = [delay_L_AvgPhase{k},delay_spikePhaseTemp_Avg{k}];
                    delay_L_FirstTs{k} = [delay_L_FirstTs{k},delay_TsTemp_First{k}];
                end
            else
                for k = 1:clusterNum
                    delay_R_FirstPhase{k} = [delay_R_FirstPhase{k},delay_spikePhaseTemp_First{k}];
                    delay_R_AvgPhase{k} = [delay_R_AvgPhase{k},delay_spikePhaseTemp_Avg{k}];
                    delay_R_FirstTs{k} = [delay_R_FirstTs{k},delay_TsTemp_First{k}];
                end
            end
        end    
        
        phaseT_Quant.(sesDirs{j}).return_L_FirstTs = return_L_FirstTs;
        phaseT_Quant.(sesDirs{j}).return_R_FirstTs = return_R_FirstTs;
        phaseT_Quant.(sesDirs{j}).delay_L_FirstTs = delay_L_FirstTs;
        phaseT_Quant.(sesDirs{j}).delay_R_FirstTs = delay_R_FirstTs;
        
        phaseT_Quant.(sesDirs{j}).return_L_FirstPhase = return_L_FirstPhase;
        phaseT_Quant.(sesDirs{j}).return_R_FirstPhase = return_R_FirstPhase;
        phaseT_Quant.(sesDirs{j}).delay_L_FirstPhase = delay_L_FirstPhase;
        phaseT_Quant.(sesDirs{j}).delay_R_FirstPhase = delay_R_FirstPhase;
        
        phaseT_Quant.(sesDirs{j}).return_L_AvgPhase = return_L_AvgPhase;
        phaseT_Quant.(sesDirs{j}).return_R_AvgPhase = return_R_AvgPhase;
        phaseT_Quant.(sesDirs{j}).delay_L_AvgPhase = delay_L_AvgPhase;
        phaseT_Quant.(sesDirs{j}).delay_R_AvgPhase = delay_R_AvgPhase;
    end

    % start calculate pairwise phase differences
    % with two phases from 2 cells having >= 360*p.cycleThres differences,
    % do not consider it
    
    % phase differences
    pairSize = nchoosek(clusterNum,2);
    
    on_return_L_firstPhaseDiff_PairWise = cell(pairSize,1);
    off_return_L_firstPhaseDiff_PairWise = cell(pairSize,1);
    on_return_R_firstPhaseDiff_PairWise = cell(pairSize,1);
    off_return_R_firstPhaseDiff_PairWise = cell(pairSize,1);
    on_delay_L_firstPhaseDiff_PairWise = cell(pairSize,1);
    off_delay_L_firstPhaseDiff_PairWise = cell(pairSize,1);
    on_delay_R_firstPhaseDiff_PairWise = cell(pairSize,1);
    off_delay_R_firstPhaseDiff_PairWise = cell(pairSize,1);
    
    on_return_L_avgPhaseDiff_PairWise = cell(pairSize,1);
    off_return_L_avgPhaseDiff_PairWise = cell(pairSize,1);
    on_return_R_avgPhaseDiff_PairWise = cell(pairSize,1);
    off_return_R_avgPhaseDiff_PairWise = cell(pairSize,1);
    on_delay_L_avgPhaseDiff_PairWise = cell(pairSize,1);
    off_delay_L_avgPhaseDiff_PairWise = cell(pairSize,1);
    on_delay_R_avgPhaseDiff_PairWise = cell(pairSize,1);
    off_delay_R_avgPhaseDiff_PairWise = cell(pairSize,1);
    
    CellInd1 = nan(pairSize,1);
    CellInd2 = nan(pairSize,1);
    CellName1 = cell(pairSize,1);
    CellName2 = cell(pairSize,1);
            
    pairCount = 0;
    for kkk = 1:length(TTListTempId)-1
        for mmm = kkk+1:length(TTListTempId)
            % start add on pairCount
            pairCount = pairCount+1;
            
            k = TTListTempId(kkk);
            m = TTListTempId(mmm);
            
            % on return L
            firstPhaseTemp1 = [phaseT_Quant.on10_1.return_L_FirstPhase{k},phaseT_Quant.on10_2.return_L_FirstPhase{k},...
                phaseT_Quant.on30_1.return_L_FirstPhase{k},phaseT_Quant.on30_2.return_L_FirstPhase{k}];
            firstPhaseTemp2 = [phaseT_Quant.on10_1.return_L_FirstPhase{m},phaseT_Quant.on10_2.return_L_FirstPhase{m},...
                phaseT_Quant.on30_1.return_L_FirstPhase{m},phaseT_Quant.on30_2.return_L_FirstPhase{m}];
            
            avgPhaseTemp1 = [phaseT_Quant.on10_1.return_L_AvgPhase{k},phaseT_Quant.on10_2.return_L_AvgPhase{k},...
                phaseT_Quant.on30_1.return_L_AvgPhase{k},phaseT_Quant.on30_2.return_L_AvgPhase{k}];
            avgPhaseTemp2 = [phaseT_Quant.on10_1.return_L_AvgPhase{m},phaseT_Quant.on10_2.return_L_AvgPhase{m},...
                phaseT_Quant.on30_1.return_L_AvgPhase{m},phaseT_Quant.on30_2.return_L_AvgPhase{m}];
            
            tsTemp1 = [phaseT_Quant.on10_1.return_L_FirstTs{k},phaseT_Quant.on10_2.return_L_FirstTs{k},...
                phaseT_Quant.on30_1.return_L_FirstTs{k},phaseT_Quant.on30_2.return_L_FirstTs{k}];
            tsTemp2 = [phaseT_Quant.on10_1.return_L_FirstTs{m},phaseT_Quant.on10_2.return_L_FirstTs{m},...
                phaseT_Quant.on30_1.return_L_FirstTs{m},phaseT_Quant.on30_2.return_L_FirstTs{m}];


            if ~isempty(firstPhaseTemp1) && ~isempty(firstPhaseTemp2)
                phaseInd = abs(firstPhaseTemp1-firstPhaseTemp2') <= 360*p.cycleThres;
                % phaseInd2 = triu(phaseInd,1);
                phaseNum = sum(sum(phaseInd));
                if phaseNum > 0
                    firstPhaseDiff = firstPhaseTemp1-firstPhaseTemp2';
                    avgPhaseDiff = avgPhaseTemp1-avgPhaseTemp2';
                    firstPhaseDiff2 = firstPhaseDiff(phaseInd);
                    avgPhaseDiff2 = avgPhaseDiff(phaseInd);
                    
                    spkTemp1 = sum(phaseInd,1);
                    spkTemp2 = sum(phaseInd,2);
                    
                    plot(tsTemp1(spkTemp1>0),-1000-kkk*50,'r*');
                    plot(tsTemp2(spkTemp2>0),-1000-mmm*50,'r*');
            
                else
                    firstPhaseDiff2 = [];
                    avgPhaseDiff2 = [];
                end
            else
                firstPhaseDiff2 = [];
                avgPhaseDiff2 = [];
            end
            
            on_return_L_firstPhaseDiff_PairWise{pairCount} = firstPhaseDiff2;
            on_return_L_avgPhaseDiff_PairWise{pairCount} = avgPhaseDiff2;
            
            % on return R
            firstPhaseTemp1 = [phaseT_Quant.on10_1.return_R_FirstPhase{k},phaseT_Quant.on10_2.return_R_FirstPhase{k},...
                phaseT_Quant.on30_1.return_R_FirstPhase{k},phaseT_Quant.on30_2.return_R_FirstPhase{k}];
            firstPhaseTemp2 = [phaseT_Quant.on10_1.return_R_FirstPhase{m},phaseT_Quant.on10_2.return_R_FirstPhase{m},...
                phaseT_Quant.on30_1.return_R_FirstPhase{m},phaseT_Quant.on30_2.return_R_FirstPhase{m}];
            
            avgPhaseTemp1 = [phaseT_Quant.on10_1.return_R_AvgPhase{k},phaseT_Quant.on10_2.return_R_AvgPhase{k},...
                phaseT_Quant.on30_1.return_R_AvgPhase{k},phaseT_Quant.on30_2.return_R_AvgPhase{k}];
            avgPhaseTemp2 = [phaseT_Quant.on10_1.return_R_AvgPhase{m},phaseT_Quant.on10_2.return_R_AvgPhase{m},...
                phaseT_Quant.on30_1.return_R_AvgPhase{m},phaseT_Quant.on30_2.return_R_AvgPhase{m}];
            
            tsTemp1 = [phaseT_Quant.on10_1.return_R_FirstTs{k},phaseT_Quant.on10_2.return_R_FirstTs{k},...
                phaseT_Quant.on30_1.return_R_FirstTs{k},phaseT_Quant.on30_2.return_R_FirstTs{k}];
            tsTemp2 = [phaseT_Quant.on10_1.return_R_FirstTs{m},phaseT_Quant.on10_2.return_R_FirstTs{m},...
                phaseT_Quant.on30_1.return_R_FirstTs{m},phaseT_Quant.on30_2.return_R_FirstTs{m}];


            if ~isempty(firstPhaseTemp1) && ~isempty(firstPhaseTemp2)
                phaseInd = abs(firstPhaseTemp1-firstPhaseTemp2') <= 360*p.cycleThres;
                % phaseInd2 = triu(phaseInd,1);
                phaseNum = sum(sum(phaseInd));
                if phaseNum > 0
                    firstPhaseDiff = firstPhaseTemp1-firstPhaseTemp2';
                    avgPhaseDiff = avgPhaseTemp1-avgPhaseTemp2';
                    firstPhaseDiff2 = firstPhaseDiff(phaseInd);
                    avgPhaseDiff2 = avgPhaseDiff(phaseInd);
                    
                    spkTemp1 = sum(phaseInd,1);
                    spkTemp2 = sum(phaseInd,2);
                    
                    plot(tsTemp1(spkTemp1>0),-1000-kkk*50,'r*');
                    plot(tsTemp2(spkTemp2>0),-1000-mmm*50,'r*');
            
                else
                    firstPhaseDiff2 = [];
                    avgPhaseDiff2 = [];
                end
            else
                firstPhaseDiff2 = [];
                avgPhaseDiff2 = [];
            end
            
            on_return_R_firstPhaseDiff_PairWise{pairCount} = firstPhaseDiff2;
            on_return_R_avgPhaseDiff_PairWise{pairCount} = avgPhaseDiff2;
            
            % off return L
            firstPhaseTemp1 = [phaseT_Quant.off10_1.return_L_FirstPhase{k},phaseT_Quant.off10_2.return_L_FirstPhase{k},...
                phaseT_Quant.off30_1.return_L_FirstPhase{k},phaseT_Quant.off30_2.return_L_FirstPhase{k}];
            firstPhaseTemp2 = [phaseT_Quant.off10_1.return_L_FirstPhase{m},phaseT_Quant.off10_2.return_L_FirstPhase{m},...
                phaseT_Quant.off30_1.return_L_FirstPhase{m},phaseT_Quant.off30_2.return_L_FirstPhase{m}];
            
            avgPhaseTemp1 = [phaseT_Quant.off10_1.return_L_AvgPhase{k},phaseT_Quant.off10_2.return_L_AvgPhase{k},...
                phaseT_Quant.off30_1.return_L_AvgPhase{k},phaseT_Quant.off30_2.return_L_AvgPhase{k}];
            avgPhaseTemp2 = [phaseT_Quant.off10_1.return_L_AvgPhase{m},phaseT_Quant.off10_2.return_L_AvgPhase{m},...
                phaseT_Quant.off30_1.return_L_AvgPhase{m},phaseT_Quant.off30_2.return_L_AvgPhase{m}];
            
            tsTemp1 = [phaseT_Quant.off10_1.return_L_FirstTs{k},phaseT_Quant.off10_2.return_L_FirstTs{k},...
                phaseT_Quant.off30_1.return_L_FirstTs{k},phaseT_Quant.off30_2.return_L_FirstTs{k}];
            tsTemp2 = [phaseT_Quant.off10_1.return_L_FirstTs{m},phaseT_Quant.off10_2.return_L_FirstTs{m},...
                phaseT_Quant.off30_1.return_L_FirstTs{m},phaseT_Quant.off30_2.return_L_FirstTs{m}];


            if ~isempty(firstPhaseTemp1) && ~isempty(firstPhaseTemp2)
                phaseInd = abs(firstPhaseTemp1-firstPhaseTemp2') <= 360*p.cycleThres;
                % phaseInd2 = triu(phaseInd,1);
                phaseNum = sum(sum(phaseInd));
                if phaseNum > 0
                    firstPhaseDiff = firstPhaseTemp1-firstPhaseTemp2';
                    avgPhaseDiff = avgPhaseTemp1-avgPhaseTemp2';
                    firstPhaseDiff2 = firstPhaseDiff(phaseInd);
                    avgPhaseDiff2 = avgPhaseDiff(phaseInd);
                    
                    spkTemp1 = sum(phaseInd,1);
                    spkTemp2 = sum(phaseInd,2);
                    
                    plot(tsTemp1(spkTemp1>0),-1000-kkk*50,'r*');
                    plot(tsTemp2(spkTemp2>0),-1000-mmm*50,'r*');
            
                else
                    firstPhaseDiff2 = [];
                    avgPhaseDiff2 = [];
                end
            else
                firstPhaseDiff2 = [];
                avgPhaseDiff2 = [];
            end
            
            off_return_L_firstPhaseDiff_PairWise{pairCount} = firstPhaseDiff2;
            off_return_L_avgPhaseDiff_PairWise{pairCount} = avgPhaseDiff2;
            
            % off return R
            firstPhaseTemp1 = [phaseT_Quant.off10_1.return_R_FirstPhase{k},phaseT_Quant.off10_2.return_R_FirstPhase{k},...
                phaseT_Quant.off30_1.return_R_FirstPhase{k},phaseT_Quant.off30_2.return_R_FirstPhase{k}];
            firstPhaseTemp2 = [phaseT_Quant.off10_1.return_R_FirstPhase{m},phaseT_Quant.off10_2.return_R_FirstPhase{m},...
                phaseT_Quant.off30_1.return_R_FirstPhase{m},phaseT_Quant.off30_2.return_R_FirstPhase{m}];
            
            avgPhaseTemp1 = [phaseT_Quant.off10_1.return_R_AvgPhase{k},phaseT_Quant.off10_2.return_R_AvgPhase{k},...
                phaseT_Quant.off30_1.return_R_AvgPhase{k},phaseT_Quant.off30_2.return_R_AvgPhase{k}];
            avgPhaseTemp2 = [phaseT_Quant.off10_1.return_R_AvgPhase{m},phaseT_Quant.off10_2.return_R_AvgPhase{m},...
                phaseT_Quant.off30_1.return_R_AvgPhase{m},phaseT_Quant.off30_2.return_R_AvgPhase{m}];
            
            tsTemp1 = [phaseT_Quant.off10_1.return_R_FirstTs{k},phaseT_Quant.off10_2.return_R_FirstTs{k},...
                phaseT_Quant.off30_1.return_R_FirstTs{k},phaseT_Quant.off30_2.return_R_FirstTs{k}];
            tsTemp2 = [phaseT_Quant.off10_1.return_R_FirstTs{m},phaseT_Quant.off10_2.return_R_FirstTs{m},...
                phaseT_Quant.off30_1.return_R_FirstTs{m},phaseT_Quant.off30_2.return_R_FirstTs{m}];

            
            if ~isempty(firstPhaseTemp1) && ~isempty(firstPhaseTemp2)
                phaseInd = abs(firstPhaseTemp1-firstPhaseTemp2') <= 360*p.cycleThres;
                % phaseInd2 = triu(phaseInd,1);
                phaseNum = sum(sum(phaseInd));
                if phaseNum > 0
                    firstPhaseDiff = firstPhaseTemp1-firstPhaseTemp2';
                    avgPhaseDiff = avgPhaseTemp1-avgPhaseTemp2';
                    firstPhaseDiff2 = firstPhaseDiff(phaseInd);
                    avgPhaseDiff2 = avgPhaseDiff(phaseInd);
                    
                    spkTemp1 = sum(phaseInd,1);
                    spkTemp2 = sum(phaseInd,2);
                    
                    plot(tsTemp1(spkTemp1>0),-1000-kkk*50,'r*');
                    plot(tsTemp2(spkTemp2>0),-1000-mmm*50,'r*');
            
                else
                    firstPhaseDiff2 = [];
                    avgPhaseDiff2 = [];
                end
            else
                firstPhaseDiff2 = [];
                avgPhaseDiff2 = [];
            end
            
            off_return_R_firstPhaseDiff_PairWise{pairCount} = firstPhaseDiff2;
            off_return_R_avgPhaseDiff_PairWise{pairCount} = avgPhaseDiff2;
            
            % -------------------------------------------------------------
            
            % on delay L
            firstPhaseTemp1 = [phaseT_Quant.on10_1.delay_L_FirstPhase{k},phaseT_Quant.on10_2.delay_L_FirstPhase{k},...
                phaseT_Quant.on30_1.delay_L_FirstPhase{k},phaseT_Quant.on30_2.delay_L_FirstPhase{k}];
            firstPhaseTemp2 = [phaseT_Quant.on10_1.delay_L_FirstPhase{m},phaseT_Quant.on10_2.delay_L_FirstPhase{m},...
                phaseT_Quant.on30_1.delay_L_FirstPhase{m},phaseT_Quant.on30_2.delay_L_FirstPhase{m}];
            
            avgPhaseTemp1 = [phaseT_Quant.on10_1.delay_L_AvgPhase{k},phaseT_Quant.on10_2.delay_L_AvgPhase{k},...
                phaseT_Quant.on30_1.delay_L_AvgPhase{k},phaseT_Quant.on30_2.delay_L_AvgPhase{k}];
            avgPhaseTemp2 = [phaseT_Quant.on10_1.delay_L_AvgPhase{m},phaseT_Quant.on10_2.delay_L_AvgPhase{m},...
                phaseT_Quant.on30_1.delay_L_AvgPhase{m},phaseT_Quant.on30_2.delay_L_AvgPhase{m}];
            
           tsTemp1 = [phaseT_Quant.on10_1.delay_L_FirstTs{k},phaseT_Quant.on10_2.delay_L_FirstTs{k},...
                phaseT_Quant.on30_1.delay_L_FirstTs{k},phaseT_Quant.on30_2.delay_L_FirstTs{k}];
            tsTemp2 = [phaseT_Quant.on10_1.delay_L_FirstTs{m},phaseT_Quant.on10_2.delay_L_FirstTs{m},...
                phaseT_Quant.on30_1.delay_L_FirstTs{m},phaseT_Quant.on30_2.delay_L_FirstTs{m}];


            if ~isempty(firstPhaseTemp1) && ~isempty(firstPhaseTemp2)
                phaseInd = abs(firstPhaseTemp1-firstPhaseTemp2') <= 360*p.cycleThres;
                % phaseInd2 = triu(phaseInd,1);
                phaseNum = sum(sum(phaseInd));
                if phaseNum > 0
                    firstPhaseDiff = firstPhaseTemp1-firstPhaseTemp2';
                    avgPhaseDiff = avgPhaseTemp1-avgPhaseTemp2';
                    firstPhaseDiff2 = firstPhaseDiff(phaseInd);
                    avgPhaseDiff2 = avgPhaseDiff(phaseInd);
                    
                    spkTemp1 = sum(phaseInd,1);
                    spkTemp2 = sum(phaseInd,2);
                    
                    plot(tsTemp1(spkTemp1>0),-1000-kkk*50,'r*');
                    plot(tsTemp2(spkTemp2>0),-1000-mmm*50,'r*');
            
                else
                    firstPhaseDiff2 = [];
                    avgPhaseDiff2 = [];
                end
            else
                firstPhaseDiff2 = [];
                avgPhaseDiff2 = [];
            end
            on_delay_L_firstPhaseDiff_PairWise{pairCount} = firstPhaseDiff2;
            on_delay_L_avgPhaseDiff_PairWise{pairCount} = avgPhaseDiff2;
            
            % on delay R
            firstPhaseTemp1 = [phaseT_Quant.on10_1.delay_R_FirstPhase{k},phaseT_Quant.on10_2.delay_R_FirstPhase{k},...
                phaseT_Quant.on30_1.delay_R_FirstPhase{k},phaseT_Quant.on30_2.delay_R_FirstPhase{k}];
            firstPhaseTemp2 = [phaseT_Quant.on10_1.delay_R_FirstPhase{m},phaseT_Quant.on10_2.delay_R_FirstPhase{m},...
                phaseT_Quant.on30_1.delay_R_FirstPhase{m},phaseT_Quant.on30_2.delay_R_FirstPhase{m}];
            
            avgPhaseTemp1 = [phaseT_Quant.on10_1.delay_R_AvgPhase{k},phaseT_Quant.on10_2.delay_R_AvgPhase{k},...
                phaseT_Quant.on30_1.delay_R_AvgPhase{k},phaseT_Quant.on30_2.delay_R_AvgPhase{k}];
            avgPhaseTemp2 = [phaseT_Quant.on10_1.delay_R_AvgPhase{m},phaseT_Quant.on10_2.delay_R_AvgPhase{m},...
                phaseT_Quant.on30_1.delay_R_AvgPhase{m},phaseT_Quant.on30_2.delay_R_AvgPhase{m}];
            
            tsTemp1 = [phaseT_Quant.on10_1.delay_R_FirstTs{k},phaseT_Quant.on10_2.delay_R_FirstTs{k},...
                phaseT_Quant.on30_1.delay_R_FirstTs{k},phaseT_Quant.on30_2.delay_R_FirstTs{k}];
            tsTemp2 = [phaseT_Quant.on10_1.delay_R_FirstTs{m},phaseT_Quant.on10_2.delay_R_FirstTs{m},...
                phaseT_Quant.on30_1.delay_R_FirstTs{m},phaseT_Quant.on30_2.delay_R_FirstTs{m}];

            
            if ~isempty(firstPhaseTemp1) && ~isempty(firstPhaseTemp2)
                phaseInd = abs(firstPhaseTemp1-firstPhaseTemp2') <= 360*p.cycleThres;
                % phaseInd2 = triu(phaseInd,1);
                phaseNum = sum(sum(phaseInd));
                if phaseNum > 0
                    firstPhaseDiff = firstPhaseTemp1-firstPhaseTemp2';
                    avgPhaseDiff = avgPhaseTemp1-avgPhaseTemp2';
                    firstPhaseDiff2 = firstPhaseDiff(phaseInd);
                    avgPhaseDiff2 = avgPhaseDiff(phaseInd);
                    
                    spkTemp1 = sum(phaseInd,1);
                    spkTemp2 = sum(phaseInd,2);
                    
                    plot(tsTemp1(spkTemp1>0),-1000-kkk*50,'r*');
                    plot(tsTemp2(spkTemp2>0),-1000-mmm*50,'r*');
            
                else
                    firstPhaseDiff2 = [];
                    avgPhaseDiff2 = [];
                end
            else
                firstPhaseDiff2 = [];
                avgPhaseDiff2 = [];
            end
            on_delay_R_firstPhaseDiff_PairWise{pairCount} = firstPhaseDiff2;
            on_delay_R_avgPhaseDiff_PairWise{pairCount} = avgPhaseDiff2;
            
            % off delay L
            firstPhaseTemp1 = [phaseT_Quant.off10_1.delay_L_FirstPhase{k},phaseT_Quant.off10_2.delay_L_FirstPhase{k},...
                phaseT_Quant.off30_1.delay_L_FirstPhase{k},phaseT_Quant.off30_2.delay_L_FirstPhase{k}];
            firstPhaseTemp2 = [phaseT_Quant.off10_1.delay_L_FirstPhase{m},phaseT_Quant.off10_2.delay_L_FirstPhase{m},...
                phaseT_Quant.off30_1.delay_L_FirstPhase{m},phaseT_Quant.off30_2.delay_L_FirstPhase{m}];
            
            avgPhaseTemp1 = [phaseT_Quant.off10_1.delay_L_AvgPhase{k},phaseT_Quant.off10_2.delay_L_AvgPhase{k},...
                phaseT_Quant.off30_1.delay_L_AvgPhase{k},phaseT_Quant.off30_2.delay_L_AvgPhase{k}];
            avgPhaseTemp2 = [phaseT_Quant.off10_1.delay_L_AvgPhase{m},phaseT_Quant.off10_2.delay_L_AvgPhase{m},...
                phaseT_Quant.off30_1.delay_L_AvgPhase{m},phaseT_Quant.off30_2.delay_L_AvgPhase{m}];
            
           tsTemp1 = [phaseT_Quant.off10_1.delay_L_FirstTs{k},phaseT_Quant.off10_2.delay_L_FirstTs{k},...
                phaseT_Quant.off30_1.delay_L_FirstTs{k},phaseT_Quant.off30_2.delay_L_FirstTs{k}];
            tsTemp2 = [phaseT_Quant.off10_1.delay_L_FirstTs{m},phaseT_Quant.off10_2.delay_L_FirstTs{m},...
                phaseT_Quant.off30_1.delay_L_FirstTs{m},phaseT_Quant.off30_2.delay_L_FirstTs{m}];

            
            if ~isempty(firstPhaseTemp1) && ~isempty(firstPhaseTemp2)
                phaseInd = abs(firstPhaseTemp1-firstPhaseTemp2') <= 360*p.cycleThres;
                % phaseInd2 = triu(phaseInd,1);
                phaseNum = sum(sum(phaseInd));
                if phaseNum > 0
                    firstPhaseDiff = firstPhaseTemp1-firstPhaseTemp2';
                    avgPhaseDiff = avgPhaseTemp1-avgPhaseTemp2';
                    firstPhaseDiff2 = firstPhaseDiff(phaseInd);
                    avgPhaseDiff2 = avgPhaseDiff(phaseInd);
                    
                    spkTemp1 = sum(phaseInd,1);
                    spkTemp2 = sum(phaseInd,2);
                    
                    plot(tsTemp1(spkTemp1>0),-1000-kkk*50,'r*');
                    plot(tsTemp2(spkTemp2>0),-1000-mmm*50,'r*');
            
                else
                    firstPhaseDiff2 = [];
                    avgPhaseDiff2 = [];
                end
            else
                firstPhaseDiff2 = [];
                avgPhaseDiff2 = [];
            end
            off_delay_L_firstPhaseDiff_PairWise{pairCount} = firstPhaseDiff2;
            off_delay_L_avgPhaseDiff_PairWise{pairCount} = avgPhaseDiff2;
            
            % off delay R
            firstPhaseTemp1 = [phaseT_Quant.off10_1.delay_R_FirstPhase{k},phaseT_Quant.off10_2.delay_R_FirstPhase{k},...
                phaseT_Quant.off30_1.delay_R_FirstPhase{k},phaseT_Quant.off30_2.delay_R_FirstPhase{k}];
            firstPhaseTemp2 = [phaseT_Quant.off10_1.delay_R_FirstPhase{m},phaseT_Quant.off10_2.delay_R_FirstPhase{m},...
                phaseT_Quant.off30_1.delay_R_FirstPhase{m},phaseT_Quant.off30_2.delay_R_FirstPhase{m}];
            
            avgPhaseTemp1 = [phaseT_Quant.off10_1.delay_R_AvgPhase{k},phaseT_Quant.off10_2.delay_R_AvgPhase{k},...
                phaseT_Quant.off30_1.delay_R_AvgPhase{k},phaseT_Quant.off30_2.delay_R_AvgPhase{k}];
            avgPhaseTemp2 = [phaseT_Quant.off10_1.delay_R_AvgPhase{m},phaseT_Quant.off10_2.delay_R_AvgPhase{m},...
                phaseT_Quant.off30_1.delay_R_AvgPhase{m},phaseT_Quant.off30_2.delay_R_AvgPhase{m}];
            
           tsTemp1 = [phaseT_Quant.off10_1.delay_R_FirstTs{k},phaseT_Quant.off10_2.delay_R_FirstTs{k},...
                phaseT_Quant.off30_1.delay_R_FirstTs{k},phaseT_Quant.off30_2.delay_R_FirstTs{k}];
            tsTemp2 = [phaseT_Quant.off10_1.delay_R_FirstTs{m},phaseT_Quant.off10_2.delay_R_FirstTs{m},...
                phaseT_Quant.off30_1.delay_R_FirstTs{m},phaseT_Quant.off30_2.delay_R_FirstTs{m}];


            
            if ~isempty(firstPhaseTemp1) && ~isempty(firstPhaseTemp2)
                phaseInd = abs(firstPhaseTemp1-firstPhaseTemp2') <= 360*p.cycleThres;
                % phaseInd2 = triu(phaseInd,1);
                phaseNum = sum(sum(phaseInd));
                if phaseNum > 0
                    firstPhaseDiff = firstPhaseTemp1-firstPhaseTemp2';
                    avgPhaseDiff = avgPhaseTemp1-avgPhaseTemp2';
                    firstPhaseDiff2 = firstPhaseDiff(phaseInd);
                    avgPhaseDiff2 = avgPhaseDiff(phaseInd);
                    
                    spkTemp1 = sum(phaseInd,1);
                    spkTemp2 = sum(phaseInd,2);
                    
                    plot(tsTemp1(spkTemp1>0),-1000-kkk*50,'r*');
                    plot(tsTemp2(spkTemp2>0),-1000-mmm*50,'r*');
            
                else
                    firstPhaseDiff2 = [];
                    avgPhaseDiff2 = [];
                end
            else
                firstPhaseDiff2 = [];
                avgPhaseDiff2 = [];
            end
            off_delay_R_firstPhaseDiff_PairWise{pairCount} = firstPhaseDiff2;
            off_delay_R_avgPhaseDiff_PairWise{pairCount} = avgPhaseDiff2;
       
            CellInd1(pairCount) = k;
            CellInd2(pairCount) = m;
            CellName1{pairCount} = ThetaLock_LR_Pairwise.tList{k};
            CellName2{pairCount} = ThetaLock_LR_Pairwise.tList{m};
    
        end
    end
    
    clear ThetaLock_LR_Pairwise        
    close all 
    fprintf('Finished position analysis for session %d\n',i);
end

end