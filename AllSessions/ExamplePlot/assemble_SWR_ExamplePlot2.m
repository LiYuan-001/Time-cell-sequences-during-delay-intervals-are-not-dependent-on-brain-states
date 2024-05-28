% example for population activity in delay area for manuscript
% example for popilation activity in delay and reward area
close all

p.savePlot = 0;
p.writeToFile = 0;

% set cofire time window, unit: sec
p.timeWindow = 200./10^3;
p.timeWindowIncrement = 20./10^3;
p.cellThreshold = 1/5;
% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'

% wavelet
p.space = 'log'; % 'log' or 'lin'  spacing of f's
p.frange = [100 250];
p.nCycles = 5; % cycles for morelet wavelet
p.dt = 0.05; 
p.nfreqs = 100;

% Read in input information
sessInfo = SessInfoImport('V:\LiYuan\Codes\Fig8MazeTreadmill_V2\Fig8Treadmill_OnOff.xlsx');

% EEG related
% parameters for import EEG from neuralynx / matlab utilities
FieldSelectionFlags = [1 0 0 0 1];
HeaderExtractionFlag = 1;
ExtractionMode = 1;
ExtractionModeVector = [];

for i = 10
    
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    %% spike data
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);

    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);

    % get valid cell ind
    % rate [0.1 10] hz on fig 8 maze
    % time cell / non-time cell
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    validNum = sum(rateLabel);
    validCellInd = find(rateLabel==1);
    
    %% assemblies data 
    % load assembly file
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_DelayOnOff_LR-25ms.mat');
    load(assemblyFile);
    
    patNum_on = CellAssembly_DelayLR.DelayOn.patNum; 
    AssmblPtrnCellIDs_on = CellAssembly_DelayLR.DelayOn.AssmblPtrnCellIDs;
    AssmblWght_on = CellAssembly_DelayLR.DelayOn.AssmblWght;
    AssmblStrength_on = CellAssembly_DelayLR.DelayOn.AssmblStrength;
    event_Time_on = CellAssembly_DelayLR.DelayOn.event_Time;
    event_strength_on = CellAssembly_DelayLR.DelayOn.event_strength;
    event_Num_on = CellAssembly_DelayLR.DelayOn.event_Num;
    
    patNum_off = CellAssembly_DelayLR.DelayOff.patNum; 
    AssmblPtrnCellIDs_off = CellAssembly_DelayLR.DelayOff.AssmblPtrnCellIDs;
    AssmblWght_off = CellAssembly_DelayLR.DelayOff.AssmblWght;
    AssmblStrength_off = CellAssembly_DelayLR.DelayOff.AssmblStrength;
    event_Time_off = CellAssembly_DelayLR.DelayOff.event_Time;
    event_strength_off = CellAssembly_DelayLR.DelayOff.event_strength;
    event_Num_off = CellAssembly_DelayLR.DelayOff.event_Num;
    
    binTime = CellAssembly_DelayLR.DelayOff.binTime;
    %% import EEG
    % load eeg file
    EEGch = 4;
    eegInd = sprintf('%s%d%s', 'CSC',EEGch,'.ncs');
    cscFile = fullfile(mainDir,eegInd);
    % extract EEG by Nlx import
    [Timestamps,Samples,Header] = Nlx2MatCSC(cscFile, FieldSelectionFlags,HeaderExtractionFlag, ExtractionMode, ExtractionModeVector);
    p.t1 = Timestamps(1);    
    % reshape EEG samples
    EEG = reshape(Samples,length(Samples(:)),1);
    % get Frequency and ADBvolts for this channel
    [p.Fs, p.ADBVolts] = ReadHeaderEEG(Header);
    % reasample timestamps;
    % unit: sec
    EEGTs = timeStampResample(Samples,Timestamps,p.Fs)./10^6;   
    % get EEG to uV
    EEG = EEG*p.ADBVolts*10^6;
    % filter sharp wave ripple band out
    EEG_SWR = fftbandpass(EEG,p.Fs,149,150,250,251); % ripples
    
    % block 7 trial 6 for treadmill on
    % block 8 trial 1 for treadmill off 
    
        
    for j = 4 % session
        m = 3; % trial
        
        
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        % load turning directions and correctness
        tInfoFile = fullfile(mainDir,sessDirs{j}, 'trialInfo.mat');
        tInfo = load(tInfoFile);
        % -----------------------------------------------------------------
        % use the mode of delay is always 10 sec or 30 sec
        % rather than identify real time in delay zone
        % to make plot cleaner
        % I can change to delay time in delay zone later
        % Li Yuan, 19-Aug-2020, UCSD
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
        
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
            rateBin1 = 0:p.timeWindowIncrement:maxT-p.timeWindow-p.timeWindowIncrement;
            rateBin2 = p.timeWindow:p.timeWindowIncrement:maxT-p.timeWindowIncrement;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
            rateBin1 = 0:p.timeWindowIncrement:maxT-p.timeWindow;
            rateBin2 = p.timeWindow:p.timeWindowIncrement:maxT;
        else
            error('Delay time is wrong')
        end
        
        
        % initiate the figure to plot spikes
        h = figure(1);
        h.Position = [100,100,1200,800];
        % plot assemblies
        h6=subplot(6,1,6);
        for k = 1:patNum_off
            colorValue = [1-k*0.1,k*0.1,0];
            startTime = delayTstart1(m);
            endTime = delayTend1_2(m);
            binInd = binTime>=startTime & binTime<=endTime;

            AssmblStrength = AssmblStrength_off(k,binInd);
            plot(binTime(binInd)-startTime,log2(AssmblStrength),'Color',colorValue);
            hold on
        end
        ylim([-2 10])
        ylabel('log2 (Strength)')
        
        
        % plot EEG
        % delay area
        [startTimeDifference,delaystartIdx] = min(abs(EEGTs-delayTstart1(m)));
        [endTimeDifference,delayendIdx] = min(abs(EEGTs-delayTend1_2(m)));
        % accumulate delay Idx from each trial
        delayIdx = (delaystartIdx:delayendIdx);
        if any([startTimeDifference,endTimeDifference] > (EEGTs(2)-EEGTs(1)))
            error('EEG could not find matching position')
        end
        
        lfp.data = EEG(delayIdx);
        lfp.timestamps = EEGTs(delayIdx);
        lfp.samplingRate = p.Fs;
        
        spec = PowerSpectrum_Wavelet(lfp,p.nCycles,p.dt,p.frange,'nfreqs',100,'space',p.space);
        
        feqInd = spec.freqs>=p.frange(1) & spec.freqs<=p.frange(end);
        
        h1 = subplot(6,1,1);
        plot(EEGTs(delayIdx)-EEGTs(delaystartIdx),EEG(delayIdx),'k');
        TITLE1 = sprintf('%s%d%s%d%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',sessDirs{j});
        TITLE2 = sprintf('%s%d','Trial-',m);
        title({TITLE1;TITLE2},'Interpreter','None');
        
        h2 = subplot(6,1,2);
        plot(EEGTs(delayIdx)-EEGTs(delaystartIdx),EEG_SWR(delayIdx),'k');
        title('SWR')
        ylim([-200 200])
        
        h3 = subplot(6,1,3);
        % plot raw LFP
        imagesc(spec.timestamps-spec.timestamps(1),log2(spec.freqs(feqInd)),spec.amp(:,feqInd)')
        title('Spectrogram')
        
        colorRange = spec.amp';
        range = [1 1.5];
        [~,mu,sig] = zscore(colorRange);
        crange = [min(mu)-range(1)*max(sig) max(mu)+range(2)*max(sig)];
        caxis(crange)
%         SpecColorRange(crange);
        colormap jet
        LogScale('y',2)
        ylabel('f (Hz)')
        axis xy
        %         hold on
        %         % plot bandpass ripple signal
        %         plot(EEGTs(delayIdx)-EEGTs(delayIdx(1)),EEG_SWR(delayIdx)/(max(EEG_SWR)/4)+10,'k');
               
        
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        binCount = length(rateBin1);
        ts_Def1_Off.spikeCount = zeros(validNum,binCount);
        ts_Def1_Off.cellLabel = zeros(validNum,binCount);
        ts_Def1_Off.spikeTime = cell(validNum,1);
        spikeCount_Off = zeros(validNum,1);
        
        for k = 1:validNum
            % get each spike time, change unit to msec from sec
            % ts unit: ms
            kk = validCellInd(k);
            tSp = Spike_Session.(sessDirs{j}){kk};
            if ~isempty(tSp)
                
                % delay definition 1 start from barrier
                ts_Delay1 = tSp(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                ts_Delay1Temp= ts_Delay1-delayTstart1(m);
                
                % time bin activate label
                % def 1, delay start at around barrier
                fireTemp = zeros(1,binCount);
                for n = 1:binCount
                    fireTemp(n) = sum(ts_Delay1Temp>rateBin1(n) & ts_Delay1Temp<rateBin2(n));
                end
                ts_Def1_Off.spikeCount(k,:) = fireTemp;
                ts_Def1_Off.cellLabel(k,:) = double(fireTemp>0);
                
                % each trial ts
                % add trials into combined ts
                if isempty(ts_Delay1)
                    ts_Delay1 = double.empty(0,1);
                end
                ts_Def1_Off.spikeTime{k} = ts_Delay1-delayTstart1(m);
                spikeCount_Off(k) = length(ts_Def1_Off.spikeTime{k});
            end
        end
        
        % calculate cofire and cross correlation for each pair
        
        cellTemp = ts_Def1_Off.cellLabel;
        tsTemp = ts_Def1_Off.spikeCount;
        cellPopMat = sum(cellTemp,1)/sum(rateLabel);
        ratePopMat = sum(tsTemp,1)/p.timeWindow;
        ratePopMat = gaussfilt(1:binCount,ratePopMat,2);
        
        h5=subplot(6,1,5);
        sigCofire = cellPopMat*sum(rateLabel)>=floor(p.cellThreshold*sum(rateLabel));
        %         plot(rateBin1+p.timeWindow/2,ratePopMat,'m')
        %         hold on
        plot(rateBin1+p.timeWindow/2,cellPopMat,'m')
        ylim([0 0.4])
%         xlabel('time (sec)')
        
        % define event by consecutive sigCofire labels
        sigTemp = diff([0,sigCofire]);
        eventStart = find(sigTemp == 1);
        eventEnd = find(sigTemp == -1);
        % consider boundary conditions
        if sigCofire(end) == 1 && sigCofire(end-1) == 0
            eventStart = eventStart(1:end-1);
        end
        if sigCofire(end) == 1 && sigCofire(end-1) == 1
            eventEnd = [eventEnd,length(sigTemp)+1];
        end
        eventEnd = eventEnd-1;
        
        if any((eventEnd-eventStart)<0)
            error('Population event detection is wrong')
        end
        
        eventStartTs = rateBin1(eventStart);
        eventEndTs = rateBin1(eventEnd);
        
        eventTspTemp = cell(length(eventStartTs),clusterNum);
        eventTsp2StartTemp = nan(length(eventStartTs),clusterNum);
            
        for kk = 1:length(eventStartTs)
            % get spikes within this event
            for k = 1:clusterNum
                
                tSp = Spike_Session.(sessDirs{j}){k};
                eventStartTemp3 = eventStartTs(kk)+delayTstart1(m);
                eventEndTemp3 = eventEndTs(kk)+delayTstart1(m);
                
                if rateLabel(k) == 1
                    spkTsTemp = tSp(tSp>=eventStartTemp3 & tSp<eventEndTemp3);
                    if ~isempty(spkTsTemp)
                        if sum(tSp>=spkTsTemp(1)-0.05& tSp<eventStartTemp3)<= 2
                            eventTspTemp{kk,k} = spkTsTemp-eventStartTemp3;
                        end
                    end
                    
                    % first spike timing to event start
                    if ~isempty(eventTspTemp{kk,k})
                        % first spike distance from start
                        eventTsp2StartTemp(kk,k) = eventTspTemp{kk,k}(1);
                    end
                end
            end
        end
        
        validInd = sum(~isnan(eventTsp2StartTemp),2)>= floor(1/5*sum(rateLabel));
        eventStartTs = eventStartTs(validInd);
        eventEndTs = eventEndTs(validInd);
        
        h4=subplot(6,1,4);
        hold on
        for kk = 1:length(eventStartTs)
            startTs = eventStartTs(kk);
            endTs = eventEndTs(kk)+p.timeWindow;
            y1 = -1*clusterNum-2;
            y2 = 10;
            h = patch([startTs endTs endTs startTs],[y1 y1 y2 y2],'y','FaceAlpha',0.4,'LineStyle','none');
        end        
        
        linkaxes([h1,h2,h3,h4,h5,h6],'x')
        xlim([0 30])
    end
    
    for j = 3
        
        % initiate the figure to plot spikes
        h = figure(2);
        h.Position = [100,100,1200,800];

        m = 7; % trial
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        % load turning directions and correctness
        tInfoFile = fullfile(mainDir,sessDirs{j}, 'trialInfo.mat');
        tInfo = load(tInfoFile);
        % -----------------------------------------------------------------
        % use the mode of delay is always 10 sec or 30 sec
        % rather than identify real time in delay zone
        % to make plot cleaner
        % I can change to delay time in delay zone later
        % Li Yuan, 19-Aug-2020, UCSD
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
        
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
            rateBin1 = 0:p.timeWindowIncrement:maxT-p.timeWindow-p.timeWindowIncrement; 
            rateBin2 = p.timeWindow:p.timeWindowIncrement:maxT-p.timeWindowIncrement;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
            rateBin1 = 0:p.timeWindowIncrement:maxT-p.timeWindow; 
            rateBin2 = p.timeWindow:p.timeWindowIncrement:maxT;
        else
            error('Delay time is wrong')
        end       
             
        % plot assemblies
        h6=subplot(6,1,6);
        for k = 1:patNum_on
            colorValue = [1-k*0.1,k*0.1,0];
            startTime = delayTstart1(m);
            endTime = delayTend1_2(m);
            binInd = binTime>=startTime & binTime<=endTime;

            AssmblStrength = AssmblStrength_on(k,binInd);
            plot(binTime(binInd)-startTime,log2(AssmblStrength),'Color',colorValue);
            hold on
        end
        ylim([-2 10])
        ylabel('log2 (Strength)')
        
        
        % plot EEG
        % delay area
        [startTimeDifference,delaystartIdx] = min(abs(EEGTs-delayTstart1(m)));
        [endTimeDifference,delayendIdx] = min(abs(EEGTs-delayTend1_2(m)));
        % accumulate delay Idx from each trial
        delayIdx = (delaystartIdx:delayendIdx);
        if any([startTimeDifference,endTimeDifference] > (EEGTs(2)-EEGTs(1)))
            error('EEG could not find matching position')
        end
        
        h1 = subplot(6,1,1);
        plot(EEGTs(delayIdx)-EEGTs(delaystartIdx),EEG(delayIdx),'k');
        TITLE1 = sprintf('%s%d%s%d%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',sessDirs{j});
        TITLE2 = sprintf('%s%d','Trial-',m);
        title({TITLE1;TITLE2},'Interpreter','None');
        
        h2 = subplot(6,1,2);
        plot(EEGTs(delayIdx)-EEGTs(delaystartIdx),EEG_SWR(delayIdx),'k');
        title('SWR')
        ylim([-200 200])
        
        lfp.data = EEG(delayIdx);
        lfp.timestamps = EEGTs(delayIdx);
        lfp.samplingRate = p.Fs;
        
        spec = PowerSpectrum_Wavelet(lfp,p.nCycles,p.dt,p.frange,'nfreqs',200,'space',p.space);
        
        feqInd = spec.freqs>=p.frange(1) & spec.freqs<=p.frange(end);
        
        h3 = subplot(6,1,3);
        % plot raw LFP
        imagesc(spec.timestamps-spec.timestamps(1),log2(spec.freqs(feqInd)),spec.amp(:,feqInd)')
        title('Spectrogram')
        caxis(crange)
        colormap jet

        LogScale('y',2)
        ylabel('f (Hz)')
        axis xy
        
        
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        binCount = length(rateBin1);
        ts_Def1_On.spikeCount = zeros(validNum,binCount);
        ts_Def1_On.cellLabel = zeros(validNum,binCount);
        ts_Def1_On.spikeTime = cell(validNum,1);
        spikeCount_On = zeros(validNum,1);
        
        for k = 1:validNum
            % get each spike time, change unit to msec from sec
            % ts unit: ms
            kk = validCellInd(k);
            tSp = Spike_Session.(sessDirs{j}){kk};
            if ~isempty(tSp)
                
                % delay definition 1 start from barrier
                ts_Delay1 = tSp(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                ts_Delay1Temp= ts_Delay1-delayTstart1(m);
                
                % time bin activate label
                % def 1, delay start at around barrier
                fireTemp = zeros(1,binCount);
                for n = 1:binCount
                    fireTemp(n) = sum(ts_Delay1Temp>rateBin1(n) & ts_Delay1Temp<rateBin2(n));
                end
                ts_Def1_On.spikeCount(k,:) = fireTemp;
                ts_Def1_On.cellLabel(k,:) = double(fireTemp>0);
                
                % each trial ts
                % add trials into combined ts
                if isempty(ts_Delay1)
                    ts_Delay1 = double.empty(0,1);
                end
                ts_Def1_On.spikeTime{k} = ts_Delay1-delayTstart1(m);
                spikeCount_On(k) = length(ts_Def1_On.spikeTime{k});
            end
        end
        
        % calculate cofire and cross correlation for each pair
   
        cellTemp = ts_Def1_On.cellLabel;
        tsTemp = ts_Def1_On.spikeCount;
        cellPopMat = sum(cellTemp,1)/30;
        ratePopMat = sum(tsTemp,1)/p.timeWindow;
        ratePopMat = gaussfilt(1:binCount,ratePopMat,2);

        h5 = subplot(6,1,5);
        sigCofire = cellPopMat>=1/5;
%         plot(rateBin1+p.timeWindow/2,ratePopMat,'m')
        plot(rateBin1+p.timeWindow/2,cellPopMat,'m')
        ylim([0 0.4])

                % define event by consecutive sigCofire labels
        sigTemp = diff([0,sigCofire]);
        eventStart = find(sigTemp == 1);
        eventEnd = find(sigTemp == -1);
        % consider boundary conditions
        if sigCofire(end) == 1 && sigCofire(end-1) == 0
            eventStart = eventStart(1:end-1);
        end
        if sigCofire(end) == 1 && sigCofire(end-1) == 1
            eventEnd = [eventEnd,length(sigTemp)+1];
        end
        eventEnd = eventEnd-1;
        
        if any((eventEnd-eventStart)<0)
            error('Population event detection is wrong')
        end
        
        eventStartTs = rateBin1(eventStart);
        eventEndTs = rateBin1(eventEnd);
        
        eventTspTemp = cell(length(eventStartTs),clusterNum);
        eventTsp2StartTemp = nan(length(eventStartTs),clusterNum);
            
        for kk = 1:length(eventStartTs)
            % get spikes within this event
            for k = 1:clusterNum
                
                tSp = Spike_Session.(sessDirs{j}){k};
                eventStartTemp3 = eventStartTs(kk)+delayTstart1(m);
                eventEndTemp3 = eventEndTs(kk)+delayTstart1(m);
                
                if rateLabel(k) == 1
                    spkTsTemp = tSp(tSp>=eventStartTemp3 & tSp<eventEndTemp3);
                    if ~isempty(spkTsTemp)
                        if sum(tSp>=spkTsTemp(1)-0.05 & tSp<eventStartTemp3)<= 2
                            eventTspTemp{kk,k} = spkTsTemp-eventStartTemp3;
                        end
                    end
                    
                    % first spike timing to event start
                    if ~isempty(eventTspTemp{kk,k})
                        % first spike distance from start
                        eventTsp2StartTemp(kk,k) = eventTspTemp{kk,k}(1);
                    end
                end
            end
        end
        
        validInd = sum(~isnan(eventTsp2StartTemp),2)>= floor(1/5*sum(rateLabel));
        eventStartTs = eventStartTs(validInd);
        eventEndTs = eventEndTs(validInd);

        h4 = subplot(6,1,4);
        hold on
        for kk = 1:length(eventStartTs)
            startTs = eventStartTs(kk);
            endTs = eventEndTs(kk)+p.timeWindow;
            y1 = -1*clusterNum-2;
            y2 = 10;
            h = patch([startTs endTs endTs startTs],[y1 y1 y2 y2],'y','FaceAlpha',0.4,'LineStyle','none');
        end
        
        linkaxes([h1,h2,h3,h4,h5,h6],'x')
        xlim([0 30])
    end
    
%     
%     [onSpike,onOrder] = sort(spikeCount_On);
%     onOrder = onOrder(onSpike>0);
%     offIdx = spikeCount_Off>0;
%     offIdx(onOrder) = 0;
%     off_CellIdx = find(offIdx == 1);
%     order = [onOrder;off_CellIdx];
    
    [~,cellorder] = sort(spikeCount_Off);
%     cellorder = 1:sum(rateLabel);
    figure(1)
    subplot(6,1,4);
    for kk = 1:length(cellorder)  
        
        k = cellorder(kk);
        if ~isempty(ts_Def1_Off.spikeTime{k})
            xPoints = [ts_Def1_Off.spikeTime{k}';ts_Def1_Off.spikeTime{k}'];
            yPoints = [-1*kk+zeros(size(ts_Def1_Off.spikeTime{k}'))-0.35;-1*kk+zeros(size(ts_Def1_Off.spikeTime{k}'))+0.35];
            plot(xPoints,yPoints,'k')
            hold on
        end
    end   %                     plot(spikeTemp1{k},-1*k,'b>');
    ylim([-25 0])
    
%     [~,order] = sort(spikeCount_On);
    figure(2)
    subplot(6,1,4)
    for kk = 1:length(cellorder)   
        k = cellorder(kk);
        if ~isempty(ts_Def1_On.spikeTime{k})
            xPoints = [ts_Def1_On.spikeTime{k}';ts_Def1_On.spikeTime{k}'];
            yPoints = [-1*kk+zeros(size(ts_Def1_On.spikeTime{k}'))-0.35;-1*kk+zeros(size(ts_Def1_On.spikeTime{k}'))+0.35];
            plot(xPoints,yPoints,'k')
            hold on
        end
    end   %                     plot(spikeTemp1{k},-1*k,'b>');
    ylim([-25 0])
    

%     patNum_off = CellAssembly_DelayLR.DelayOff.patNum; 
%     AssmblPtrnCellIDs_off = CellAssembly_DelayLR.DelayOff.AssmblPtrnCellIDs;
%     AssmblWght_onff = CellAssembly_DelayLR.DelayOff.AssmblWght;
    for m = 1:patNum_off
            
        figure(3)
        subplot(1,patNum_off,m)
        colorValue = [1-m*0.1,m*0.1,0];
        AssemblyWeight = AssmblWght_off(:,m);
        AssemblyWeight = AssemblyWeight(cellorder);
        [~,AssmblPtrnCellIDs] = ismember(AssmblPtrnCellIDs_off{m},cellorder);
        stem(AssemblyWeight,'k')
        hold on
        stem(AssmblPtrnCellIDs,AssemblyWeight(AssmblPtrnCellIDs),'Color',colorValue)
        set(gca,'XDir','reverse')
        xlim([0.5 sum(rateLabel)+0.5])
        view([90 -90])
        ylim([-1 1])
        title('off assemblies')
        
        sim_Temp = [];
        pattern_Off = AssmblWght_off(:,m);
        for mm = 1:patNum_on
            pattern_On = AssmblWght_on(:,mm);
            sim_Temp(mm) = dot(pattern_Off,pattern_On)/(norm(pattern_Off)*norm(pattern_On));   
        end
        [max_sim,maxInd] = max(sim_Temp);  
        TEXT = sprintf('%s%d%s%1.3f','Off:',maxInd,' ',max_sim);
        text(2,0,TEXT);
        
        
        figure(1)
        subplot(6,1,4);
        for kk = 1:length(AssmblPtrnCellIDs)
            
            k = AssmblPtrnCellIDs(kk);
            mm = AssmblPtrnCellIDs_off{m}(kk);
            if ~isempty(ts_Def1_Off.spikeTime{mm})
                xPoints = [ts_Def1_Off.spikeTime{mm}';ts_Def1_Off.spikeTime{mm}'];
                yPoints = [-1*k+zeros(size(ts_Def1_Off.spikeTime{mm}'))-0.35;-1*k+zeros(size(ts_Def1_Off.spikeTime{mm}'))+0.35];
                plot(xPoints,yPoints,'Color',colorValue)
                hold on
            end
        end
    end

    for m = 1:patNum_on
            
        figure(4)
        subplot(1,patNum_on,m)
        colorValue = [1-m*0.1,m*0.1,0];
        AssemblyWeight = AssmblWght_on(:,m);
        AssemblyWeight = AssemblyWeight(cellorder);
        [~,AssmblPtrnCellIDs] = ismember(AssmblPtrnCellIDs_on{m},cellorder);
        stem(AssemblyWeight,'k')
        hold on
        stem(AssmblPtrnCellIDs,AssemblyWeight(AssmblPtrnCellIDs),'Color',colorValue)
        set(gca,'XDir','reverse')
        xlim([0.5 sum(rateLabel)+0.5])
        view([90 -90])
        ylim([-1 1])
        title('on assemblies')
        
        sim_Temp = [];
        pattern_On = AssmblWght_on(:,m);
        for mm = 1:patNum_off
            pattern_Off = AssmblWght_off(:,mm);
            sim_Temp(mm) = dot(pattern_On,pattern_Off)/(norm(pattern_On)*norm(pattern_Off));   
        end
        [max_sim,maxInd] = max(sim_Temp);  
        TEXT = sprintf('%s%d%s%1.3f','Off:',maxInd,' ',max_sim);
        text(2,0,TEXT);
        
        figure(2)
        subplot(6,1,4);
        for kk = 1:length(AssmblPtrnCellIDs)
            
            k = AssmblPtrnCellIDs(kk);
            mm = AssmblPtrnCellIDs_on{m}(kk);
            if ~isempty(ts_Def1_On.spikeTime{mm})
                xPoints = [ts_Def1_On.spikeTime{mm}';ts_Def1_On.spikeTime{mm}'];
                yPoints = [-1*k+zeros(size(ts_Def1_On.spikeTime{mm}'))-0.35;-1*k+zeros(size(ts_Def1_On.spikeTime{mm}'))+0.35];
                plot(xPoints,yPoints,'Color',colorValue)
                hold on
            end
        end
    end
    
end
