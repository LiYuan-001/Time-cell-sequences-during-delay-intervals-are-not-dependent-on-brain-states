% example for population activity in delay area for manuscript
% example for popilation activity in delay and reward area
close all

p.savePlot = 0;
p.writeToFile = 0;

% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'

% Read in input information
sessInfo = SessInfoImport('V:\LiYuan\Codes\Fig8MazeTreadmill_V2\Fig8Treadmill_OnOff.xlsx');

% EEG related
% parameters for import EEG from neuralynx / matlab utilities
FieldSelectionFlags = [1 0 0 0 1];
HeaderExtractionFlag = 1;
ExtractionMode = 1;
ExtractionModeVector = [];

for i = 5

%     if p.savePlot
%         % directory for plot figures
%         % generate a folder for each rat eah day under the current folder
%         savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures',sessInfo(i).animal,'-day',sessInfo(i).day,'\Cell pair cofire');
%         if ~exist(savedir, 'dir')
%             mkdir(savedir);
%         end
%     end

    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);

    % load assembly file
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_WholeSes-25ms.mat');
    load(assemblyFile);
    patNum = CellAssembly_WholeSes.patNum ;
    AssmblPtrnCellIDs = CellAssembly_WholeSes.AssmblPtrnCellIDs;
    AssmblWght = CellAssembly_WholeSes.AssmblWght;
    AssmblStrength = CellAssembly_WholeSes.AssmblStrength;
    event_Time = CellAssembly_WholeSes.event_Time;
    event_strength = CellAssembly_WholeSes.event_strength;
    event_Num = CellAssembly_WholeSes.event_Num;
    bin_Time = CellAssembly_WholeSes.binTime;
    
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
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
    rateLabelID = find(rateLabel==1);
    
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    % import EEG
    % load eeg file
    EEGch = 1;
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
    [evtOnIdx.R,evtOffIdx.R,eegLabel.R] = eventFinder_SWR(EEG,EEG_SWR,p.Fs,3,0.05,0.020,0);
    swr_Ts = EEGTs(eegLabel.R==1);
    
    % block 7 trial 6 for treadmill on
    % block 8 trial 1 for treadmill off 
    
        
    for j = 6 % session
        m = 6; % trial
        % initiate the figure to plot spikes
        h = figure(1);
        h.Position = [100,100,1200,800];
        
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
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
        else
            error('Delay time is wrong')
        end
        
        % plot EEG
        % delay area
        [startTimeDifference,delaystartIdx] = min(abs(EEGTs-delayTstart1(m)));
        [endTimeDifference,delayendIdx] = min(abs(EEGTs-delayTend1_2(m)));
        % accumulate delay Idx from each trial
        delayIdx = (delaystartIdx:delayendIdx);
        if any([startTimeDifference,endTimeDifference] > (EEGTs(2)-EEGTs(1)))
            error('EEG could not find matching position')
        end
        
        h1 = subplot(5,1,1);
        plot(EEGTs(delayIdx)-EEGTs(delaystartIdx),EEG(delayIdx),'k');
        TITLE1 = sprintf('%s%d%s%d%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',sessDirs{j});
        TITLE2 = sprintf('%s%d','Trial-',m);
        title({TITLE1;TITLE2},'Interpreter','None');
        
        h2 = subplot(5,1,2);
        plot(EEGTs(delayIdx)-EEGTs(delaystartIdx),EEG_SWR(delayIdx),'k');
        title('SWR')
        ylim([-200 200])
        
        bin_Time_Idx = bin_Time>=Fig8DelayZonePos.delayPos1.startT(m) & bin_Time<delayTend1_2(m);
        bin_Time_Temp = bin_Time(bin_Time_Idx);
        bin_Time_SWR = zeros(1,sum(bin_Time_Idx));
        bin_Time_Temp_2 = bin_Time_Temp-bin_Time_Temp(1);
        
        for mm = 1:length(bin_Time_Temp)-1
            if sum(swr_Ts>=bin_Time_Temp(mm) & swr_Ts<bin_Time_Temp(mm+1)) > 0
                bin_Time_SWR(mm) = 1;
            end
        end
        hold on
        plot(bin_Time_Temp_2+0.1,bin_Time_SWR*30-150,'k')
        linkaxes([h1,h2],'x')
        
        h3 = subplot(5,1,[3,4]);
        for mm = rateLabelID
            tSp = Spike_Session.(sessDirs{j}){mm};
            ts_Delay1 = tSp(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
            ts_Delay1Temp = ts_Delay1-delayTstart1(m);
            
            xPoints = [ts_Delay1Temp';ts_Delay1Temp'];
            yPoints = [-1*mm+zeros(size(ts_Delay1Temp'))-0.35;-1*mm+zeros(size(ts_Delay1Temp'))+0.35];
            plot(xPoints,yPoints,'Color',[0.3,0.3,0.3])
            hold on
        end
            
        for k = 1:patNum
            plotColor = [1-1/patNum*k,0.1,1/patNum*k];
            cellInd = rateLabelID(AssmblPtrnCellIDs{k});
            for mm = cellInd
                tSp = Spike_Session.(sessDirs{j}){mm};
                ts_Delay1 = tSp(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                ts_Delay1Temp = ts_Delay1-delayTstart1(m);
                
                xPoints = [ts_Delay1Temp';ts_Delay1Temp'];
                yPoints = [-1*mm+zeros(size(ts_Delay1Temp'))-0.35;-1*mm+zeros(size(ts_Delay1Temp'))+0.35];
                plot(xPoints,yPoints,'Color',plotColor)
                hold on
            end
        end
        
        h4 = subplot(5,1,5);
        for k = 1:patNum
            plotColor = [1-1/patNum*k,0.1,1/patNum*k];
            event_tSp = event_Time{k};
            strength = AssmblStrength(k,:);
            bin_Delay1 = bin_Time(bin_Time>delayTstart1(m) & bin_Time<delayTend1_2(m));
            event_ts_Delay1Temp = bin_Delay1-delayTstart1(m);
            plot(event_ts_Delay1Temp,strength(bin_Time>delayTstart1(m) & bin_Time<delayTend1_2(m)),'Color',plotColor);
            hold on
        end
        linkaxes([h1,h2,h3,h4],'x')
    end
    
    for j = 5 % session
        m = 2; % trial
        % initiate the figure to plot spikes
        h = figure(2);
        h.Position = [100,100,1200,800];
        
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
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
        else
            error('Delay time is wrong')
        end
        
        % plot EEG
        % delay area
        [startTimeDifference,delaystartIdx] = min(abs(EEGTs-delayTstart1(m)));
        [endTimeDifference,delayendIdx] = min(abs(EEGTs-delayTend1_2(m)));
        % accumulate delay Idx from each trial
        delayIdx = (delaystartIdx:delayendIdx);
        if any([startTimeDifference,endTimeDifference] > (EEGTs(2)-EEGTs(1)))
            error('EEG could not find matching position')
        end
        
        h1 = subplot(5,1,1);
        plot(EEGTs(delayIdx)-EEGTs(delaystartIdx),EEG(delayIdx),'k');
        TITLE1 = sprintf('%s%d%s%d%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',sessDirs{j});
        TITLE2 = sprintf('%s%d','Trial-',m);
        title({TITLE1;TITLE2},'Interpreter','None');
        
        h2 = subplot(5,1,2);
        plot(EEGTs(delayIdx)-EEGTs(delaystartIdx),EEG_SWR(delayIdx),'k');
        title('SWR')
        ylim([-200 200])
        
        bin_Time_Idx = bin_Time>=Fig8DelayZonePos.delayPos1.startT(m) & bin_Time<delayTend1_2(m);
        bin_Time_Temp = bin_Time(bin_Time_Idx);
        bin_Time_SWR = zeros(1,sum(bin_Time_Idx));
        bin_Time_Temp_2 = bin_Time_Temp-bin_Time_Temp(1);
        
        for mm = 1:length(bin_Time_Temp)-1
            if sum(swr_Ts>=bin_Time_Temp(mm) & swr_Ts<bin_Time_Temp(mm+1)) > 0
                bin_Time_SWR(mm) = 1;
            end
        end
        hold on
        plot(bin_Time_Temp_2+0.1,bin_Time_SWR*30-150,'k')

        
        h3 = subplot(5,1,[3,4]);
        for mm = rateLabelID
            tSp = Spike_Session.(sessDirs{j}){mm};
            ts_Delay1 = tSp(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
            ts_Delay1Temp = ts_Delay1-delayTstart1(m);
            
            xPoints = [ts_Delay1Temp';ts_Delay1Temp'];
            yPoints = [-1*mm+zeros(size(ts_Delay1Temp'))-0.35;-1*mm+zeros(size(ts_Delay1Temp'))+0.35];
            plot(xPoints,yPoints,'Color',[0.3,0.3,0.3])
            hold on
        end
            
        for k = 1:patNum
            plotColor = [1-1/patNum*k,0.1,1/patNum*k];
            cellInd = rateLabelID(AssmblPtrnCellIDs{k});
            for mm = cellInd
                tSp = Spike_Session.(sessDirs{j}){mm};
                ts_Delay1 = tSp(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                ts_Delay1Temp = ts_Delay1-delayTstart1(m);
                
                xPoints = [ts_Delay1Temp';ts_Delay1Temp'];
                yPoints = [-1*mm+zeros(size(ts_Delay1Temp'))-0.35;-1*mm+zeros(size(ts_Delay1Temp'))+0.35];
                plot(xPoints,yPoints,'Color',plotColor)
                hold on
            end
        end
        
        h4 = subplot(5,1,5);
        for k = 1:patNum
            plotColor = [1-1/patNum*k,0.1,1/patNum*k];
            event_tSp = event_Time{k};
            strength = AssmblStrength(k,:);
            bin_Delay1 = bin_Time(bin_Time>delayTstart1(m) & bin_Time<delayTend1_2(m));
            event_ts_Delay1Temp = bin_Delay1-delayTstart1(m);
            plot(event_ts_Delay1Temp,strength(bin_Time>delayTstart1(m) & bin_Time<delayTend1_2(m)),'Color',plotColor);
            hold on
        end
        linkaxes([h1,h2,h3,h4],'x')
    end
    
    
    %% 
    % delay and reward area compare
    % delay area
    
    % load reactivation delay specific file
    reactDelayFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_Delay.mat');
    load(reactDelayFile);
    
    for j = 8 % session
        eventStartTs = DelayPopFire_Delay.(sessDirs{j}).eventStartTs;
        eventEndTs = DelayPopFire_Delay.(sessDirs{j}).eventEndTs;
        
        for m = 4
            
            ts_Def1_off.spikeTime = cell(clusterNum,1);
            [startTimeDifference,delaystartIdx] = min(abs(EEGTs-(eventStartTs(m)-0.3)));
            [endTimeDifference,delayendIdx] = min(abs(EEGTs-(eventEndTs(m)+0.3)));
            % accumulate delay Idx from each trial
            delayIdx = (delaystartIdx:delayendIdx);
            if any([startTimeDifference,endTimeDifference] > (EEGTs(2)-EEGTs(1)))
                error('EEG could not find matching position')
            end
            
            lfp.data = EEG(delayIdx);
            lfp.timestamps = EEGTs(delayIdx);
            lfp.samplingRate = p.Fs;
            
            figure
            h1 = subplot(2,1,1);
            plot(EEGTs(delaystartIdx:delayendIdx)-(eventStartTs(m)),EEG_SWR(delayIdx),'k');
            title('SWR')
            ylim([-200 200])
            title('reward pop event')
            
            h2 = subplot(2,1,2);
            startTs = -0.1;
            endTs = eventEndTs(m)-eventStartTs(m)+0.1;
            y1 = -1*clusterNum-2;
            y2 = 10;
            h = patch([startTs endTs endTs startTs],[y1 y1 y2 y2],'y','FaceAlpha',0.4,'LineStyle','none');
            hold on
            
            for k = 1:clusterNum
                % get each spike time, change unit to msec from sec
                % ts unit: ms
                tSp = Spike_Session.(sessDirs{j}){k};
                if rateLabel(k) == 1
                    if ~isempty(tSp)
                        
                        % delay definition 1 start from barrier
                        ts_Delay1 = tSp(tSp>(eventStartTs(m)-0.3) & tSp<(eventEndTs(m)+0.3));
                        % each trial ts
                        % add trials into combined ts
                        if isempty(ts_Delay1)
                            ts_Delay1 = double.empty(0,1);
                        end
                        ts_Def1_off.spikeTime{k} = ts_Delay1-(eventStartTs(m));
                    end
                    
                end
            end
            order = 1:42;
            subplot(2,1,2)
            for kk = 1:length(order)
                k = order(kk);
                if ~isempty(ts_Def1_off.spikeTime{k})
                    xPoints = [ts_Def1_off.spikeTime{k}';ts_Def1_off.spikeTime{k}'];
                    yPoints = [-1*kk+zeros(size(ts_Def1_off.spikeTime{k}'))-0.35;-1*kk+zeros(size(ts_Def1_off.spikeTime{k}'))+0.35];
                    plot(xPoints,yPoints,'k')
                    hold on
                end
            end   %                     plot(spikeTemp1{k},-1*k,'b>');
            ylim([-45 0])
            linkaxes([h1,h2],'x')
        end
    end
    
    
    
    % load reactivation file
    reactFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_WholeSes.mat');
    load(reactFile);
    % reward area
    for j = 4
        
        validInd = DelayPopFire_WholeSes.(sessDirs{j}).eventLoc == 5;
        eventStartTs = DelayPopFire_WholeSes.(sessDirs{j}).eventStartTs(validInd);
        eventEndTs = DelayPopFire_WholeSes.(sessDirs{j}).eventEndTs(validInd);
        
        for m = 3
            
            ts_Def1_reward.spikeTime = cell(clusterNum,1);
            [startTimeDifference,delaystartIdx] = min(abs(EEGTs-(eventStartTs(m)-0.3)));
            [endTimeDifference,delayendIdx] = min(abs(EEGTs-(eventEndTs(m)+0.3)));
            % accumulate delay Idx from each trial
            delayIdx = (delaystartIdx:delayendIdx);
            if any([startTimeDifference,endTimeDifference] > (EEGTs(2)-EEGTs(1)))
                error('EEG could not find matching position')
            end
            
            lfp.data = EEG(delayIdx);
            lfp.timestamps = EEGTs(delayIdx);
            lfp.samplingRate = p.Fs;
            
            figure
            h1 = subplot(2,1,1);
            plot((EEGTs(delaystartIdx:delayendIdx)-eventStartTs(m)),EEG_SWR(delayIdx),'k');
            title('SWR')
            ylim([-200 200])
            title('reward pop event')
            
            h2 = subplot(2,1,2);
            startTs = 0;
            endTs = eventEndTs(m)-eventStartTs(m)+0.1;
            y1 = -1*clusterNum-2;
            y2 = 10;
            h = patch([startTs endTs endTs startTs],[y1 y1 y2 y2],'y','FaceAlpha',0.4,'LineStyle','none');
            hold on
            
            for k = 1:clusterNum
                % get each spike time, change unit to msec from sec
                % ts unit: ms
                tSp = Spike_Session.(sessDirs{j}){k};
                if rateLabel(k) == 1
                    if ~isempty(tSp)
                        
                        % delay definition 1 start from barrier
                        ts_Delay1 = tSp(tSp>(eventStartTs(m)-0.3) & tSp<(eventEndTs(m)+0.3));
                        % each trial ts
                        % add trials into combined ts
                        if isempty(ts_Delay1)
                            ts_Delay1 = double.empty(0,1);
                        end
                        ts_Def1_reward.spikeTime{k} = ts_Delay1-(eventStartTs(m));
                    end
                    
                end
            end
            order = 1:42;
            subplot(2,1,2)
            for kk = 1:length(order)
                k = order(kk);
                if ~isempty(ts_Def1_reward.spikeTime{k})
                    xPoints = [ts_Def1_reward.spikeTime{k}';ts_Def1_reward.spikeTime{k}'];
                    yPoints = [-1*kk+zeros(size(ts_Def1_reward.spikeTime{k}'))-0.35;-1*kk+zeros(size(ts_Def1_reward.spikeTime{k}'))+0.35];
                    plot(xPoints,yPoints,'k')
                    hold on
                end
            end   %                     plot(spikeTemp1{k},-1*k,'b>');
            ylim([-45 0])
            linkaxes([h1,h2],'x')
        end
    end
end
