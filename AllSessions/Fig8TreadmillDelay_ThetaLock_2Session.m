% define spike train as Sia's manuscript
% calculate the theta precession of the spike trians
% vs the behavioral time scale
% 2 slopes will be calculated
% a. slope by cell: slope of all spike trains being detected from this cell
% b. slope by train: slope of each spike train
function Fig8TreadmillDelay_ThetaLock_2Session(inFile,AnalyzeSes)

% set parameter
p.theta = [6 10];
% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'
p.spikeFs = 30000;

p.savePlot = 1;
p.writeToFile = 1;

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
MaxH = 0.4;
xDetail = 0:360;
MaxSpike = 3000;    % max spike number for raster
MaxT = 0.4;


% Read in input information
sessInfo = SessInfoImport(inFile);  

% parameters for import EEG from neuralynx / matlab utilities
FieldSelectionFlags = [1 0 0 0 1];
HeaderExtractionFlag = 1;
ExtractionMode = 1;
ExtractionModeVector = [];
    
for i = AnalyzeSes(1:end)
    close all
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\ThetaLock_2Blocks');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end
    
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
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    % initiate the data 
    ThetaLock_2Session.rat = sessInfo(i).animal;
    ThetaLock_2Session.day = sessInfo(i).day;
    
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    sessDirs = sessInfo(i).sessDirs;
    
    % combine matching sessions together
    if length(sessDirs) == 8
        sessName2 = {'on10','off10','on30','off30'};
    elseif contains(sessDirs{1},'on')
        sessName2 = {'on10','on30'};
    else
        sessName2 = {'off10','off30'};
    end
         
    % load spikes
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);    
    tList = cell(1,length(TList));
    for k = 1:clusterNum
        tList{k} = TList{k}(1:end-2);
    end    
    ThetaLock_2Session.tFile = tList;
    
    % load eeg file
    eegInd = sprintf('%s%d%s', 'CSC',sessInfo(i).EEGch2,'.ncs');
    cscFile = fullfile(mainDir,eegInd);
    % extract EEG by Nlx import
    [Timestamps,Samples,Header] = Nlx2MatCSC(cscFile, FieldSelectionFlags,HeaderExtractionFlag, ExtractionMode, ExtractionModeVector);
    p.t1 = Timestamps(1);    
    % reshape EEG samples
    EEG=reshape(Samples,length(Samples(:)),1);
    % get Frequency and ADBvolts for this channel
    [p.Fs, p.ADBVolts,p.Invert] = ReadHeaderEEG2(Header);
    % reasample timestamps;
    % unit: sec
    EEGTs = timeStampResample(Samples,Timestamps,p.Fs)./10^6;   

    % get EEG to uV
    EEG = EEG*p.ADBVolts*10^6*p.Invert;
    
    % get theta range and theta phase
    EEGt = fftbandpass(EEG,p.Fs,p.theta(1)-1,p.theta(1),p.theta(2),p.theta(2)+1);
    [phaseT,AmpT] = thetaPhase2(EEGt);
    
    % detect spike trains from all spikes in each cell
    
    for j = 1:length(sessDirs)
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        
         % use the mode of delay is always 10 sec or 30 sec
        % rather than identify real time in delay zone
        % to make plot cleaner
        % I can change to delay time in delay zone later
        % Li Yuan, 19-Aug-2020, UCSD
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
        trialNum = size(delayTstart1,2);
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
        else
            error('Delay time is wrong')
        end

        ts_DelayCombined = cell(clusterNum,1);        
        for k = 1:clusterNum
            % get each spike time
            tSp = Spike_Session.(sessDirs{j}){k};
            
            % get spike inside delay zone and calculate rate on time
            for m = 1:trialNum             
                % delay definition 2
                ts_Delay = tSp(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                ts_DelayCombined{k} = [ts_DelayCombined{k};ts_Delay];               
            end
        end
        
        % get session timestamp and extract EEG in this time
        for k = 1:clusterNum
            % get each spike time
            tSp = Spike_Session.(sessDirs{j}){k};
            
            lock = phaseLock(tSp,phaseT,EEGTs,p);            
            lock_Delay = phaseLock(ts_DelayCombined{k},phaseT,EEGTs,p);
            
            % write whole session thetalock to file
            ThetaLock_2Session.(sessDirs{j}).thetaLock_session.ppcT(k) = lock.ppcT;
            ThetaLock_2Session.(sessDirs{j}).thetaLock_session.rT(k) = lock.rT;
            if p.subSmaple
                ThetaLock_2Session.(sessDirs{j}).thetaLock_session.D_rT(k) = lock.D_rT;
            end
            ThetaLock_2Session.(sessDirs{j}).thetaLock_session.angT(k) = lock.angT;
            ThetaLock_2Session.(sessDirs{j}).thetaLock_session.pT(k) = lock.pT;
            ThetaLock_2Session.(sessDirs{j}).thetaLock_session.numT(k) = lock.numT;
            
            % write delay area thetalock to file
            ThetaLock_2Session.(sessDirs{j}).thetaLock_delay.ppcT(k) = lock_Delay.ppcT;
            ThetaLock_2Session.(sessDirs{j}).thetaLock_delay.rT(k) = lock_Delay.rT;
            if p.subSmaple
                ThetaLock_2Session.(sessDirs{j}).thetaLock_delay.D_rT(k) = lock_Delay.D_rT;
            end
            ThetaLock_2Session.(sessDirs{j}).thetaLock_delay.angT(k) = lock_Delay.angT;
            ThetaLock_2Session.(sessDirs{j}).thetaLock_delay.pT(k) = lock_Delay.pT;
            ThetaLock_2Session.(sessDirs{j}).thetaLock_delay.numT(k) = lock_Delay.numT;
            
            
            % normalized spike count with 30deg (or 45deg) bins
            if p.degBin==30
                x=15:30:345;
            else
                x=22.5:45:360;
            end
            nT = hist(lock.spikePhaseT,x);
            nT_Delay = hist(lock_Delay.spikePhaseT,x);
            
            if p.mode==1
                nT = nT./lock.numT;
                nT_Delay = nT_Delay./lock_Delay.numT;
            elseif p.mode==2
                disp('Output: Original spike counts')
            end
            
            h = figure(k);
            if j == 1
                h.Position = [50 50 1800 900];
            end
            
            % plot for whol session
            subplot(4,8,j)
            CircularPlot(lock.numT,lock.spikeRadPhaseT,lock.rT,lock.angT,p,MaxT)        
            TITLE1 = strcat(TList{k}(1:end-2), '-',sessDirs{j});
            TITLE2 = 'Theta Lock';
            title({TITLE1;TITLE2},'Interpreter','None');
            
            subplot(4,8,j+8)
            % Histgram-----------------------------------------------------------------       
            PlotThetaHist(nT,lock.spikePhaseT,lock.spikeRadPhaseT,xDetail,p,MinH,MaxH,lock.pT,MaxSpike,x)
            xlabel('Theta phase (deg)'); ylabel('Normalized spike counts')
            
            subplot(4,8,j+16)
            % plot for delay area theta lock
            CircularPlot(lock_Delay.numT,lock_Delay.spikeRadPhaseT,lock_Delay.rT,lock_Delay.angT,p,MaxT)        
%             TITLE1 = strcat(TList{k}(1:end-2), '-',SessDirs{j});
            TITLE2 = 'Delay ZONE';
            title({TITLE2},'Interpreter','None');
            
            subplot(4,8,j+24)
            % Histgram-----------------------------------------------------------------       
            PlotThetaHist(nT_Delay,lock_Delay.spikePhaseT,lock_Delay.spikeRadPhaseT,xDetail,p,MinH,MaxH,lock_Delay.pT,MaxSpike,x)
            xlabel('Theta phase (deg)'); ylabel('Normalized spike counts')
            
        end         
    end
    
    for sesGroup = 1:4
        ts_Fig8Combined = cell(clusterNum,1);
        ts_DelayCombined = cell(clusterNum,1);
        
        for j = [sesGroup,sesGroup+4]
            % load analyzed positions
            delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
            load(delayFile);
            
            % use the mode of delay is always 10 sec or 30 sec
            % rather than identify real time in delay zone
            % to make plot cleaner
            % I can change to delay time in delay zone later
            % Li Yuan, 19-Aug-2020, UCSD
            delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
            delayTend1 = Fig8DelayZonePos.delayPos1.endT;
            trialNum = size(delayTstart1,2);
            if contains(sessDirs{j},'10')
                maxT = 10;
                delayTend1_2 = delayTstart1+maxT;
            elseif contains(sessDirs{j},'30')
                maxT = 30;
                delayTend1_2 = delayTstart1+maxT;
            else
                error('Delay time is wrong')
            end
            
            for k = 1:clusterNum
                % get each spike time
                tSp = Spike_Session.(sessDirs{j}){k};
                ts_Fig8Combined{k} = [ts_Fig8Combined{k};tSp];
                % get spike inside delay zone and calculate rate on time
                for m = 1:trialNum
                    % delay definition 2
                    ts_Delay = tSp(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                    ts_DelayCombined{k} = [ts_DelayCombined{k};ts_Delay];
                end
            end
        end
        
        for k = 1:clusterNum
            % get each spike time
            tSp = ts_Fig8Combined{k};
            
            lock = phaseLock(tSp,phaseT,EEGTs,p);           
            lock_Delay = phaseLock(ts_DelayCombined{k},phaseT,EEGTs,p);
            
            % write whole session thetalock to file
            ThetaLock_2Session.(sessName2{sesGroup}).thetaLock_session.ppcT(k) = lock.ppcT;
            ThetaLock_2Session.(sessName2{sesGroup}).thetaLock_session.rT(k) = lock.rT;
            if p.subSmaple
                ThetaLock_2Session.(sessName2{sesGroup}).thetaLock_session.D_rT(k) = lock.D_rT;
            end
            ThetaLock_2Session.(sessName2{sesGroup}).thetaLock_session.angT(k) = lock.angT;
            ThetaLock_2Session.(sessName2{sesGroup}).thetaLock_session.pT(k) = lock.pT;
            ThetaLock_2Session.(sessName2{sesGroup}).thetaLock_session.numT(k) = lock.numT;
            
            % write delay area thetalock to file
            ThetaLock_2Session.(sessName2{sesGroup}).thetaLock_delay.ppcT(k) = lock_Delay.ppcT;
            ThetaLock_2Session.(sessName2{sesGroup}).thetaLock_delay.rT(k) = lock_Delay.rT;
            if p.subSmaple
                ThetaLock_2Session.(sessName2{sesGroup}).thetaLock_delay.D_rT(k) = lock_Delay.D_rT;
            end
            ThetaLock_2Session.(sessName2{sesGroup}).thetaLock_delay.angT(k) = lock_Delay.angT;
            ThetaLock_2Session.(sessName2{sesGroup}).thetaLock_delay.pT(k) = lock_Delay.pT;
            ThetaLock_2Session.(sessName2{sesGroup}).thetaLock_delay.numT(k) = lock_Delay.numT;
            
            
            % normalized spike count with 30deg (or 45deg) bins
            if p.degBin==30
                x=15:30:345;
            else
                x=22.5:45:360;
            end
            nT = hist(lock.spikePhaseT,x);
            nT_Delay = hist(lock_Delay.spikePhaseT,x);
            
            if p.mode==1
                nT = nT./lock.numT;
                nT_Delay = nT_Delay./lock_Delay.numT;
            elseif p.mode==2
                disp('Output: Original spike counts')
            end
            
            h = figure(k+clusterNum);
            if sesGroup == 1
                h.Position = [50 50 1200 900];
            end
            
            % plot for whol session
            subplot(4,4,sesGroup)
            CircularPlot(lock.numT,lock.spikeRadPhaseT,lock.rT,lock.angT,p,MaxT)        
            TITLE1 = strcat(TList{k}(1:end-2), '-',sessName2{sesGroup});
            TITLE2 = 'Theta Lock';
            title({TITLE1;TITLE2},'Interpreter','None');
            
            subplot(4,4,sesGroup+4)
            % Histgram-----------------------------------------------------------------       
            PlotThetaHist(nT,lock.spikePhaseT,lock.spikeRadPhaseT,xDetail,p,MinH,MaxH,lock.pT,MaxSpike,x)
            xlabel('Theta phase (deg)'); ylabel('Normalized spike counts')
            
            subplot(4,4,sesGroup+8)
            % plot for delay area theta lock
            CircularPlot(lock_Delay.numT,lock_Delay.spikeRadPhaseT,lock_Delay.rT,lock_Delay.angT,p,MaxT)        
%             TITLE1 = strcat(TList{k}(1:end-2), '-',SessDirs{j});
            TITLE2 = 'Delay ZONE';
            title({TITLE2},'Interpreter','None');
            
            subplot(4,4,sesGroup+12)
            % Histgram-----------------------------------------------------------------       
            PlotThetaHist(nT_Delay,lock_Delay.spikePhaseT,lock_Delay.spikeRadPhaseT,xDetail,p,MinH,MaxH,lock_Delay.pT,MaxSpike,x)
            xlabel('Theta phase (deg)'); ylabel('Normalized spike counts')
            
        end        
        
    end
    
    if p.savePlot
        for k = 1:clusterNum
            figure(k)
            figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',tList{k},'-ThetaLock');
            print(figName,'-dpng','-r300');
            figure(k+clusterNum )
            figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',tList{k},'-ThetaLock_2Blocks');
            print(figName,'-dpng','-r300');
        end
    end
    
    if p.writeToFile
        save(fullfile(savedir2,'ThetaLock_2Session.mat'), 'ThetaLock_2Session');
    end
    fprintf('Finished theta Lock analysis for session %d\n',i);
    clear ThetaLock_2Session
end
end