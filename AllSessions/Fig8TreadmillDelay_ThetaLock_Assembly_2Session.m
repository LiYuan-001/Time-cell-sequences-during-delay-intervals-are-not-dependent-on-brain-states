% theta lock for assemblies
% Li Yuan, Apr-26-2023, UCSD
function Fig8TreadmillDelay_ThetaLock_Assembly_2Session(inFile,AnalyzeSes)

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
    
    % load cluster file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;    
    tList = SpikeProp.tList(rateLabel);
    
    % load assembly file
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_DelayOnOff_LR-25ms.mat');
    load(assemblyFile);
    
    
    % initiate the data 
    ThetaLock_Assembly_2Session.rat = sessInfo(i).animal;
    ThetaLock_Assembly_2Session.day = sessInfo(i).day;
    
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
    patNum_on = CellAssembly_DelayLR.DelayOn.patNum; 
    clusterNum = patNum_on;
    
    AssmblPtrnCellIDs_on = CellAssembly_DelayLR.DelayOn.AssmblPtrnCellIDs;
    AssmblWght_on = CellAssembly_DelayLR.DelayOn.AssmblWght;
    AssmblStrength_on = CellAssembly_DelayLR.DelayOn.AssmblStrength;
    event_Time_on = CellAssembly_DelayLR.DelayOn.event_Time;
    event_strength_on = CellAssembly_DelayLR.DelayOn.event_strength;
    event_Num_on = CellAssembly_DelayLR.DelayOn.event_Num;
    
    ThetaLock_Assembly_2Session.AssemblyNum = patNum_on;
    ThetaLock_Assembly_2Session.AssmblPtrnCellIDs_on = AssmblPtrnCellIDs_on;
    
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
    
        % get session timestamp and extract EEG in this time
        for k = 1:clusterNum
            % get each spike time
            tSp = event_Time_on{k};
            
            % get spike inside delay zone and calculate rate on time
            for m = 1:trialNum             
                % delay definition 2
                ts_Delay = tSp(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                ts_DelayCombined{k} = [ts_DelayCombined{k},ts_Delay];               
            end
            
            lock_Delay = phaseLock(ts_DelayCombined{k},phaseT,EEGTs,p);
            
            % write delay area thetalock to file
            ThetaLock_Assembly_2Session.(sessDirs{j}).thetaLock_delay.ppcT(k) = lock_Delay.ppcT;
            ThetaLock_Assembly_2Session.(sessDirs{j}).thetaLock_delay.rT(k) = lock_Delay.rT;
            if p.subSmaple
                ThetaLock_Assembly_2Session.(sessDirs{j}).thetaLock_delay.D_rT(k) = lock_Delay.D_rT;
            end
            ThetaLock_Assembly_2Session.(sessDirs{j}).thetaLock_delay.angT(k) = lock_Delay.angT;
            ThetaLock_Assembly_2Session.(sessDirs{j}).thetaLock_delay.pT(k) = lock_Delay.pT;
            ThetaLock_Assembly_2Session.(sessDirs{j}).thetaLock_delay.numT(k) = lock_Delay.numT;
            
            
            % normalized spike count with 30deg (or 45deg) bins
            if p.degBin==30
                x=15:30:345;
            else
                x=22.5:45:360;
            end
            
            nT_Delay = hist(lock_Delay.spikePhaseT,x);
            
            if p.mode==1
                nT_Delay = nT_Delay./lock_Delay.numT;
            elseif p.mode==2
                disp('Output: Original spike counts')
            end
            
            h = figure(k);
            if j == 1
                h.Position = [50 50 1800 900];
            end
            
            subplot(4,9,9)
            TEXT = tList(AssmblPtrnCellIDs_on{k});
            text(0,0,TEXT,'Interpreter','None');
            axis tight
            axis off
        
            subplot(4,9,j)
            % plot for delay area theta lock
            CircularPlot(lock_Delay.numT,lock_Delay.spikeRadPhaseT,lock_Delay.rT,lock_Delay.angT,p,MaxT)        
%             TITLE1 = strcat(TList{k}(1:end-2), '-',SessDirs{j});
            TITLE1 = sprintf('%s%d%s%s','On_Assembly: ',k,' ',sessDirs{j});
            TITLE2 = 'Delay ZONE Theta Lock';
            title({TITLE1;TITLE2},'Interpreter','None');
            
            subplot(4,9,j+9)
            % Histgram-----------------------------------------------------------------       
            PlotThetaHist(nT_Delay,lock_Delay.spikePhaseT,lock_Delay.spikeRadPhaseT,xDetail,p,MinH,MaxH,lock_Delay.pT,MaxSpike,x)
            xlabel('Theta phase (deg)'); ylabel('Normalized spike counts')
            
        end         
    end
    
    for sesGroup = 1:4
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
                tSp = event_Time_on{k};
                % get spike inside delay zone and calculate rate on time
                for m = 1:trialNum
                    % delay definition 2
                    ts_Delay = tSp(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                    ts_DelayCombined{k} = [ts_DelayCombined{k},ts_Delay];
                end
            end
        end
        
        for k = 1:clusterNum
            % get each spike time     
            lock_Delay = phaseLock(ts_DelayCombined{k},phaseT,EEGTs,p);
            % write delay area thetalock to file
            ThetaLock_Assembly_2Session.(sessName2{sesGroup}).thetaLock_delay.ppcT(k) = lock_Delay.ppcT;
            ThetaLock_Assembly_2Session.(sessName2{sesGroup}).thetaLock_delay.rT(k) = lock_Delay.rT;
            if p.subSmaple
                ThetaLock_Assembly_2Session.(sessName2{sesGroup}).thetaLock_delay.D_rT(k) = lock_Delay.D_rT;
            end
            ThetaLock_Assembly_2Session.(sessName2{sesGroup}).thetaLock_delay.angT(k) = lock_Delay.angT;
            ThetaLock_Assembly_2Session.(sessName2{sesGroup}).thetaLock_delay.pT(k) = lock_Delay.pT;
            ThetaLock_Assembly_2Session.(sessName2{sesGroup}).thetaLock_delay.numT(k) = lock_Delay.numT;
            
            
            % normalized spike count with 30deg (or 45deg) bins
            if p.degBin==30
                x=15:30:345;
            else
                x=22.5:45:360;
            end
            nT_Delay = hist(lock_Delay.spikePhaseT,x);
            
            if p.mode==1
                nT_Delay = nT_Delay./lock_Delay.numT;
            elseif p.mode==2
                disp('Output: Original spike counts')
            end
            
            h = figure(k);
            
            subplot(4,4,sesGroup+8)
            % plot for delay area theta lock
            CircularPlot(lock_Delay.numT,lock_Delay.spikeRadPhaseT,lock_Delay.rT,lock_Delay.angT,p,MaxT)        
%             TITLE1 = strcat(TList{k}(1:end-2), '-',SessDirs{j});
            TITLE2 = sessName2{sesGroup};
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
            figName = sprintf('%s%s%d%s%d%s%d%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-OnAssembly',k,'-ThetaLock');
            print(figName,'-dpng','-r300');
        end
    end
    
    if p.writeToFile
        save(fullfile(savedir2,'ThetaLock_Assembly_2Session.mat'), 'ThetaLock_Assembly_2Session');
    end
    fprintf('Finished theta Lock analysis for session %d\n',i);
    clear ThetaLock_Assembly_2Session
end
end