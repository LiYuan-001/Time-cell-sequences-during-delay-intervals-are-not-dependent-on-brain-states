% Calculate cell firing rate in the delay area woth and wothout theta
% oscillations
% Li Yuan, Jan-26-2022, UCSD
function Fig8TreadmillDelayTheta_CellRate(inFile,AnalyzeSes)

p.savePlot = 0;
p.writeToFile = 1;

% set parameter
p.delta = [1 4];
p.theta = [6 10];

% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'
p.spikeFs = 30000;

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
    ThetaCellRate.rat = sessInfo(i).animal;
    ThetaCellRate.day = sessInfo(i).day;
    
 
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);

    for k = 1:clusterNum
        ThetaCellRate.tList{k} = TList{k}(1:end-2);
    end
            
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    
    %% get theta oscillation periods
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
    

    %% calculate rate for each block
    % 
    for j = 1:length(sessDirs)
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        % load path zone time
        pathZoneFile = fullfile(mainDir,sessDirs{j}, 'PathZone.mat');
        load(pathZoneFile);
        % load turning directions and correctness
        tInfoFile = fullfile(mainDir,sessDirs{j}, 'trialInfo.mat');
        tInfo = load(tInfoFile);
        
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
%         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
%         delayTend2 = Fig8DelayZonePos.delayPos2.endT;
        
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
        
        
        delay_SpNum = zeros(trialNum,clusterNum);        
        delay_SpNum_Theta = zeros(trialNum,clusterNum);
        delay_SpNum_NonTheta = zeros(trialNum,clusterNum);
            
        thetaTime = zeros(trialNum,1);
        non_thetaTime = zeros(trialNum,1);
        time_All = zeros(trialNum,1);
        
        % each trial    
        for m = 1:trialNum
            %--------------------------------------------------------------
            % delay area spike under theta
            thetaStartTemp1 = [];
            thetaStartTemp2 = [];
            thetaEndTemp1 = [];
            thetaEndTemp2 = [];
            % return + base spikes under theta
            % condition 1: theta start after location start
            thetaInd = tsStart >= delayTstart1(m) & tsStart <= delayTend1_2(m);
            if sum(thetaInd) > 0
                thetaStartTemp1 = tsStart(thetaInd);
                thetaEndTemp1 = tsEnd(thetaInd);
                if thetaEndTemp1(end) > delayTend1_2(m)
                    thetaEndTemp1(end) = delayTend1_2(m);
                end
            end
            % condition2: theta start before location start
            if any(tsStart < delayTstart1(m) & tsEnd > delayTstart1(m))
                thetaInd = (tsStart < delayTstart1(m) & tsEnd > delayTstart1(m));
                thetaStartTemp2 = delayTstart1(m);
                thetaEndTemp2 = tsEnd(thetaInd);                
                if thetaEndTemp2(end) > delayTend1_2(m)
                    thetaEndTemp2(end) = delayTend1_2(m);
                end
            end
            
            thetaStart = sort([thetaStartTemp1;thetaStartTemp2]);
            thetaEnd = sort([thetaEndTemp1;thetaEndTemp2]);
            
            time_All(m) = maxT;
            thetaTime(m) = sum(thetaEnd - thetaStart);
            non_thetaTime(m) = maxT - thetaTime(m);
            
            for k = 1:clusterNum
                tSp = Spike_Session.(sessDirs{j}){k};
                spkInd_all = tSp>=delayTstart1(m) & tSp<=delayTend1_2(m);
                
                spikeInd_Theta = [];
                if ~isempty(thetaStart) && sum(spkInd_all)
                    for kk = 1:length(thetaStart)
                        spikeInd = tSp>=thetaStart(kk) & tSp<=thetaEnd(kk);
                        spikeInd_Theta = [spikeInd_Theta;spikeInd];
                    end
                    delay_SpNum(m,k) = sum(spkInd_all);  
                    delay_SpNum_Theta(m,k) = sum(spikeInd_Theta);
                    delay_SpNum_NonTheta(m,k) = sum(spkInd_all) - sum(spikeInd_Theta);
                end
            end          
        end       
        
        ThetaCellRate.(sessDirs{j}).delay_SpNum = delay_SpNum;
        ThetaCellRate.(sessDirs{j}).delay_SpNum_Theta = delay_SpNum_Theta;
        ThetaCellRate.(sessDirs{j}).delay_SpNum_NonTheta = delay_SpNum_NonTheta;

        ThetaCellRate.(sessDirs{j}).time_All = time_All;
        ThetaCellRate.(sessDirs{j}).thetaTime = thetaTime;
        ThetaCellRate.(sessDirs{j}).non_thetaTime = non_thetaTime;
        
        rate_All = sum(delay_SpNum,1)./sum(time_All);
        rate_Theta = sum(delay_SpNum_Theta,1)./sum(thetaTime);
        rate_NonTheta = sum(delay_SpNum_NonTheta,1)./sum(non_thetaTime);
        
        ThetaCellRate.(sessDirs{j}).rate_All = rate_All;
        ThetaCellRate.(sessDirs{j}).rate_Theta = rate_Theta;
        ThetaCellRate.(sessDirs{j}).rate_NonTheta = rate_NonTheta;
    end
    
    if p.writeToFile == 1
        save(fullfile(savedir2,'ThetaCellRate.mat'), 'ThetaCellRate');
    end
    clear ThetaCellRate
    close all 
    fprintf('Finished position analysis for session %d\n',i);
end

end