function Fig8Treadmill_Delay_PairwiseCorr(inFile,AnalyzeSes)

close all
p.savePlot = 0;
p.writeToFile = 1;

p.binWidth = 20./10^3; % unit sec
% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'
p.spikeFs = 30000;

% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    sessDirs = sessInfo(i).sessDirs;
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load population file
    reactFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_WholeSes.mat');
    load(reactFile);
    % load reactivation delay specific file
    reactDelayFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_Delay.mat');
    load(reactDelayFile);
    
    % initiate the data
    CellPairwiseCorr.rat = sessInfo(i).animal;
    CellPairwiseCorr.day = sessInfo(i).day;
    CellPairwiseCorr.binWidth = p.binWidth;
    
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    
    for j = 1:length(sessDirs)        
        % delay area
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        % load path zone time
        pathZoneFile = fullfile(mainDir,sessDirs{j}, 'PathZone.mat');
        load(pathZoneFile);
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
        trialNum = size(delayTstart1,2);
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
            rateBin1 = 0:p.binWidth:maxT-p.binWidth;
            rateBin2 = p.binWidth:p.binWidth:maxT;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
            rateBin1 = 0:p.binWidth:maxT-p.binWidth;
            rateBin2 = p.binWidth:p.binWidth:maxT;
        else
            error('Delay time is wrong')
        end
        
        binCount = length(rateBin1);
        ts_Def1.spikeCount = zeros(clusterNum,trialNum,binCount);
        
        for k = 1:clusterNum
            % get each spike time, change unit to msec from sec
            % ts unit: ms
            
            tSp = Spike_Session.(sessDirs{j}){k};
            if ~isempty(tSp)
                % only consider pupative pyramidal cells
                if rateLabel(k) == 1
                    for m = 1:trialNum
                        
                        % delay definition 1 start from barrier
                        ts_Delay1 = tSp(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                        ts_Delay1Temp= ts_Delay1-delayTstart1(m);
                        
                        % time bin activate label
                        % def 1, delay start at afloor barrier
                        fireTemp = zeros(1,binCount);
                        for n = 1:binCount
                            fireTemp(n) = sum(ts_Delay1Temp>rateBin1(n) & ts_Delay1Temp<rateBin2(n));
                        end
                        % number of spikes in each time bin
                        ts_Def1.spikeCount(k,m,:) = fireTemp;
                        
                    end
                end
            end
        end
        
        pairSize = nchoosek(clusterNum,2);
        corrTemp = zeros(pairSize,trialNum);
        pairCount = 0;
        for kk = 1:clusterNum-1
            for mm = kk+1:clusterNum
                pairCount = pairCount + 1;
                for m = 1:trialNum
                    seq1 = squeeze(ts_Def1.spikeCount(kk,m,:));
                    seq2 = squeeze(ts_Def1.spikeCount(mm,m,:));
                    corrTemp(pairCount,m) = corr(seq1,seq2,'Rows','pairwise');
                end
            end
        end
        
        CellPairwiseCorr.(sessDirs{j}).delay_pairWiseCorr = corrTemp;
        
        
        %% delay area
        rewardTime = [PathZone.posStartT.Reward,PathZone.posEndT.Reward]; % reward
        corrTemp = zeros(pairSize,trialNum);
        for m = 1:trialNum
            pairCount = 0;
            rateBin1 = rewardTime(m,1):p.binWidth:rewardTime(m,2)-p.binWidth;
            rateBin2 = rewardTime(m,1)+p.binWidth:p.binWidth:rewardTime(m,2);
            binCount = length(rateBin1);
            for kk = 1:clusterNum-1
                for mm = kk+1:clusterNum
                    pairCount = pairCount + 1;
                    % delay definition 1 start from barrier
                    tSp1 = Spike_Session.(sessDirs{j}){kk};
                    tSp2 = Spike_Session.(sessDirs{j}){mm};
                    % time bin activate label
                    % def 1, delay start at afloor barrier
                    fireTemp1 = zeros(1,binCount);
                    fireTemp2 = zeros(1,binCount);
                    for n = 1:binCount
                        fireTemp1(n) = sum(tSp1>rateBin1(n) & tSp1<rateBin2(n));
                        fireTemp2(n) = sum(tSp2>rateBin1(n) & tSp2<rateBin2(n));
                    end
                    % number of spikes in each time bin
                    corrTemp(pairCount,m) = corr(fireTemp1',fireTemp2','Rows','pairwise');
                end
            end
        end
        CellPairwiseCorr.(sessDirs{j}).reward_pairWiseCorr = corrTemp;
    end
    
    if p.writeToFile == 1
        save(fullfile(savedir2,'CellPairwiseCorr.mat'), 'CellPairwiseCorr');
    end
    clear CellPairwiseCorr
    fprintf('Finished analysis for session %d\n',i)
end
end
