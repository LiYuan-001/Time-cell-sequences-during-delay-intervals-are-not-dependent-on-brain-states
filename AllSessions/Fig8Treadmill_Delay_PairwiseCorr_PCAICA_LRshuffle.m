function Fig8Treadmill_Delay_PairwiseCorr_PCAICA_LRshuffle(inFile,AnalyzeSes)

close all
p.savePlot = 0;
p.writeToFile = 1;

p.binWidth = 25./10^3; % unit sec
% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'
p.spikeFs = 30000;

opts.Patterns.method = 'ICA';
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.number_of_iterations = 500;

% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    sessDirs = sessInfo(i).sessDirs;
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
%     % load population file
%     reactFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_WholeSes.mat');
%     load(reactFile);
    % load reactivation delay specific file
%     reactDelayFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_Delay.mat');
%     load(reactDelayFile);
%     
    % initiate the data
    CellAseembly.rat = sessInfo(i).animal;
    CellAseembly.day = sessInfo(i).day;
    CellAseembly.binWidth = p.binWidth;
    
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    
    % zscore spike train for whole session first
    sessStart = zeros(1,length(sessDirs));
    sessEnd = zeros(1,length(sessDirs));
    tsAllblocks = cell(clusterNum,1);
    for j = 1:length(sessDirs) 
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        sessStart(j) = pathData.t(1);
        sessEnd(j) = pathData.t(end);
        for k = 1:clusterNum
            tSpTemp = Spike_Session.(sessDirs{j}){k};
            tsAllblocks{k} = [tsAllblocks{k};tSpTemp];
        end
        
    end
    startT = min(sessStart);
    endT = max(sessEnd);
    rateBin1 = startT:p.binWidth:endT-p.binWidth;
    rateBin2 = startT+p.binWidth:p.binWidth:endT;
    binCount = length(rateBin1);
    spkTrain = zeros(clusterNum,binCount);
    for k = 1:clusterNum
        if rateLabel(k) == 1
            tSpTemp = sort(tsAllblocks{k});
            fireTemp = zeros(1,binCount);
            for n = 1:binCount
                fireTemp(n) = sum(tSpTemp>rateBin1(n) & tSpTemp<=rateBin2(n));
            end
%             spkTrain(k,:) = zscore(fireTemp);
            spkTrain(k,:) = fireTemp;
        end
    end
    
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
        if strcmp(sessDirs{j}(end-3:end-2),'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
        elseif strcmp(sessDirs{j}(end-3:end-2),'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
        else
            error('Delay time is wrong')
        end
        
        seqTotal = [];
        turnLabel = zeros(trialNum,1);
        for m = 1:trialNum
            delayBinIdx = rateBin1>=delayTstart1(m) & rateBin2<=delayTend1_2(m);
            seqTotal = cat(3,seqTotal,spkTrain(rateLabel,delayBinIdx));
            if strcmp(tInfo.direction{m},'R')
                turnLabel(m) = 1;
            end
        end
        
        % calculate original assembly and strength
        seqTotal_L = [];
        seqTotal_R = [];
        for m = 1:trialNum
            delayBinIdx = rateBin1>=delayTstart1(m) & rateBin2<=delayTend1_2(m);
            if strcmp(tInfo.direction{m},'L')
                seqTotal_L = [seqTotal_L,spkTrain(rateLabel,delayBinIdx)];
            elseif strcmp(tInfo.direction{m},'R')
                seqTotal_R = [seqTotal_R,spkTrain(rateLabel,delayBinIdx)];
            else
                error('Direction label error')
            end         
        end
        
        %%
        z_seqTotal_L = zscore(seqTotal_L')';
        z_seqTotal_R = zscore(seqTotal_R')';
        
        % first : get asseembly on Left and Right
        % second: get assembly strength L-L,L-R,R-R,R-L
        % 
        % get assemblies
        % Left delay
        AssemblyTemplates = assembly_patterns(seqTotal_L,opts);
        AssemblyTemplates2 = AssemblyTemplates;
        % define if it is a cell assembly
        cellID = 1:sum(rateLabel);
        nPatterns = size(AssemblyTemplates,2);
        AssmblPtrnCellIDs_L=cell(1,nPatterns);
        validInd = zeros(1,nPatterns);
        assmblPtrnWgtsTemp = [];
        for pIdx = 1:nPatterns
            assmblPtrnWgtsTemp(:,pIdx) = AssemblyTemplates(:,pIdx)./norm(AssemblyTemplates(:,pIdx));
            if max(assmblPtrnWgtsTemp(:,pIdx)) < abs(min(assmblPtrnWgtsTemp(:,pIdx)))
                assmblPtrnWgtsTemp(:,pIdx) = -assmblPtrnWgtsTemp(:,pIdx);
                AssemblyTemplates2(:,pIdx) = -AssemblyTemplates(:,pIdx);
            end
            memberInds = find(assmblPtrnWgtsTemp(:,pIdx) > ...
                (mean(assmblPtrnWgtsTemp(:,pIdx)) + 1*std(assmblPtrnWgtsTemp(:,pIdx))));
            if length(memberInds)>1
                validInd(pIdx) = 1;
            end
            AssmblPtrnCellIDs_L{pIdx} = cellID(memberInds);
        end
        AssmblPtrnCellIDs_L = AssmblPtrnCellIDs_L(validInd==1);
        AssmblWght_L = assmblPtrnWgtsTemp(:,validInd==1); 
        
        AssmblStrength_LL = assembly_activity(AssmblWght_L,seqTotal_L);
        AssmblStrength_LR = assembly_activity(AssmblWght_L,seqTotal_R);


        % Right delay
        AssemblyTemplates = assembly_patterns(seqTotal_R,opts);
        AssemblyTemplates2 = AssemblyTemplates;
        nPatterns = size(AssemblyTemplates,2);
        AssmblPtrnCellIDs_R=cell(1,nPatterns);
        validInd = zeros(1,nPatterns);
        assmblPtrnWgtsTemp = [];
        for pIdx = 1:nPatterns
            assmblPtrnWgtsTemp(:,pIdx) = AssemblyTemplates(:,pIdx)./norm(AssemblyTemplates(:,pIdx));
            if max(assmblPtrnWgtsTemp(:,pIdx)) < abs(min(assmblPtrnWgtsTemp(:,pIdx)))
                assmblPtrnWgtsTemp(:,pIdx) = -assmblPtrnWgtsTemp(:,pIdx);
                AssemblyTemplates2(:,pIdx) = -AssemblyTemplates(:,pIdx);
            end
            memberInds = find(assmblPtrnWgtsTemp(:,pIdx) > ...
                (mean(assmblPtrnWgtsTemp(:,pIdx)) + 1*std(assmblPtrnWgtsTemp(:,pIdx))));
            if length(memberInds)>1
                validInd(pIdx) = 1;
            end
            AssmblPtrnCellIDs_R{pIdx} = cellID(memberInds);
        end
        AssmblPtrnCellIDs_R = AssmblPtrnCellIDs_R(validInd==1);
        
        AssmblWght_R = assmblPtrnWgtsTemp(:,validInd==1);  
        AssmblStrength_RR = assembly_activity(AssmblWght_R,seqTotal_R);
        AssmblStrength_RL = assembly_activity(AssmblWght_R,seqTotal_L);

        % write down the results for the delay area
        delay.seqTotal_L = seqTotal_L;
        delay.seqTotal_R = seqTotal_R;
        delay.AssmblWght_L = AssmblWght_L;
        delay.AssmblWght_R = AssmblWght_R;
        delay.AssmblPtrnCellIDs_L = AssmblPtrnCellIDs_L;
        delay.AssmblPtrnCellIDs_R = AssmblPtrnCellIDs_R;
        delay.AssmblStrength_LL = AssmblStrength_LL;
        delay.AssmblStrength_LR = AssmblStrength_LR;
        delay.AssmblStrength_RR = AssmblStrength_RR;
        delay.AssmblStrength_RL = AssmblStrength_RL;
        CellAseembly.(sessDirs{j}).delay = delay;
        
        
        % calculate the strength of the original assembly in L/R shuffled
        % activity patterns
        % shuffle for p.shuffle times
        % get 95% of the shuffled strength and compare with the original
        % strength
        
        
        
        
        
        
        
        
        
        
        
        
        
        %% reward area
        rewardTime = [PathZone.posStartT.Reward,PathZone.posEndT.Reward]; % reward
        seqTotal_L = [];
        seqTotal_R = [];
        for m = 1:trialNum
            rewardBinIdx = rateBin1>=rewardTime(m,1) & rateBin2<=rewardTime(m,2);
            if strcmp(tInfo.direction{m},'L')
                seqTotal_L = [seqTotal_L,spkTrain(rateLabel,rewardBinIdx)];
            elseif strcmp(tInfo.direction{m},'R')
                seqTotal_R = [seqTotal_R,spkTrain(rateLabel,rewardBinIdx)];
            else
                error('Direction label error')
            end
        end
        
        % first : get asseembly on Left and Right
        % second: get assembly strength L-L,L-R,R-R,R-L
        % 
        % get assemblies
        % Left reward
        AssemblyTemplates = assembly_patterns(seqTotal_L,opts);
        AssemblyTemplates2 = AssemblyTemplates;
        % define if it is a cell assembly
        cellID = 1:sum(rateLabel);
        nPatterns = size(AssemblyTemplates,2);
        AssmblPtrnCellIDs_L=cell(1,nPatterns);
        validInd = zeros(1,nPatterns);
        assmblPtrnWgtsTemp = [];
        for pIdx = 1:nPatterns
            assmblPtrnWgtsTemp(:,pIdx) = AssemblyTemplates(:,pIdx)./norm(AssemblyTemplates(:,pIdx));
            if max(assmblPtrnWgtsTemp(:,pIdx)) < abs(min(assmblPtrnWgtsTemp(:,pIdx)))
                assmblPtrnWgtsTemp(:,pIdx) = -assmblPtrnWgtsTemp(:,pIdx);
                AssemblyTemplates2(:,pIdx) = -AssemblyTemplates(:,pIdx);
            end
            memberInds = find(assmblPtrnWgtsTemp(:,pIdx) > ...
                (mean(assmblPtrnWgtsTemp(:,pIdx)) + 1*std(assmblPtrnWgtsTemp(:,pIdx))));
            if length(memberInds)>1
                validInd(pIdx) = 1;
            end
            AssmblPtrnCellIDs_L{pIdx} = cellID(memberInds);
        end
        AssmblPtrnCellIDs_L = AssmblPtrnCellIDs_L(validInd==1);
        AssmblWght_L = assmblPtrnWgtsTemp(:,validInd==1);
        
        AssmblStrength_LL = assembly_activity(AssmblWght_L,seqTotal_L);
        AssmblStrength_LR = assembly_activity(AssmblWght_L,seqTotal_R);
        
 
        % Right reward
        AssemblyTemplates = assembly_patterns(seqTotal_R,opts);
        AssemblyTemplates2 = AssemblyTemplates;
        nPatterns = size(AssemblyTemplates,2);
        AssmblPtrnCellIDs_R=cell(1,nPatterns);
        validInd = zeros(1,nPatterns);
        assmblPtrnWgtsTemp = [];
        for pIdx = 1:nPatterns
            assmblPtrnWgtsTemp(:,pIdx) = AssemblyTemplates(:,pIdx)./norm(AssemblyTemplates(:,pIdx));
            if max(assmblPtrnWgtsTemp(:,pIdx)) < abs(min(assmblPtrnWgtsTemp(:,pIdx)))
                assmblPtrnWgtsTemp(:,pIdx) = -assmblPtrnWgtsTemp(:,pIdx);
                AssemblyTemplates2(:,pIdx) = -AssemblyTemplates(:,pIdx);
            end
            memberInds = find(assmblPtrnWgtsTemp(:,pIdx) > ...
                (mean(assmblPtrnWgtsTemp(:,pIdx)) + 1*std(assmblPtrnWgtsTemp(:,pIdx))));
            if length(memberInds)>1
                validInd(pIdx) = 1;
            end
            AssmblPtrnCellIDs_R{pIdx} = cellID(memberInds);
        end
        AssmblPtrnCellIDs_R = AssmblPtrnCellIDs_R(validInd==1);
        
        AssmblWght_R = assmblPtrnWgtsTemp(:,validInd==1); 
        AssmblStrength_RR = assembly_activity(AssmblWght_R,seqTotal_R);
        AssmblStrength_RL = assembly_activity(AssmblWght_R,seqTotal_L);
        

        % write down the results for the delay area
        reward.seqTotal_L = seqTotal_L;
        reward.seqTotal_R = seqTotal_R;
        reward.AssmblWght_L = AssmblWght_L;
        reward.AssmblWght_R = AssmblWght_R;
        reward.AssmblPtrnCellIDs_L = AssmblPtrnCellIDs_L;
        reward.AssmblPtrnCellIDs_R = AssmblPtrnCellIDs_R;
        reward.AssmblStrength_LL = AssmblStrength_LL;
        reward.AssmblStrength_LR = AssmblStrength_LR;
        reward.AssmblStrength_RR = AssmblStrength_RR;
        reward.AssmblStrength_RL = AssmblStrength_RL;
        CellAseembly.(sessDirs{j}).reward = reward;
        
    end  
    
    if p.writeToFile == 1
        save(fullfile(savedir2,'CellAseembly_DelayReward.mat'), 'CellAseembly');
    end
    clear CellAseembly
    fprintf('Finished analysis for session %d\n',i)
end
end
