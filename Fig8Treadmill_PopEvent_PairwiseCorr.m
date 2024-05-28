function Fig8Treadmill_PopEvent_PairwiseCorr(inFile,AnalyzeSes)

close all
p.savePlot = 0;
p.writeToFile = 1;

p.binWidth = 20./10^3; % unit sec

% % Read in input information
sessInfo = SessInfoImport(inFile);
dayCount = 0;

for i = AnalyzeSes(1:end)
    
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
    PopEvent_CellPairwiseCorr.rat = sessInfo(i).animal;
    PopEvent_CellPairwiseCorr.day = sessInfo(i).day;
    PopEvent_CellPairwiseCorr.binWidth = p.binWidth;
    
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    if sum(rateLabel) >= 20
        dayCount = dayCount + 1;
        for j = 1:length(sessDirs)
            
            % pop event combine
            % reward area
            % 1 return, 2 delay, 3 stem, 4 choice, 5 reward
            % quantify spike attend event ratio
            % quantify spike rate in all events combined
            rewardIdx = DelayPopFire_WholeSes.(sessDirs{j}).eventLoc==5;
            rewardInd = find(rewardIdx==1);
            
            eventNum = sum(rewardIdx);
            eventStartTs = DelayPopFire_WholeSes.(sessDirs{j}).eventStartTs(rewardIdx);
            eventEndTs = DelayPopFire_WholeSes.(sessDirs{j}).eventEndTs(rewardIdx);
            
            pairSize = nchoosek(clusterNum,2);
            rewardpop_CellPairCorr = nan(pairSize,eventNum);
            if eventNum>0
                for n = 1:eventNum
                    % calculate spk count in each time bin in the event
                    binTemp = 0:p.binWidth:eventEndTs(n)-eventStartTs(n);
                    spkCountTemp = zeros(clusterNum,length(binTemp)-1);                    
                    for k = 1:clusterNum
                        if rateLabel(k) == 1
                            spkTemp = DelayPopFire_WholeSes.(sessDirs{j}).eventTsp{rewardInd(n),k};
                            countTemp = histcounts(spkTemp,binTemp);
                            spkCountTemp(k,:) = countTemp;
                        end
                    end
                    % calculate pairwise correlation
                    pairCount = 0;
                    for kk = 1:clusterNum-1
                        for mm = kk+1:clusterNum
                            pairCount = pairCount + 1;
                            rewardpop_CellPairCorr(pairCount,n) = corr(spkCountTemp(kk,:)',spkCountTemp(mm,:)','Rows','pairwise');
                        end
                    end
                end
            end   

            PopEvent_CellPairwiseCorr.(sessDirs{j}).reward_EventNum = eventNum;
            PopEvent_CellPairwiseCorr.(sessDirs{j}).reward_DirectionLabel = DelayPopFire_WholeSes.(sessDirs{j}).eventDirectionLabel(rewardIdx);
            PopEvent_CellPairwiseCorr.(sessDirs{j}).reward_pop_CellPairCorr = rewardpop_CellPairCorr;
            
            
            % delay area
            eventNum = length(DelayPopFire_Delay.(sessDirs{j}).eventStartTs);
            eventStartTs = DelayPopFire_Delay.(sessDirs{j}).eventStartTs;
            eventEndTs = DelayPopFire_Delay.(sessDirs{j}).eventEndTs;
            
            pairSize = nchoosek(clusterNum,2);
            delaypop_CellPairCorr = nan(pairSize,eventNum);
            if eventNum>0
                for n = 1:eventNum
                    % calculate spk count in each time bin in the event
                    binTemp = 0:p.binWidth:eventEndTs(n)-eventStartTs(n);
                    spkCountTemp = zeros(clusterNum,length(binTemp)-1);                    
                    for k = 1:clusterNum
                        if rateLabel(k) == 1
                            spkTemp = DelayPopFire_Delay.(sessDirs{j}).eventTsp{n,k};
                            countTemp = histcounts(spkTemp,binTemp);
                            spkCountTemp(k,:) = countTemp;
                        end
                    end
                    % calculate pairwise correlation
                    pairCount = 0;
                    for kk = 1:clusterNum-1
                        for mm = kk+1:clusterNum
                            pairCount = pairCount + 1;
                            delaypop_CellPairCorr(pairCount,n) = corr(spkCountTemp(kk,:)',spkCountTemp(mm,:)','Rows','pairwise');
                        end
                    end
                end
            end   
            PopEvent_CellPairwiseCorr.(sessDirs{j}).delay_EventNum = eventNum;
            PopEvent_CellPairwiseCorr.(sessDirs{j}).delay_DirectionLabel = DelayPopFire_Delay.(sessDirs{j}).eventTurnDir;
            PopEvent_CellPairwiseCorr.(sessDirs{j}).delay_pop_CellPairCorr = delaypop_CellPairCorr;
        end       

        if p.writeToFile == 1
            save(fullfile(savedir2,'PopEvent_CellPairwiseCorr.mat'), 'PopEvent_CellPairwiseCorr');
        end
        clear PopEvent_CellPairwiseCorr
    end
    fprintf('Finished analysis for session %d\n',i)
end

end