function Fig8Treadmill_ActiveCell_PopEvent_Quant(inFile,AnalyzeSes)

close all
p.linROI = [32, 46];
p.linChoice = [55, 64];
p.xTickLin = [16 40 50 60 68];
p.xTickLabel = {'Return','Dl','St','Ch','Rw'};

p.savePlot = 0;
p.writeToFile = 0;

p.avgRateThres = 1;

avgRate = [];
allClusterNum = 0;
delayMeanRate_Def1 = [];
% initiate the map
map_All = [];
dayCount  = 0;

% pop event
rewardSpkEventRatio = [];
delaySpkEventRatio = [];

rewardSpkEventRate = [];
delaySpkEventRate = [];
        
% % Read in input information
sessInfo = SessInfoImport(inFile);

if p.savePlot
    % directory for plot figures
    % generate a folder for each rat eah day under the current folder
    savedir = sprintf('%s%s',cd,'\Figures\Quantification Figures\Time cell num-dist quant');
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end
end

for i = AnalyzeSes(1:end)
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load delayFire map
    delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap.mat');
    load(delayFile);
    % load population file
    reactFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_WholeSes.mat');
    load(reactFile);
    % load reactivation delay specific file
    reactDelayFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_Delay.mat');
    load(reactDelayFile);
    
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    sessDirs = sessInfo(i).sessDirs;
    for j = 1:length(sessDirs)
            % load time map
            if i == AnalyzeSes(1)
                timeMap_Def1.(sessDirs{j}) = [];
                popEvent.(sessDirs{j}).Reward = [];
                popEvent.(sessDirs{j}).Delay = [];
            end
    end
    
    if sum(rateLabel) >= 20
        dayCount = dayCount + 1;
        allClusterNum = allClusterNum + clusterNum;
        avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
    
        % get each phase names (no delay etc)

        delayMeanTemp_Def1 = zeros(clusterNum,length(sessDirs));
        
        spikeMap = cell(clusterNum,1);
        occupMap = cell(clusterNum,1);
        
        for j = 1:length(sessDirs)
            
            cellCount = 0;
            for k = 1:clusterNum
                cellCount = cellCount+1;
                
                timeMap_Def1.(sessDirs{j}) = [timeMap_Def1.(sessDirs{j});DelayFire.(sessDirs{j}).spikeRate1_Combined_Smooth{k}];
                delayMeanTemp_Def1(cellCount,j) = mean(DelayFire.(sessDirs{j}).spikeRate1_Combined{k});
            end
            
            % load spatial map
            % load analyzed map per trial
            load(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'));
            
            trialNum = length(ratesByECLR.ECLR);
            % SIMPLE WAY: put all trials together without considering
            % differences
            % Accurate way: match left to left and right to right (return,
            % base, choice, reward)
            for n = 1:trialNum
                for k = 1:clusterNum
                    spikeMap{k} = [spikeMap{k};ratesByECLR.spikemap{k,n}];
                    occupMap{k} = [occupMap{k};ratesByECLR.occupmap{k,n}];
                end
            end
            
            %% pop event combine
            % reward area
            % 1 return, 2 delay, 3 stem, 4 choice, 5 reward
            % quantify spike attend event ratio
            % quantify spike rate in all events combined
            rewardIdx = DelayPopFire_WholeSes.(sessDirs{j}).eventLoc==5;
            rewardInd = find(rewardIdx==1);
            
            eventNum = sum(rewardIdx);
            eventStartTs = DelayPopFire_WholeSes.(sessDirs{j}).eventStartTs(rewardIdx);
            eventEndTs = DelayPopFire_WholeSes.(sessDirs{j}).eventEndTs(rewardIdx);
            
            spkRatio = zeros(clusterNum,1);
            spkCount = zeros(clusterNum,1);
            spkRate = zeros(clusterNum,1);
            
            if eventNum>0
                for n = 1:length(rewardInd)
                    eventIndTemp = rewardInd(n);
                    for k = 1:clusterNum
                        spkCount(k) = spkCount(k) + length(DelayPopFire_WholeSes.(sessDirs{j}).eventTsp{eventIndTemp,k});
                        if length(DelayPopFire_WholeSes.(sessDirs{j}).eventTsp{eventIndTemp,k})>=1
                            spkRatio(k) = spkRatio(k) + 1;
                        end
                    end
                end
                spkRate = spkCount./(sum(eventEndTs-eventStartTs));
            end
         
            reward.eventNum = eventNum;
            reward.eventStartTs = eventStartTs;
            reward.eventEndTs = eventEndTs;
            reward.spkRatio = spkRatio;
            reward.spkCount = spkCount;
            reward.spkRate = spkRate;           
            popEvent.(sessDirs{j}).Reward = reward;
            
            % delay area
            % quantify spike attend event ratio
            % quantify spike rate in all events combined
            

            eventNum = length(DelayPopFire_Delay.(sessDirs{j}).eventStartTs);
            eventStartTs = DelayPopFire_Delay.(sessDirs{j}).eventStartTs;
            eventEndTs = DelayPopFire_Delay.(sessDirs{j}).eventEndTs;
            
            spkRatio = zeros(clusterNum,1);
            spkCount = zeros(clusterNum,1);
            spkRate = zeros(clusterNum,1);
            
            if eventNum>0
                for n = 1:eventNum
                    for k = 1:clusterNum
                        spkCount(k) = spkCount(k) + length(DelayPopFire_Delay.(sessDirs{j}).eventTsp{n,k});
                        if length(DelayPopFire_Delay.(sessDirs{j}).eventTsp{n,k})>=1
                            spkRatio(k) = spkRatio(k) + 1;
                        end
                    end
                end
                spkRate = spkCount./(sum(eventEndTs-eventStartTs));
            end
         
            delay.eventNum = eventNum;
            delay.eventStartTs = eventStartTs;
            delay.eventEndTs = eventEndTs;
            delay.spkRatio = spkRatio;
            delay.spkCount = spkCount;
            delay.spkRate = spkRate;           
            popEvent.(sessDirs{j}).Delay = delay;
            
        end
        
        rewardEventNumTemp = 0;
        delayEventNumTemp = 0;
        
        rewardEventLengthTemp = 0;
        delayEventLengthTemp = 0;
        
        rewardSpkEventTemp = zeros(clusterNum,1);
        delaySpkEventTemp = zeros(clusterNum,1);
        
        rewardSpkCountTemp = zeros(clusterNum,1);
        delaySpkCountTemp = zeros(clusterNum,1);
        
        for j = 1:length(sessDirs)
            rewardEventNumTemp = rewardEventNumTemp + popEvent.(sessDirs{j}).Reward.eventNum;
            delayEventNumTemp = delayEventNumTemp + popEvent.(sessDirs{j}).Delay.eventNum;
            
            rewardEventLengthTemp = rewardEventLengthTemp + sum(popEvent.(sessDirs{j}).Reward.eventEndTs-popEvent.(sessDirs{j}).Reward.eventStartTs);
            delayEventLengthTemp = delayEventLengthTemp + sum(popEvent.(sessDirs{j}).Delay.eventEndTs-popEvent.(sessDirs{j}).Delay.eventStartTs);
            
            rewardSpkEventTemp = rewardSpkEventTemp + popEvent.(sessDirs{j}).Reward.spkRatio;
            delaySpkEventTemp = delaySpkEventTemp + popEvent.(sessDirs{j}).Delay.spkRatio;
            
            rewardSpkCountTemp = rewardSpkCountTemp + popEvent.(sessDirs{j}).Reward.spkCount;
            delaySpkCountTemp = delaySpkCountTemp + popEvent.(sessDirs{j}).Delay.spkCount;
        end
        
        rewardSpkEventRatio = [rewardSpkEventRatio;rewardSpkEventTemp./rewardEventNumTemp];
        delaySpkEventRatio = [delaySpkEventRatio;delaySpkEventTemp./delayEventNumTemp];
        
        rewardSpkEventRate = [rewardSpkEventRate;rewardSpkCountTemp./rewardEventLengthTemp];
        delaySpkEventRate = [delaySpkEventRate;delaySpkCountTemp./delayEventLengthTemp];
        
        
        
        %% delay avg rate        
        delayMeanRate_Def1 = [delayMeanRate_Def1;delayMeanTemp_Def1];
        
        % all maps
        for k = 1:clusterNum
            spikeMap_AllTrial = sum(spikeMap{k},1);
            occupMap_AllTrial = sum(occupMap{k},1);
            map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
            map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
            map_All = [map_All;map_Temp_2];
        end
    end
end

pyrIdx = avgRate>0.1 & avgRate<5;

%% delay definition 1: delay aligned by barrier
% as long as it is over average rate threshold in one session
avgRate_pyr = avgRate(pyrIdx);
map_Pyr = map_All(pyrIdx,:);
delayMeanRate_Def1 = delayMeanRate_Def1(pyrIdx,:);

rewardSpkEventRatio = rewardSpkEventRatio(pyrIdx);
delaySpkEventRatio = delaySpkEventRatio(pyrIdx);

rewardSpkEventRate = rewardSpkEventRate(pyrIdx);
delaySpkEventRate = delaySpkEventRate(pyrIdx);

% get delay active cells
delay_onIdx = (sum(delayMeanRate_Def1>p.avgRateThres,2))>0;

% get maze active cells
mazePeakRate = max(map_Pyr(:,[1:p.linROI(1)-2,p.linROI(2)+2:end]),[],2);
mazeMeanRate = mean(map_Pyr(:,[1:p.linROI(1)-2,p.linROI(2)+2:end]),2);
mazeMeanRate_SD = std(map_Pyr(:,[1:p.linROI(1)-2,p.linROI(2)+2:end]),0,2);

maze_onIdx_1 = mazePeakRate > 3;
maze_onIdx_2 = mazePeakRate > (mazeMeanRate+2*mazeMeanRate_SD);
maze_onIdx = maze_onIdx_1 & maze_onIdx_2;

% plot all cells
figure(1)
[~,peakIdx] = max(map_Pyr,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = map_Pyr(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);

imagesc(cellMapTemp1_Sort_Norm)
colormap(jet)
set(gca, 'xtick', p.xTickLin);
set(gca, 'xticklabels', p.xTickLabel);
hold on; vertplot(p.linROI, 0, size(cellMapTemp1_Sort,1), 'r--');
vertplot(p.linChoice, 0, size(cellMapTemp1_Sort,1), 'r--');


figure(2)
maze_On_map = map_Pyr(maze_onIdx,:);
[~,peakIdx] = max(maze_On_map,[],2);
% consider cells outside reward area
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = maze_On_map(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);

imagesc(cellMapTemp1_Sort)
colormap(jet)
set(gca, 'xtick', p.xTickLin);
set(gca, 'xticklabels', p.xTickLabel);
hold on; vertplot(p.linROI, 0, size(cellMapTemp1_Sort,1), 'r--');
vertplot(p.linChoice, 0, size(cellMapTemp1_Sort,1), 'r--');
title('All maze on cells');
caxis([0 max(max(cellMapTemp1_Sort))/4])

figure(3)
delay_On_map = map_Pyr(delay_onIdx,:);
[~,peakIdx] = max(delay_On_map,[],2);
% consider cells outside reward area
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = delay_On_map(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);

imagesc(cellMapTemp1_Sort)
colormap(jet)
set(gca, 'xtick', p.xTickLin);
set(gca, 'xticklabels', p.xTickLabel);
hold on; vertplot(p.linROI, 0, size(cellMapTemp1_Sort,1), 'r--');
vertplot(p.linChoice, 0, size(cellMapTemp1_Sort,1), 'r--');
title('All Delay on cells');
caxis([0 max(max(cellMapTemp1_Sort))/4])

mazeCell_RewardPopRatio = rewardSpkEventRatio(maze_onIdx);
mazeCell_DelayPopRatio = delaySpkEventRatio(maze_onIdx);
mazeCell_RewardPopRate= rewardSpkEventRate(maze_onIdx);
mazeCell_DelayPopRate = delaySpkEventRate(maze_onIdx);

delayCell_RewardPopRatio = rewardSpkEventRatio(delay_onIdx);
delayCell_DelayPopRatio = delaySpkEventRatio(delay_onIdx);
delayCell_RewardPopRate= rewardSpkEventRate(delay_onIdx);
delayCell_DelayPopRate = delaySpkEventRate(delay_onIdx);

figure(2)
[~,peakIdx] = max(maze_On_map,[],2);
[~,peak_Sort] = sort(peakIdx);
ratioTemp1 = movmean(mazeCell_RewardPopRatio(peak_Sort),10)*20;
plot(74+ratioTemp1,1:length(mazeCell_RewardPopRatio),'k')
ratioTemp2 = movmean(mazeCell_DelayPopRatio(peak_Sort),10)*20;
plot(74+ratioTemp2,1:length(mazeCell_RewardPopRatio),'r')
xlim([0 90])

figure(3)
[~,peakIdx] = max(delay_On_map,[],2);
[~,peak_Sort] = sort(peakIdx);
ratioTemp1 = movmean(delayCell_RewardPopRatio(peak_Sort),10)*20;
plot(74+ratioTemp1,1:length(delayCell_RewardPopRatio),'k')
ratioTemp2 = movmean(delayCell_DelayPopRatio(peak_Sort),10)*20;
plot(74+ratioTemp2,1:length(delayCell_DelayPopRatio),'r')
xlim([0 90])


figure(4)
Violin(delayCell_RewardPopRatio,1,'ViolinColor',[0.2,0.2,0.2],'ShowData',false);
Violin(mazeCell_RewardPopRatio,2,'ViolinColor',[0.2,0.2,0.2],'ShowData',false);

Violin(delayCell_DelayPopRatio,4,'ViolinColor',[1,0.2,0.2],'ShowData',false);
Violin(mazeCell_DelayPopRatio,5,'ViolinColor',[1,0.2,0.2],'ShowData',false);
set(gca, 'XTick', [1,2,4,5], 'XTickLabel', {'Delay on','Maze on','Delay on','Maze on'});
title('Pop attending ratio')

figure(5)
Violin(delayCell_RewardPopRate,1,'ViolinColor',[0.2,0.2,0.2],'ShowData',false);
Violin(mazeCell_RewardPopRate,2,'ViolinColor',[0.2,0.2,0.2],'ShowData',false);

Violin(delayCell_DelayPopRate,4,'ViolinColor',[1,0.2,0.2],'ShowData',false);
Violin(mazeCell_DelayPopRate,5,'ViolinColor',[1,0.2,0.2],'ShowData',false);
set(gca, 'XTick', [1,2,4,5], 'XTickLabel', {'Delay on','Maze on','Delay on','Maze on'});
title('Pop attending rate')

%% remove cells which have peak in the reward area
figure(6)
maze_On_map = map_Pyr(maze_onIdx,:);
[~,peakIdx] = max(maze_On_map,[],2);
non_RewardIdx = peakIdx < 64;
non_Reward_peakIdx = peakIdx(non_RewardIdx);
maze_nonReward_map = maze_On_map(non_RewardIdx,:);
% consider cells outside reward area
[~,peak_Sort] = sort(non_Reward_peakIdx);

cellMapTemp1_Sort = maze_nonReward_map(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);

imagesc(cellMapTemp1_Sort)
colormap(jet)
set(gca, 'xtick', p.xTickLin);
set(gca, 'xticklabels', p.xTickLabel);
hold on; vertplot(p.linROI, 0, size(cellMapTemp1_Sort,1), 'r--');
vertplot(p.linChoice, 0, size(cellMapTemp1_Sort,1), 'r--');
title('All maze on cells');
caxis([0 max(max(cellMapTemp1_Sort))/4])

mazeCell_nonReward_RewardPopRatio = mazeCell_RewardPopRatio(non_RewardIdx);
mazeCell_nonReward__DelayPopRatio = mazeCell_DelayPopRatio(non_RewardIdx);
mazeCell_nonReward__RewardPopRate= mazeCell_RewardPopRate(non_RewardIdx);
mazeCell_nonReward__DelayPopRate = mazeCell_DelayPopRate(non_RewardIdx);

figure(6)
ratioTemp1 = movmean(mazeCell_nonReward_RewardPopRatio(peak_Sort),10)*20;
plot(74+ratioTemp1,1:length(mazeCell_nonReward_RewardPopRatio),'k')
ratioTemp2 = movmean(mazeCell_nonReward__DelayPopRatio(peak_Sort),10)*20;
plot(74+ratioTemp2,1:length(mazeCell_nonReward__DelayPopRatio),'r')
xlim([0 90])

figure(7)
Violin(mazeCell_nonReward_RewardPopRatio,1,'ViolinColor',[0.2,0.2,0.2]);
Violin(delayCell_RewardPopRatio,2,'ViolinColor',[0.2,0.2,0.2]);

Violin(mazeCell_nonReward__DelayPopRatio,4,'ViolinColor',[1,0.2,0.2]);
Violin(delayCell_DelayPopRatio,5,'ViolinColor',[1,0.2,0.2]);

figure(8)
Violin(mazeCell_nonReward__RewardPopRate,1,'ViolinColor',[0.2,0.2,0.2]);
Violin(delayCell_RewardPopRate,2,'ViolinColor',[0.2,0.2,0.2]);

Violin(mazeCell_nonReward__DelayPopRate,4,'ViolinColor',[1,0.2,0.2]);
Violin(delayCell_DelayPopRate,5,'ViolinColor',[1,0.2,0.2]);


% get exclusive cells
delay_On_Ex_Idx = delay_onIdx & ~maze_onIdx;
maze_On_Ex_Idx = ~delay_onIdx & maze_onIdx;

figure(2)
maze_On_Ex_map = map_Pyr(maze_On_Ex_Idx,:);
[~,peakIdx] = max(maze_On_Ex_map,[],2);
maze_On_nonReward_Idx = peakIdx <= 64;

% consider cells outside reward area
maze_On_nonReward_map = maze_On_Ex_map(maze_On_nonReward_Idx,:);
[~,peak_Sort] = sort(peakIdx(maze_On_nonReward_Idx));
cellMapTemp1_Sort = maze_On_nonReward_map(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);

imagesc(cellMapTemp1_Sort_Norm)
colormap(jet)
set(gca, 'xtick', p.xTickLin);
set(gca, 'xticklabels', p.xTickLabel);
hold on; vertplot(p.linROI, 0, size(cellMapTemp1_Sort,1), 'r--');
vertplot(p.linChoice, 0, size(cellMapTemp1_Sort,1), 'r--');

figure(3)
delay_On_Ex_map = map_Pyr(delay_On_Ex_Idx,:);
[~,peakIdx] = max(delay_On_Ex_map,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = delay_On_Ex_map(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);

imagesc(cellMapTemp1_Sort_Norm)
colormap(jet)
set(gca, 'xtick', p.xTickLin);
set(gca, 'xticklabels', p.xTickLabel);
hold on; vertplot(p.linROI, 0, size(cellMapTemp1_Sort,1), 'r--');
vertplot(p.linChoice, 0, size(cellMapTemp1_Sort,1), 'r--');


%% Pop event combine
% delay_On_Ex_Idx = delay_onIdx & ~maze_onIdx;
% maze_On_Ex_Idx = ~delay_onIdx & maze_onIdx;
mazeCell_RewardPopRatio = rewardSpkEventRatio(maze_On_Ex_Idx);
mazeCell_DelayPopRatio = delaySpkEventRatio(maze_On_Ex_Idx);
mazeCell_RewardPopRate= rewardSpkEventRate(maze_On_Ex_Idx);
mazeCell_DelayPopRate = delaySpkEventRate(maze_On_Ex_Idx);

mazeCell_RewardPopRatio = mazeCell_RewardPopRatio(maze_On_nonReward_Idx);
mazeCell_DelayPopRatio = mazeCell_DelayPopRatio(maze_On_nonReward_Idx);
mazeCell_RewardPopRate= mazeCell_RewardPopRate(maze_On_nonReward_Idx);
mazeCell_DelayPopRate = mazeCell_DelayPopRate(maze_On_nonReward_Idx);


delayCell_RewardPopRatio = rewardSpkEventRatio(delay_On_Ex_Idx);
delayCell_DelayPopRatio = delaySpkEventRatio(delay_On_Ex_Idx);
delayCell_RewardPopRate= rewardSpkEventRate(delay_On_Ex_Idx);
delayCell_DelayPopRate = delaySpkEventRate(delay_On_Ex_Idx);

figure(4)
Violin(mazeCell_RewardPopRatio,1)
Violin(mazeCell_DelayPopRatio,2)
Violin(delayCell_RewardPopRatio,4)
Violin(delayCell_DelayPopRatio,5)

figure(5)
Violin(mazeCell_RewardPopRate,1)
Violin(mazeCell_DelayPopRate,2)
Violin(delayCell_RewardPopRate,4)
Violin(delayCell_DelayPopRate,5)

% figure(6)
% Violin(avgRate_pyr(maze_On_Ex_Idx((maze_On_nonReward_Idx))),1)
% Violin(avgRate_pyr(delay_On_Ex_Idx),2)


figure(2)
[~,peakIdx] = max(maze_On_nonReward_map,[],2);
[~,peak_Sort] = sort(peakIdx);
plot(74+mazeCell_RewardPopRatio(peak_Sort)*30,1:length(mazeCell_RewardPopRatio),'k')
plot(74+mazeCell_DelayPopRatio(peak_Sort)*30,1:length(mazeCell_RewardPopRatio),'r')
xlim([0 100])

figure(3)
[~,peakIdx] = max(delay_On_Ex_map,[],2);
[~,peak_Sort] = sort(peakIdx);
plot(74+delayCell_RewardPopRatio(peak_Sort)*30,1:length(delayCell_RewardPopRatio),'k')
plot(74+delayCell_DelayPopRatio(peak_Sort)*30,1:length(delayCell_DelayPopRatio),'r')
xlim([0 100])

if p.writeToFile == 1
    save(fullfile(savedir,'timecellQuant_Timebin.mat'), 'timecellQuant_Timebin');
end
clear timecellQuant

end