function Fig8TreadmillArmChoice_Burst(inFile,AnalyzeSes,shuffleTimes)

%  Li Yuan, UCSD, 22-Mar-2021
%  Calculate arm selectivity based on average rate
% -------------------------------------------------------------------------
% set parameters for analysis

p.savePlot = 1;
p.writeToFile = 1;

% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    close all
    % selectivity for the whole map
    Fig8TreadmillArmChoice_Burst.session = sessInfo(i).mainDir;
    Fig8TreadmillArmChoice_Burst.rat = sessInfo(i).animal;
    Fig8TreadmillArmChoice_Burst.day = sessInfo(i).day;
    
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\ArmChoice_Burst');
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
    
    display(['About to do session ' sessInfo(i).mainDir]);
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    % load spikes
    % load burst and single spike
    burstFile = fullfile(sessInfo(i).mainDir,'Cell Property','BurstActivity.mat');
    load(burstFile); 
    TList = BurstActivity.TList;
    clusterNum = length(TList);
    if length(unique(TList))~=clusterNum
        error('TTList file has repeated clusters')
    end
    
    if prod(~isempty(char(sessDirs{1})))
        %% burst 
        map_1D = [];
        for sesGroup = 1:4
            sesNumSum = 0;
            posLabel = [];
            turnLabel = [];
            time_All = [];

            behaveType = strsplit(sessDirs{sesGroup},'_');
            behaveType = behaveType{1};
            
            for j = [sesGroup,sesGroup+4]
                
                % load maps for each cluster
                pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'PathZone.mat');
                load(pathZoneFile);
                delayFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
                load(delayFile);
                
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
                        
                if exist(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'), 'file')==2
                    map_1D{j} = load(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'));
                    
                    for m = 1:trialNum
                        if map_1D{j}.ratesByECLR.ECLR(m) == 1
                            posLabel(sesNumSum+m) = 2;
                            turnLabel(sesNumSum+m) = 1;
                        elseif map_1D{j}.ratesByECLR.ECLR(m) == 2
                            posLabel(sesNumSum+m) = 1;
                            turnLabel(sesNumSum+m) = 1;
                        elseif map_1D{j}.ratesByECLR.ECLR(m) == 3
                            posLabel(sesNumSum+m) = 1;
                            turnLabel(sesNumSum+m) = 0;
                        elseif map_1D{j}.ratesByECLR.ECLR(m) == 4
                            posLabel(sesNumSum+m) = 2;
                            turnLabel(sesNumSum+m) = 0;
                        else
                            error('Turning label ERROR');
                        end
                                                
                        
                        % return
                        timeTemp = PathZone.posEndT.Return(m)-PathZone.posStartT.Return(m);
                        time_All(sesNumSum+m,1) = timeTemp;
                        % delay
                        timeTemp = delayTend1_2(m)-delayTstart1(m);
                        time_All(sesNumSum+m,2) = timeTemp;
                        % stem
                        timeTemp = PathZone.posEndT.Center(m)-(Fig8DelayZonePos.delayPos1.endT(m)+1/30);
                        time_All(sesNumSum+m,3) = timeTemp;                        
                    end
                    sesNumSum = sesNumSum + sum(map_1D{j}.ratesByECLR.valid);
                end               
                
            end
            
            
            % find Left / Right turn trials
            leftIdx = find(turnLabel==1); % idx on all trials
            rightIdx = find(turnLabel==0); % idx on all trials
            
            
            for k = 1:clusterNum                                   
                sessionTemp = 1;
                % initiate variables
                % return, delay, stem
                spike_All = nan(sesNumSum,3);
                for j = [sesGroup,sesGroup+4]
                    % load path
                    pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'PathZone.mat');
                    load(pathZoneFile);
                    delayFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
                    load(delayFile);
                
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
                
                    % get spikes
                    %% burst
                    if ~isempty(BurstActivity.(sessDirs{j}).BurstProp{k})
                        tSp = BurstActivity.(sessDirs{j}).BurstProp{k}.BurststartTs;
                    else
                        tSp = [];
                    end
                    
                    for m = 1:sum(map_1D{j}.ratesByECLR.valid)
                        % get spike count
                        % return
                        spikeCount = sum(tSp>=PathZone.posStartT.Return(m) & tSp<PathZone.posEndT.Return(m));
                        spike_All(sessionTemp,1) = spikeCount;
                        % delay
                        spikeCount = sum(tSp>=delayTstart1(m) & tSp<delayTend1_2(m));
                        spike_All(sessionTemp,2) = spikeCount;
                        % stem
                        spikeCount = sum(tSp>=(Fig8DelayZonePos.delayPos1.endT(m)+1/30) & tSp<PathZone.posEndT.Center(m));
                        spike_All(sessionTemp,3) = spikeCount;
                        sessionTemp = sessionTemp+1;
                    end
                end
                
%                 % LR comparison
%                 rateL = nansum(nansum(spike_All(leftIdx,1)))./sum(sum(time_All(leftIdx,1)));
%                 rateR= nansum(nansum(spike_All(rightIdx,1)))./sum(sum(time_All(rightIdx,1)));
%                 LRrateDiff = (rateL-rateR)/(rateL+rateR);
%                 
                % delay zone (LR turning comparison)
                delayRateLAll = nansum(nansum(spike_All(leftIdx,2)))./sum(sum(time_All(leftIdx,2)));
                delayRateRAll = nansum(nansum(spike_All(rightIdx,2)))./sum(sum(time_All(rightIdx,2)));
                delayRateLRDiffAll = (delayRateLAll-delayRateRAll)/(delayRateLAll+delayRateRAll);
                
                % stem zone (LR turning comparison)
                stemRateLAll = nansum(nansum(spike_All(leftIdx,3)))./sum(sum(time_All(leftIdx,3)));
                stemRateRAll = nansum(nansum(spike_All(rightIdx,3)))./sum(sum(time_All(rightIdx,3)));
                stemRateLRDiffAll = (stemRateLAll-stemRateRAll)/(stemRateLAll+stemRateRAll);
     

%                 Fig8TreadmillCenterRateDiff.(behaveType).LRrateDiff(k) = LRrateDiff;
                Fig8TreadmillArmChoice_Burst.(behaveType).Burst.delayRateLRDiffAll(k) = delayRateLRDiffAll;
                Fig8TreadmillArmChoice_Burst.(behaveType).Burst.stemRateLRDiffAll(k) = stemRateLRDiffAll;
                
                
                % initiate shuffle quantification               
%                 shf_LRrateDiffAll = nan(shuffleTimes,1);
                shf_delayRateLRDiffAll = nan(shuffleTimes,1);
                shf_stemRateLRDiffAll = nan(shuffleTimes,1);
                
                
                for n = 1:shuffleTimes
                    
                    % all the followings are switch left & right
                    % all trial shulffle
                    % randi 0 = left, 1 = right
                    shfOrder = logical(randi([0,1],[sesNumSum,1]));
%                     shf_LrateAll = nansum(nansum(spike_All(~shfOrder,1)))./sum(sum(time_All(~shfOrder,1)));
%                     shf_RrateAll = nansum(nansum(spike_All(shfOrder,1)))./sum(sum(time_All(shfOrder,1)));
%                     shf_LRrateDiffAll(n) = (shf_LrateAll-shf_RrateAll)/(shf_LrateAll+shf_RrateAll);
%                     
                    shf_delayRateLall = nansum(nansum(spike_All(~shfOrder,2)))./nansum(nansum(time_All(~shfOrder,2)));
                    shf_delayRateRall = nansum(nansum(spike_All(shfOrder,2)))./nansum(nansum(time_All(shfOrder,2)));
                    shf_delayRateLRDiffAll(n) = (shf_delayRateLall-shf_delayRateRall)/(shf_delayRateLall+shf_delayRateRall);
                    
                    shf_stemRateLall = nansum(nansum(spike_All(~shfOrder,3)))./nansum(nansum(time_All(~shfOrder,3)));
                    shf_stemRateRall = nansum(nansum(spike_All(shfOrder,3)))./nansum(nansum(time_All(shfOrder,3)));
                    shf_stemRateLRDiffAll(n) = (shf_stemRateLall-shf_stemRateRall)/(shf_stemRateLall+shf_stemRateRall);
                    
                end 
                
                h = figure(k);
                if sesGroup == 1
                    h.Position = [100 100 1200 2000];
                end
                
                % D(L) vs D(R)
                subplot(2,4,sesGroup)
                tempSort = sort(shf_delayRateLRDiffAll(~isnan(shf_delayRateLRDiffAll)));
                if length(tempSort) > shuffleTimes/10
                    if delayRateLRDiffAll >= 0
                        shf_delayRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.975));
                        shf_delayRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.995));
                    else
                        shf_delayRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.025));
                        shf_delayRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.005));
                    end
                    Violin(shf_delayRateLRDiffAll,1);
                    hold on
                    plot([0.7,1.3],ones(1,2).*shf_delayRateLRDiffAll95,'r-')
                    plot(1,delayRateLRDiffAll,'r*')
                    hold off
                    Fig8TreadmillArmChoice_Burst.(behaveType).Burst.shf_delayRateLRDiffAll(k) = 1;
                else
                    shf_delayRateLRDiffAll95 = NaN;
                    shf_delayRateLRDiffAll99 = NaN;
                    text(0.1,0.5,'Comparison not exist','FontSize',8);
                    Fig8TreadmillArmChoice_Burst.(behaveType).Burst.shf_delayRateLRDiffAll(k) = 0;
                end
                Fig8TreadmillArmChoice_Burst.(behaveType).Burst.shf_delayRateLRDiffAll95(k) = shf_delayRateLRDiffAll95;
                Fig8TreadmillArmChoice_Burst.(behaveType).Burst.shf_delayRateLRDiffAll99(k) = shf_delayRateLRDiffAll99;
                TITLE1 = sprintf('%s%d%s%d%s%s','Arm Choice-Burst-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k}(1:end-2));
                TITLE2 = behaveType;
                title({TITLE1;TITLE2},'Interpreter','None')
                ylabel('Delay (L vs R)')
                xticklabels({''});
                yticklabels({''});
                
                % S(L) vs S(R)
                subplot(2,4,4+sesGroup)
                tempSort = sort(shf_stemRateLRDiffAll(~isnan(shf_stemRateLRDiffAll)));
                if length(tempSort) > shuffleTimes/10
                    if stemRateLRDiffAll >= 0
                        shf_stemRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.975));
                        shf_stemRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.995));
                    else
                        shf_stemRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.025));
                        shf_stemRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.005));
                    end
                    Violin(shf_stemRateLRDiffAll,1);
                    hold on
                    plot([0.7,1.3],ones(1,2).*shf_stemRateLRDiffAll95,'r-')
                    plot(1,stemRateLRDiffAll,'r*')
                    hold off
                    Fig8TreadmillArmChoice_Burst.(behaveType).Burst.shf_stemRateLRDiffAll(k) = 1;
                else
                    shf_stemRateLRDiffAll95 = NaN;
                    shf_stemRateLRDiffAll99 = NaN;
                    text(0.1,0.5,'Comparison not exist','FontSize',8);
                    Fig8TreadmillArmChoice_Burst.(behaveType).Burst.shf_stemRateLRDiffAll(k) = 0;
                end
                Fig8TreadmillArmChoice_Burst.(behaveType).Burst.shf_stemRateLRDiffAll95(k) = shf_stemRateLRDiffAll95;
                Fig8TreadmillArmChoice_Burst.(behaveType).Burst.shf_stemRateLRDiffAll99(k) = shf_stemRateLRDiffAll99;
                ylabel('Stem (L vs R)')
                xticklabels({''});
                yticklabels({''});
 
            end       
        end
    end
    
    % save figure and save .mat file
    if p.savePlot
        for k = 1:clusterNum
            figure(k)
            figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',TList{k}(1:end-2),'-ArmChoice_Burst');
            print(figName,'-dpng','-r300');
        end
    end
    close all
    
     if prod(~isempty(char(sessDirs{1})))
        %% single spike
        map_1D = [];
        for sesGroup = 1:4
            sesNumSum = 0;
            posLabel = [];
            turnLabel = [];
            time_All = [];

            behaveType = strsplit(sessDirs{sesGroup},'_');
            behaveType = behaveType{1};
            
            for j = [sesGroup,sesGroup+4]
                
                % load maps for each cluster
                pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'PathZone.mat');
                load(pathZoneFile);
                delayFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
                load(delayFile);
                
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
                        
                if exist(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'), 'file')==2
                    map_1D{j} = load(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'));
                    
                    for m = 1:trialNum
                        if map_1D{j}.ratesByECLR.ECLR(m) == 1
                            posLabel(sesNumSum+m) = 2;
                            turnLabel(sesNumSum+m) = 1;
                        elseif map_1D{j}.ratesByECLR.ECLR(m) == 2
                            posLabel(sesNumSum+m) = 1;
                            turnLabel(sesNumSum+m) = 1;
                        elseif map_1D{j}.ratesByECLR.ECLR(m) == 3
                            posLabel(sesNumSum+m) = 1;
                            turnLabel(sesNumSum+m) = 0;
                        elseif map_1D{j}.ratesByECLR.ECLR(m) == 4
                            posLabel(sesNumSum+m) = 2;
                            turnLabel(sesNumSum+m) = 0;
                        else
                            error('Turning label ERROR');
                        end
                                                
                        
                        % return
                        timeTemp = PathZone.posEndT.Return(m)-PathZone.posStartT.Return(m);
                        time_All(sesNumSum+m,1) = timeTemp;
                        % delay
                        timeTemp = delayTend1_2(m)-delayTstart1(m);
                        time_All(sesNumSum+m,2) = timeTemp;
                        % stem
                        timeTemp = PathZone.posEndT.Center(m)-(Fig8DelayZonePos.delayPos1.endT(m)+1/30);
                        time_All(sesNumSum+m,3) = timeTemp;                        
                    end
                    sesNumSum = sesNumSum + sum(map_1D{j}.ratesByECLR.valid);
                end               
                
            end
            
            
            % find Left / Right turn trials
            leftIdx = find(turnLabel==1); % idx on all trials
            rightIdx = find(turnLabel==0); % idx on all trials
            
            
            for k = 1:clusterNum                                   
                sessionTemp = 1;
                % initiate variables
                % return, delay, stem
                spike_All = nan(sesNumSum,3);
                for j = [sesGroup,sesGroup+4]
                    % load path
                    pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'PathZone.mat');
                    load(pathZoneFile);
                    delayFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
                    load(delayFile);
                
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
                
                    % get spikes
                    %% burst
                    if ~isempty(BurstActivity.(sessDirs{j}).SinglespikeProp{k})
                        tSp = BurstActivity.(sessDirs{j}).SinglespikeProp{k}.SinglespikeTs;
                    else
                        tSp = [];
                    end
                    
                    for m = 1:sum(map_1D{j}.ratesByECLR.valid)
                        % get spike count
                        % return
                        spikeCount = sum(tSp>=PathZone.posStartT.Return(m) & tSp<PathZone.posEndT.Return(m));
                        spike_All(sessionTemp,1) = spikeCount;
                        % delay
                        spikeCount = sum(tSp>=delayTstart1(m) & tSp<delayTend1_2(m));
                        spike_All(sessionTemp,2) = spikeCount;
                        % stem
                        spikeCount = sum(tSp>=(Fig8DelayZonePos.delayPos1.endT(m)+1/30) & tSp<PathZone.posEndT.Center(m));
                        spike_All(sessionTemp,3) = spikeCount;
                        sessionTemp = sessionTemp+1;
                    end
                end
                
%                 % LR comparison
%                 rateL = nansum(nansum(spike_All(leftIdx,1)))./sum(sum(time_All(leftIdx,1)));
%                 rateR= nansum(nansum(spike_All(rightIdx,1)))./sum(sum(time_All(rightIdx,1)));
%                 LRrateDiff = (rateL-rateR)/(rateL+rateR);
%                 
                % delay zone (LR turning comparison)
                delayRateLAll = nansum(nansum(spike_All(leftIdx,2)))./sum(sum(time_All(leftIdx,2)));
                delayRateRAll = nansum(nansum(spike_All(rightIdx,2)))./sum(sum(time_All(rightIdx,2)));
                delayRateLRDiffAll = (delayRateLAll-delayRateRAll)/(delayRateLAll+delayRateRAll);
                
                % stem zone (LR turning comparison)
                stemRateLAll = nansum(nansum(spike_All(leftIdx,3)))./sum(sum(time_All(leftIdx,3)));
                stemRateRAll = nansum(nansum(spike_All(rightIdx,3)))./sum(sum(time_All(rightIdx,3)));
                stemRateLRDiffAll = (stemRateLAll-stemRateRAll)/(stemRateLAll+stemRateRAll);
     

%                 Fig8TreadmillCenterRateDiff.(behaveType).LRrateDiff(k) = LRrateDiff;
                Fig8TreadmillArmChoice_Burst.(behaveType).SingleSpk.delayRateLRDiffAll(k) = delayRateLRDiffAll;
                Fig8TreadmillArmChoice_Burst.(behaveType).SingleSpk.stemRateLRDiffAll(k) = stemRateLRDiffAll;
                
                
                % initiate shuffle quantification               
%                 shf_LRrateDiffAll = nan(shuffleTimes,1);
                shf_delayRateLRDiffAll = nan(shuffleTimes,1);
                shf_stemRateLRDiffAll = nan(shuffleTimes,1);
                
                
                for n = 1:shuffleTimes
                    
                    % all the followings are switch left & right
                    % all trial shulffle
                    % randi 0 = left, 1 = right
                    shfOrder = logical(randi([0,1],[sesNumSum,1]));
%                     shf_LrateAll = nansum(nansum(spike_All(~shfOrder,1)))./sum(sum(time_All(~shfOrder,1)));
%                     shf_RrateAll = nansum(nansum(spike_All(shfOrder,1)))./sum(sum(time_All(shfOrder,1)));
%                     shf_LRrateDiffAll(n) = (shf_LrateAll-shf_RrateAll)/(shf_LrateAll+shf_RrateAll);
%                     
                    shf_delayRateLall = nansum(nansum(spike_All(~shfOrder,2)))./nansum(nansum(time_All(~shfOrder,2)));
                    shf_delayRateRall = nansum(nansum(spike_All(shfOrder,2)))./nansum(nansum(time_All(shfOrder,2)));
                    shf_delayRateLRDiffAll(n) = (shf_delayRateLall-shf_delayRateRall)/(shf_delayRateLall+shf_delayRateRall);
                    
                    shf_stemRateLall = nansum(nansum(spike_All(~shfOrder,3)))./nansum(nansum(time_All(~shfOrder,3)));
                    shf_stemRateRall = nansum(nansum(spike_All(shfOrder,3)))./nansum(nansum(time_All(shfOrder,3)));
                    shf_stemRateLRDiffAll(n) = (shf_stemRateLall-shf_stemRateRall)/(shf_stemRateLall+shf_stemRateRall);
                    
                end 
                
                h = figure(k);
                if sesGroup == 1
                    h.Position = [100 100 1200 2000];
                end
                
                % D(L) vs D(R)
                subplot(2,4,sesGroup)
                tempSort = sort(shf_delayRateLRDiffAll(~isnan(shf_delayRateLRDiffAll)));
                if length(tempSort) > shuffleTimes/10
                    if delayRateLRDiffAll >= 0
                        shf_delayRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.975));
                        shf_delayRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.995));
                    else
                        shf_delayRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.025));
                        shf_delayRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.005));
                    end
                    Violin(shf_delayRateLRDiffAll,1);
                    hold on
                    plot([0.7,1.3],ones(1,2).*shf_delayRateLRDiffAll95,'r-')
                    plot(1,delayRateLRDiffAll,'r*')
                    hold off
                    Fig8TreadmillArmChoice_Burst.(behaveType).SingleSpk.shf_delayRateLRDiffAll(k) = 1;
                else
                    shf_delayRateLRDiffAll95 = NaN;
                    shf_delayRateLRDiffAll99 = NaN;
                    text(0.1,0.5,'Comparison not exist','FontSize',8);
                    Fig8TreadmillArmChoice_Burst.(behaveType).SingleSpk.shf_delayRateLRDiffAll(k) = 0;
                end
                Fig8TreadmillArmChoice_Burst.(behaveType).SingleSpk.shf_delayRateLRDiffAll95(k) = shf_delayRateLRDiffAll95;
                Fig8TreadmillArmChoice_Burst.(behaveType).SingleSpk.shf_delayRateLRDiffAll99(k) = shf_delayRateLRDiffAll99;
                TITLE1 = sprintf('%s%d%s%d%s%s','Arm Choice-SingleSpk-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k}(1:end-2));
                TITLE2 = behaveType;
                title({TITLE1;TITLE2},'Interpreter','None')
                ylabel('Delay (L vs R)')
                xticklabels({''});
                yticklabels({''});
                
                % S(L) vs S(R)
                subplot(2,4,4+sesGroup)
                tempSort = sort(shf_stemRateLRDiffAll(~isnan(shf_stemRateLRDiffAll)));
                if length(tempSort) > shuffleTimes/10
                    if stemRateLRDiffAll >= 0
                        shf_stemRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.975));
                        shf_stemRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.995));
                    else
                        shf_stemRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.025));
                        shf_stemRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.005));
                    end
                    Violin(shf_stemRateLRDiffAll,1);
                    hold on
                    plot([0.7,1.3],ones(1,2).*shf_stemRateLRDiffAll95,'r-')
                    plot(1,stemRateLRDiffAll,'r*')
                    hold off
                    Fig8TreadmillArmChoice_Burst.(behaveType).SingleSpk.shf_stemRateLRDiffAll(k) = 1;
                else
                    shf_stemRateLRDiffAll95 = NaN;
                    shf_stemRateLRDiffAll99 = NaN;
                    text(0.1,0.5,'Comparison not exist','FontSize',8);
                    Fig8TreadmillArmChoice_Burst.(behaveType).SingleSpk.shf_stemRateLRDiffAll(k) = 0;
                end
                Fig8TreadmillArmChoice_Burst.(behaveType).SingleSpk.shf_stemRateLRDiffAll95(k) = shf_stemRateLRDiffAll95;
                Fig8TreadmillArmChoice_Burst.(behaveType).SingleSpk.shf_stemRateLRDiffAll99(k) = shf_stemRateLRDiffAll99;
                ylabel('Stem (L vs R)')
                xticklabels({''});
                yticklabels({''});
 
            end       
        end
    end
    
    % save figure and save .mat file
    if p.savePlot
        for k = 1:clusterNum
            figure(k)
            figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',TList{k}(1:end-2),'-ArmChoice_SingleSpk');
            print(figName,'-dpng','-r300');
        end
    end
    close all
    
    
    if p.writeToFile
        save(fullfile(savedir2,'Fig8TreadmillArmChoice_Burst.mat'), 'Fig8TreadmillArmChoice_Burst');
    end
    clear Fig8TreadmillArmChoice_Burst map_1D
    fprintf('Finished Arm choice analysis for session %d\n',i);
end

end
