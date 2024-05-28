function Fig8TreadmillArmChoice_EarlyLate(inFile,AnalyzeSes,shuffleTimes)

%  Li Yuan, UCSD, 22-Mar-2021
%  Calculate arm selectivity based on average rate
% -------------------------------------------------------------------------
% set parameters for analysis

p.savePlot = 0;
p.writeToFile = 1;
% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'

% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    close all
    % selectivity for the whole map
    Fig8TreadmillArmChoice.session = sessInfo(i).mainDir;
    Fig8TreadmillArmChoice.rat = sessInfo(i).animal;
    Fig8TreadmillArmChoice.day = sessInfo(i).day;
    
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\ArmChoice');
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
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);
    Fig8TreadmillArmChoice.TList = TList;
    
    if prod(~isempty(char(sessDirs{1})))           
        for j = 1:length(sessDirs)
            
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
            
            posLabel = zeros(trialNum,1);
            turnLabel = zeros(trialNum,1);
            time_All = zeros(trialNum,1);
            
            if exist(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'), 'file')==2
                map_1D{j} = load(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'));
                
                for m = 1:trialNum
                    if map_1D{j}.ratesByECLR.ECLR(m) == 1
                        posLabel(m) = 2;
                        turnLabel(m) = 1;
                    elseif map_1D{j}.ratesByECLR.ECLR(m) == 2
                        posLabel(m) = 1;
                        turnLabel(m) = 1;
                    elseif map_1D{j}.ratesByECLR.ECLR(m) == 3
                        posLabel(m) = 1;
                        turnLabel(m) = 0;
                    elseif map_1D{j}.ratesByECLR.ECLR(m) == 4
                        posLabel(m) = 2;
                        turnLabel(m) = 0;
                    else
                        error('Turning label ERROR');
                    end
                    
                    
                    % return
                    timeTemp = PathZone.posEndT.Return(m)-PathZone.posStartT.Return(m);
                    time_All(m,1) = timeTemp;
                    % delay
                    timeTemp = delayTend1_2(m)-delayTstart1(m);
                    time_All(m,2) = timeTemp;
                    % stem
                    timeTemp = PathZone.posEndT.Center(m)-(Fig8DelayZonePos.delayPos1.endT(m)+1/30);
                    time_All(m,3) = timeTemp;
                    % choice and reward
                    timeTemp = PathZone.posEndT.Choice(m)-PathZone.posStartT.Choice(m);
                    time_All(m,4) = timeTemp;
                end
                % find left / right return arm
                leftReturnIdx = find(posLabel==1); % idx on all trials
                rightReturnIdx = find(posLabel==2); % idx on all trials
                % find Left / Right turn trials
                leftTurnIdx = find(turnLabel==1); % idx on all trials
                rightTurnIdx = find(turnLabel==0); % idx on all trials
                
                for k = 1:clusterNum
                    % initiate variables
                    % return, delay, stem
                    spike_All = nan(trialNum,4);
                    % get spikes
                    tSp = Spike_Session.(sessDirs{j}){k};
                    
                    for m = 1:sum(map_1D{j}.ratesByECLR.valid)
                        % get spike count
                        % return
                        spikeCount = sum(tSp>=PathZone.posStartT.Return(m) & tSp<PathZone.posEndT.Return(m));
                        spike_All(m,1) = spikeCount;
                        % delay
                        spikeCount = sum(tSp>=delayTstart1(m) & tSp<delayTend1_2(m));
                        spike_All(m,2) = spikeCount;
                        % stem
                        spikeCount = sum(tSp>=(Fig8DelayZonePos.delayPos1.endT(m)+1/30) & tSp<PathZone.posEndT.Center(m));
                        spike_All(m,3) = spikeCount;
                        % choice
                        spikeCount = sum(tSp>(PathZone.posEndT.Center(m)) & tSp<PathZone.posEndT.Choice(m));
                        spike_All(m,4) = spikeCount;
                        
                    end
                    
                    %                 % LR comparison
                    %                 rateL = nansum(nansum(spike_All(leftIdx,1)))./sum(sum(time_All(leftIdx,1)));
                    %                 rateR= nansum(nansum(spike_All(rightIdx,1)))./sum(sum(time_All(rightIdx,1)));
                    %                 LRrateDiff = (rateL-rateR)/(rateL+rateR);
                    %
                    % delay zone (LR turning comparison)
                    returnRateLAll = nansum(nansum(spike_All(leftReturnIdx,1)))./sum(sum(time_All(leftReturnIdx,1)));
                    returnRateRAll = nansum(nansum(spike_All(rightReturnIdx,1)))./sum(sum(time_All(rightReturnIdx,1)));
                    returnRateLRDiffAll = (returnRateLAll-returnRateRAll)./(returnRateLAll+returnRateRAll);
                    
                    % delay zone (LR turning comparison)
                    delayRateLAll = nansum(nansum(spike_All(leftTurnIdx,2)))./sum(sum(time_All(leftTurnIdx,2)));
                    delayRateRAll = nansum(nansum(spike_All(rightTurnIdx,2)))./sum(sum(time_All(rightTurnIdx,2)));
                    delayRateLRDiffAll = (delayRateLAll-delayRateRAll)/(delayRateLAll+delayRateRAll);
                    
                    % stem zone (LR turning comparison)
                    stemRateLAll = nansum(nansum(spike_All(leftTurnIdx,3)))./sum(sum(time_All(leftTurnIdx,3)));
                    stemRateRAll = nansum(nansum(spike_All(rightTurnIdx,3)))./sum(sum(time_All(rightTurnIdx,3)));
                    stemRateLRDiffAll = (stemRateLAll-stemRateRAll)/(stemRateLAll+stemRateRAll);
                    
                    % choice
                    ChoRateLAll = nansum(nansum(spike_All(leftTurnIdx,4)))./sum(sum(time_All(leftTurnIdx,4)));
                    ChoRateRAll = nansum(nansum(spike_All(rightTurnIdx,4)))./sum(sum(time_All(rightTurnIdx,4)));
                    ChoRateLRDiffAll = (ChoRateLAll-ChoRateRAll)./(ChoRateLAll+ChoRateRAll);
                    
                    
                    %                 Fig8TreadmillCenterRateDiff.(sessDirs{j}).LRrateDiff(k) = LRrateDiff;
                    Fig8TreadmillArmChoice.(sessDirs{j}).returnRateLRDiffAll(k) = returnRateLRDiffAll;
                    Fig8TreadmillArmChoice.(sessDirs{j}).stemRateLRDiffAll(k) = stemRateLRDiffAll;
                    Fig8TreadmillArmChoice.(sessDirs{j}).delayRateLRDiffAll(k) = delayRateLRDiffAll;
                    Fig8TreadmillArmChoice.(sessDirs{j}).ChoRateLRDiffAll(k) = ChoRateLRDiffAll;
                    
                    % initiate shuffle quantification
                    %                 shf_LRrateDiffAll = nan(shuffleTimes,1);
                    shf_delayRateLRDiffAll = nan(shuffleTimes,1);
                    shf_stemRateLRDiffAll = nan(shuffleTimes,1);
                    shf_returnRateLRDiffAll = nan(shuffleTimes,1);
                    shf_ChoRateLRDiffAll = nan(shuffleTimes,1);
                    
                    for n = 1:shuffleTimes
                        
                        % all the followings are switch left & right
                        % all trial shulffle
                        % randi 0 = left, 1 = right
                        shfOrder = logical(randi([0,1],[trialNum,1]));
                        %                     shf_LrateAll = nansum(nansum(spike_All(~shfOrder,1)))./sum(sum(time_All(~shfOrder,1)));
                        %                     shf_RrateAll = nansum(nansum(spike_All(shfOrder,1)))./sum(sum(time_All(shfOrder,1)));
                        %                     shf_LRrateDiffAll(n) = (shf_LrateAll-shf_RrateAll)/(shf_LrateAll+shf_RrateAll);
                        %
                        shf_returnRateLall = nansum(nansum(spike_All(~shfOrder,1)))./nansum(nansum(time_All(~shfOrder,1)));
                        shf_returnRateRall = nansum(nansum(spike_All(shfOrder,1)))./nansum(nansum(time_All(shfOrder,1)));
                        shf_returnRateLRDiffAll(n) = (shf_returnRateLall-shf_returnRateRall)./(shf_returnRateLall+shf_returnRateRall);
                        
                        shf_delayRateLall = nansum(nansum(spike_All(~shfOrder,2)))./nansum(nansum(time_All(~shfOrder,2)));
                        shf_delayRateRall = nansum(nansum(spike_All(shfOrder,2)))./nansum(nansum(time_All(shfOrder,2)));
                        shf_delayRateLRDiffAll(n) = (shf_delayRateLall-shf_delayRateRall)/(shf_delayRateLall+shf_delayRateRall);
                        
                        shf_stemRateLall = nansum(nansum(spike_All(~shfOrder,3)))./nansum(nansum(time_All(~shfOrder,3)));
                        shf_stemRateRall = nansum(nansum(spike_All(shfOrder,3)))./nansum(nansum(time_All(shfOrder,3)));
                        shf_stemRateLRDiffAll(n) = (shf_stemRateLall-shf_stemRateRall)/(shf_stemRateLall+shf_stemRateRall);
                        
                        shf_ChoRateLall = nansum(nansum(spike_All(~shfOrder,4)))./nansum(nansum(time_All(~shfOrder,4)));
                        shf_ChoRateRall = nansum(nansum(spike_All(shfOrder,4)))./nansum(nansum(time_All(shfOrder,4)));
                        shf_ChoRateLRDiffAll(n) = (shf_ChoRateLall-shf_ChoRateRall)./(shf_ChoRateLall+shf_ChoRateRall);
                        
                    end
                    
                    h = figure(k);
                    if j == 1
                        h.Position = [100 100 1800 1200];
                    end
                    
                    % D(L) vs D(R)
                    subplot(4,8,j)
                    tempSort = sort(shf_delayRateLRDiffAll(~isnan(shf_delayRateLRDiffAll)));
                    if length(tempSort) > shuffleTimes/10
                        if delayRateLRDiffAll >= 0
                            shf_delayRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.975));
%                             shf_delayRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.995));
                        else
                            shf_delayRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.025));
%                             shf_delayRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.005));
                        end
                        Violin(shf_delayRateLRDiffAll,1);
                        hold on
                        plot([0.7,1.3],ones(1,2).*shf_delayRateLRDiffAll95,'r-')
                        plot(1,delayRateLRDiffAll,'r*')
                        hold off
                        Fig8TreadmillArmChoice.(sessDirs{j}).shf_delayRateLRDiffAll(k) = 1;
                    else
                        shf_delayRateLRDiffAll95 = NaN;
                        shf_delayRateLRDiffAll99 = NaN;
                        text(0.1,0.5,'Comparison not exist','FontSize',8);
                        Fig8TreadmillArmChoice.(sessDirs{j}).shf_delayRateLRDiffAll(k) = 0;
                    end
                    Fig8TreadmillArmChoice.(sessDirs{j}).shf_delayRateLRDiffAll95(k) = shf_delayRateLRDiffAll95;
%                     Fig8TreadmillArmChoice.(sessDirs{j}).shf_delayRateLRDiffAll99(k) = shf_delayRateLRDiffAll99;
                    TITLE1 = sprintf('%s%d%s%d%s%s','Arm Choice-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k}(1:end-2));
                    TITLE2 = sessDirs{j};
                    title({TITLE1;TITLE2},'Interpreter','None')
                    ylabel('Delay (L vs R)')
                    xticklabels({''});
                    yticklabels({''});
                    
                    % Return(L) vs Return(R)
                    subplot(4,8,8+j)
                    tempSort = sort(shf_returnRateLRDiffAll(~isnan(shf_returnRateLRDiffAll)));
                    if length(tempSort) > shuffleTimes/10
                        if returnRateLRDiffAll >= 0
                            shf_returnRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.95));
%                             shf_returnRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.99));
                        else
                            shf_returnRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.05));
%                             shf_returnRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.01));
                        end
                        Violin(shf_returnRateLRDiffAll,1);
                        hold on
                        plot([0.7,1.3],ones(1,2).*shf_returnRateLRDiffAll95,'r-')
                        plot(1,returnRateLRDiffAll,'r*')
                        hold off
                        Fig8TreadmillArmChoice.(sessDirs{j}).shf_returnRateLRDiffAll(k) = 1;
                    else
                        shf_returnRateLRDiffAll95 = NaN;
                        shf_returnRateLRDiffAll99 = NaN;
                        text(0.1,0.5,'Comparison not exist','FontSize',8);
                        Fig8TreadmillArmChoice.(sessDirs{j}).shf_returnRateLRDiffAll(k) = 0;
                    end
                    Fig8TreadmillArmChoice.(sessDirs{j}).shf_returnRateLRDiffAll95(k) = shf_returnRateLRDiffAll95;
%                     Fig8TreadmillArmChoice.(sessDirs{j}).shf_returnRateLRDiffAll99(k) = shf_returnRateLRDiffAll99;
                    ylabel('Return+Base (L vs R)')
                    
                    
                    % S(L) vs S(R)
                    subplot(4,8,16+j)
                    tempSort = sort(shf_stemRateLRDiffAll(~isnan(shf_stemRateLRDiffAll)));
                    if length(tempSort) > shuffleTimes/10
                        if stemRateLRDiffAll >= 0
                            shf_stemRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.975));
%                             shf_stemRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.995));
                        else
                            shf_stemRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.025));
%                             shf_stemRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.005));
                        end
                        Violin(shf_stemRateLRDiffAll,1);
                        hold on
                        plot([0.7,1.3],ones(1,2).*shf_stemRateLRDiffAll95,'r-')
                        plot(1,stemRateLRDiffAll,'r*')
                        hold off
                        Fig8TreadmillArmChoice.(sessDirs{j}).shf_stemRateLRDiffAll(k) = 1;
                    else
                        shf_stemRateLRDiffAll95 = NaN;
                        shf_stemRateLRDiffAll99 = NaN;
                        text(0.1,0.5,'Comparison not exist','FontSize',8);
                        Fig8TreadmillArmChoice.(sessDirs{j}).shf_stemRateLRDiffAll(k) = 0;
                    end
                    Fig8TreadmillArmChoice.(sessDirs{j}).shf_stemRateLRDiffAll95(k) = shf_stemRateLRDiffAll95;
%                     Fig8TreadmillArmChoice.(sessDirs{j}).shf_stemRateLRDiffAll99(k) = shf_stemRateLRDiffAll99;
                    ylabel('Stem (L vs R)')
                    xticklabels({''});
                    yticklabels({''});
                    
                    
                    % Choice (L) vs (R)
                    subplot(4,8,24+j)
                    tempSort = sort(shf_ChoRateLRDiffAll(~isnan(shf_ChoRateLRDiffAll)));
                    if length(tempSort) > shuffleTimes/10
                        if ChoRateLRDiffAll >= 0
                            shf_ChoRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.95));
%                             shf_ChoRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.99));
                        else
                            shf_ChoRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.05));
%                             shf_ChoRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.01));
                        end
                        Violin(shf_ChoRateLRDiffAll,1);
                        hold on
                        plot([0.7,1.3],ones(1,2).*shf_ChoRateLRDiffAll95,'r-')
                        plot(1,ChoRateLRDiffAll,'r*')
                        hold off
                        Fig8TreadmillArmChoice.(sessDirs{j}).shf_ChoRateLRDiffAll(k) = 1;
                    else
                        shf_ChoRateLRDiffAll95 = NaN;
                        shf_ChoRateLRDiffAll99 = NaN;
                        text(0.1,0.5,'Comparison not exist','FontSize',8);
                        Fig8TreadmillArmChoice.(sessDirs{j}).shf_ChoRateLRDiffAll(k) = 0;
                    end
                    Fig8TreadmillArmChoice.(sessDirs{j}).shf_ChoRateLRDiffAll95(k) = shf_ChoRateLRDiffAll95;
%                     Fig8TreadmillArmChoice.(sessDirs{j}).shf_ChoRateLRDiffAll99(k) = shf_ChoRateLRDiffAll99;
                    ylabel('Choice (L vs R)')
                    
                    
                end
                
            end
            
        end
    end
    
    % save figure and save .mat file
    if p.savePlot
        for k = 1:clusterNum
            figure(k)
            figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',TList{k}(1:end-2),'-ArmChoice');
            print(figName,'-dpng','-r300');
        end
    end
    
    if p.writeToFile
        save(fullfile(savedir2,'Fig8TreadmillArmChoice.mat'), 'Fig8TreadmillArmChoice');
    end
    clear Fig8TreadmillArmChoice map_1D
    close all
    fprintf('Finished Arm choice analysis for session %d\n',i);
end

end
