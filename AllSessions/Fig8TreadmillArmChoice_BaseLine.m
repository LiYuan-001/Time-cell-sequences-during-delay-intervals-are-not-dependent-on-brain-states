function Fig8TreadmillArmChoice_BaseLine(inFile,AnalyzeSes,shuffleTimes,shuffleTime2)

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
    Fig8TreadmillArmChoice_Baseline.session = sessInfo(i).mainDir;
    Fig8TreadmillArmChoice_Baseline.rat = sessInfo(i).animal;
    Fig8TreadmillArmChoice_Baseline.day = sessInfo(i).day;
    
%     if p.savePlot
%         % directory for plot figures
%         % generate a folder for each rat eah day under the current folder
%         savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\ArmChoice');
%         if ~exist(savedir, 'dir')
%             mkdir(savedir);
%         end
%     end
    
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
    Fig8TreadmillArmChoice_Baseline.TList = TList;
    
    if prod(~isempty(char(sessDirs{1})))
        
        for sesGroup = 1:4
            sesNumSum = 0;
            posLabel = [];
            turnLabel = [];
            time_All = [];
            regionStart_All = [];
            regionEnd_All = [];

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
                        regionStart_All(sesNumSum+m,1) = PathZone.posStartT.Return(m);
                        regionEnd_All(sesNumSum+m,1) = PathZone.posEndT.Return(m);
                        
                        % delay
                        timeTemp = delayTend1_2(m)-delayTstart1(m);
                        time_All(sesNumSum+m,2) = timeTemp;
                        regionStart_All(sesNumSum+m,2) = delayTstart1(m);
                        regionEnd_All(sesNumSum+m,2) = delayTend1_2(m);
                        
                        % stem
                        timeTemp = PathZone.posEndT.Center(m)-(Fig8DelayZonePos.delayPos1.endT(m)+1/30);
                        time_All(sesNumSum+m,3) = timeTemp;     
                        regionStart_All(sesNumSum+m,3) = Fig8DelayZonePos.delayPos1.endT(m)+1/30;
                        regionEnd_All(sesNumSum+m,3) = PathZone.posEndT.Center(m);
                        
                        % choice
                        timeTemp = PathZone.posEndT.Choice(m)-PathZone.posStartT.Choice(m);
                        time_All(sesNumSum+m,4) = timeTemp;
                        regionStart_All(sesNumSum+m,4) = PathZone.posStartT.Choice(m);
                        regionEnd_All(sesNumSum+m,4) = PathZone.posEndT.Choice(m);
                        
                    end
                    sesNumSum = sesNumSum + sum(map_1D{j}.ratesByECLR.valid);
                end               
                
            end
            
            % find left / right return arm
            leftReturnIdx = find(posLabel==1); % idx on all trials
            rightReturnIdx = find(posLabel==2); % idx on all trials        
            % find Left / Right turn trials
            leftTurnIdx = find(turnLabel==1); % idx on all trials
            rightTurnIdx = find(turnLabel==0); % idx on all trials
   
            
            for k = 1:clusterNum                                   
                sessionTemp = 1;
                % initiate variables
                % return, delay, stem
                spike_All = nan(sesNumSum,4);
                spike_All_Random = nan(shuffleTime2,sesNumSum,4);
                for j = [sesGroup,sesGroup+4]
                    % load path
                    pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'PathZone.mat');
                    load(pathZoneFile);
                    delayFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
                    load(delayFile);
                    pathDataFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'pathdata.mat');
                    pathData = load(pathDataFile);
                    
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
                    tSp = Spike_Session.(sessDirs{j}){k};
                    
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
                        % choice
                        spikeCount = sum(tSp>(PathZone.posEndT.Center(m)) & tSp<PathZone.posEndT.Choice(m));
                        spike_All(sessionTemp,4) = spikeCount;
                        
                        sessionTemp = sessionTemp+1;
                    end             
                end
                
                % generate shuffleTime2 times randome distributed
                % shuffle spike set
                spike_AllSum = sum(spike_All,1);
                
                for nn = 1:length(spike_AllSum)
                    steps = regionEnd_All(:,nn) - regionStart_All(:,nn);
                    cum_steps = cumsum(steps); % Cumulative widths 
                    
                    for ii = 1:shuffleTime2
                        tsp_temp = rand(1,spike_AllSum(nn))*cum_steps(end);
                        
                        cum_steps2 = [0;cum_steps];
                        tsp_temp2 = tsp_temp;
                        for kk = 1:length(cum_steps2)-1
                            ind_temp = tsp_temp>=cum_steps2(kk) & tsp_temp<cum_steps2(kk+1);
                            tsp_temp2(ind_temp) = tsp_temp2(ind_temp) + cum_steps2(kk);
                            spike_All_Random(ii,kk,nn) = sum(ind_temp);
                        end
                    end
                    
                end
                
%                 % LR comparison
%                 rateL = nansum(nansum(spike_All(leftIdx,1)))./sum(sum(time_All(leftIdx,1)));
%                 rateR= nansum(nansum(spike_All(rightIdx,1)))./sum(sum(time_All(rightIdx,1)));
%                 LRrateDiff = (rateL-rateR)/(rateL+rateR);
%                 
                % delay zone (LR turning comparison)
                returnRateLAll = nansum(spike_All_Random(:,leftReturnIdx,1),2)./sum(time_All(leftReturnIdx,1));
                returnRateRAll = nansum(spike_All_Random(:,rightReturnIdx,1),2)./sum(time_All(rightReturnIdx,1));
                returnRateLRDiffAll = (returnRateLAll-returnRateRAll)./(returnRateLAll+returnRateRAll);
                
                % delay zone (LR turning comparison)
                delayRateLAll = nansum(spike_All_Random(:,leftTurnIdx,2),2)./sum(time_All(leftTurnIdx,2));
                delayRateRAll = nansum(spike_All_Random(:,rightTurnIdx,2),2)./sum(time_All(rightTurnIdx,2));
                delayRateLRDiffAll = (delayRateLAll-delayRateRAll)./(delayRateLAll+delayRateRAll);
                
                % stem zone (LR turning comparison)
                stemRateLAll = nansum(spike_All_Random(:,leftTurnIdx,3),2)./sum(time_All(leftTurnIdx,3));
                stemRateRAll = nansum(spike_All_Random(:,rightTurnIdx,3),2)./sum(time_All(rightTurnIdx,3));
                stemRateLRDiffAll = (stemRateLAll-stemRateRAll)./(stemRateLAll+stemRateRAll);
                
                % choice
                ChoRateLAll = nansum(spike_All_Random(:,leftTurnIdx,4),2)./sum(time_All(leftTurnIdx,4));
                ChoRateRAll = nansum(spike_All_Random(:,rightTurnIdx,4),2)./sum(time_All(rightTurnIdx,4));
                ChoRateLRDiffAll = (ChoRateLAll-ChoRateRAll)./(ChoRateLAll+ChoRateRAll);
     

%                 Fig8TreadmillCenterRateDiff.(behaveType).LRrateDiff(k) = LRrateDiff;
                Fig8TreadmillArmChoice_Baseline.(behaveType).returnRateLRDiffAll(:,k) = returnRateLRDiffAll;
                Fig8TreadmillArmChoice_Baseline.(behaveType).stemRateLRDiffAll(:,k) = stemRateLRDiffAll;
                Fig8TreadmillArmChoice_Baseline.(behaveType).delayRateLRDiffAll(:,k) = delayRateLRDiffAll;
                Fig8TreadmillArmChoice_Baseline.(behaveType).ChoRateLRDiffAll(:,k) = ChoRateLRDiffAll;
                
                % initiate shuffle quantification               
%                 shf_LRrateDiffAll = nan(shuffleTimes,1);

                for ii = 1:shuffleTime2
                    
                    shf_delayRateLRDiffAll = nan(shuffleTimes,1);
                    shf_stemRateLRDiffAll = nan(shuffleTimes,1);
                    shf_returnRateLRDiffAll = nan(shuffleTimes,1);
                    shf_ChoRateLRDiffAll = nan(shuffleTimes,1);


                    for n = 1:shuffleTimes

                        % all the followings are switch left & right
                        % all trial shulffle
                        % randi 0 = left, 1 = right
                        shfOrder = logical(randi([0,1],[sesNumSum,1]));
                        %                     shf_LrateAll = nansum(nansum(spike_All(~shfOrder,1)))./sum(sum(time_All(~shfOrder,1)));
                        %                     shf_RrateAll = nansum(nansum(spike_All(shfOrder,1)))./sum(sum(time_All(shfOrder,1)));
                        %                     shf_LRrateDiffAll(n) = (shf_LrateAll-shf_RrateAll)/(shf_LrateAll+shf_RrateAll);
                        %
                        shf_returnRateLall = nansum(spike_All_Random(ii,~shfOrder,1))./nansum(time_All(~shfOrder,1));
                        shf_returnRateRall =  nansum(spike_All_Random(ii,shfOrder,1))./nansum(time_All(shfOrder,1));
                        shf_returnRateLRDiffAll(n) = (shf_returnRateLall-shf_returnRateRall)./(shf_returnRateLall+shf_returnRateRall);

                        shf_delayRateLall =  nansum(spike_All_Random(ii,~shfOrder,2))./nansum(time_All(~shfOrder,2));
                        shf_delayRateRall =  nansum(spike_All_Random(ii,shfOrder,2))./nansum(time_All(shfOrder,2));
                        shf_delayRateLRDiffAll(n) = (shf_delayRateLall-shf_delayRateRall)./(shf_delayRateLall+shf_delayRateRall);

                        shf_stemRateLall =  nansum(spike_All_Random(ii,~shfOrder,3))./nansum(time_All(~shfOrder,3));
                        shf_stemRateRall =  nansum(spike_All_Random(ii,shfOrder,3))./nansum(time_All(shfOrder,3));
                        shf_stemRateLRDiffAll(n) = (shf_stemRateLall-shf_stemRateRall)./(shf_stemRateLall+shf_stemRateRall);

                        shf_ChoRateLall =  nansum(spike_All_Random(ii,~shfOrder,4))./nansum(time_All(~shfOrder,4));
                        shf_ChoRateRall =  nansum(spike_All_Random(ii,shfOrder,4))./nansum(time_All(shfOrder,4));
                        shf_ChoRateLRDiffAll(n) = (shf_ChoRateLall-shf_ChoRateRall)./(shf_ChoRateLall+shf_ChoRateRall);
                        
                    end


                    % D(L) vs D(R)
                    tempSort = sort(shf_delayRateLRDiffAll(~isnan(shf_delayRateLRDiffAll)));
                    if length(tempSort) > shuffleTimes/10
                        if delayRateLRDiffAll(ii) >= 0
                            shf_delayRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.975));
                            shf_delayRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.995));
                        else
                            shf_delayRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.025));
                            shf_delayRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.005));
                        end
                        Fig8TreadmillArmChoice_Baseline.(behaveType).shf_delayRateLRDiffAll(ii,k) = 1;
                    else
                        shf_delayRateLRDiffAll95 = NaN;
                        shf_delayRateLRDiffAll99 = NaN;
                        Fig8TreadmillArmChoice_Baseline.(behaveType).shf_delayRateLRDiffAll(ii,k) = 0;
                    end
                    Fig8TreadmillArmChoice_Baseline.(behaveType).shf_delayRateLRDiffAll95(ii,k) = shf_delayRateLRDiffAll95;
                    Fig8TreadmillArmChoice_Baseline.(behaveType).shf_delayRateLRDiffAll99(ii,k) = shf_delayRateLRDiffAll99;

                    % Return(L) vs Return(R)
                    tempSort = sort(shf_returnRateLRDiffAll(~isnan(shf_returnRateLRDiffAll)));
                    if length(tempSort) > shuffleTimes/10
                        if returnRateLRDiffAll(ii) >= 0
                            shf_returnRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.95));
                            shf_returnRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.99));
                        else
                            shf_returnRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.05));
                            shf_returnRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.01));
                        end

                        Fig8TreadmillArmChoice_Baseline.(behaveType).shf_returnRateLRDiffAll(ii,k) = 1;
                    else
                        shf_returnRateLRDiffAll95 = NaN;
                        shf_returnRateLRDiffAll99 = NaN;
                        Fig8TreadmillArmChoice_Baseline.(behaveType).shf_returnRateLRDiffAll(ii,k) = 0;
                    end
                    Fig8TreadmillArmChoice_Baseline.(behaveType).shf_returnRateLRDiffAll95(ii,k) = shf_returnRateLRDiffAll95;
                    Fig8TreadmillArmChoice_Baseline.(behaveType).shf_returnRateLRDiffAll99(ii,k) = shf_returnRateLRDiffAll99;


                    % S(L) vs S(R)
                    tempSort = sort(shf_stemRateLRDiffAll(~isnan(shf_stemRateLRDiffAll)));
                    if length(tempSort) > shuffleTimes/10
                        if stemRateLRDiffAll(ii) >= 0
                            shf_stemRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.975));
                            shf_stemRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.995));
                        else
                            shf_stemRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.025));
                            shf_stemRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.005));
                        end
                        Fig8TreadmillArmChoice_Baseline.(behaveType).shf_stemRateLRDiffAll(ii,k) = 1;
                    else
                        shf_stemRateLRDiffAll95 = NaN;
                        shf_stemRateLRDiffAll99 = NaN;
                        Fig8TreadmillArmChoice_Baseline.(behaveType).shf_stemRateLRDiffAll(ii,k) = 0;
                    end
                    Fig8TreadmillArmChoice_Baseline.(behaveType).shf_stemRateLRDiffAll95(ii,k) = shf_stemRateLRDiffAll95;
                    Fig8TreadmillArmChoice_Baseline.(behaveType).shf_stemRateLRDiffAll99(ii,k) = shf_stemRateLRDiffAll99;


                    % Choice (L) vs (R)
                    tempSort = sort(shf_ChoRateLRDiffAll(~isnan(shf_ChoRateLRDiffAll)));
                    if length(tempSort) > shuffleTimes/10
                        if ChoRateLRDiffAll(ii) >= 0
                            shf_ChoRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.95));
                            shf_ChoRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.99));
                        else
                            shf_ChoRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.05));
                            shf_ChoRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.01));
                        end
                        Fig8TreadmillArmChoice_Baseline.(behaveType).shf_ChoRateLRDiffAll(ii,k) = 1;
                    else
                        shf_ChoRateLRDiffAll95 = NaN;
                        shf_ChoRateLRDiffAll99 = NaN;
                        Fig8TreadmillArmChoice_Baseline.(behaveType).shf_ChoRateLRDiffAll(ii,k) = 0;
                    end
                    Fig8TreadmillArmChoice_Baseline.(behaveType).shf_ChoRateLRDiffAll95(ii,k) = shf_ChoRateLRDiffAll95;
                    Fig8TreadmillArmChoice_Baseline.(behaveType).shf_ChoRateLRDiffAll99(ii,k) = shf_ChoRateLRDiffAll99;
                end
            end
        end
    end
%     % save figure and save .mat file
%     if p.savePlot
%         for k = 1:clusterNum
%             figure(k)
%             figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',TList{k}(1:end-2),'-ArmChoice');
%             print(figName,'-dpng','-r300');
%         end
%     end
    
    if p.writeToFile
        save(fullfile(savedir2,'Fig8TreadmillArmChoice_Baseline.mat'), 'Fig8TreadmillArmChoice_Baseline');
    end
    clear Fig8TreadmillArmChoice_Baseline map_1D
    close all
    fprintf('Finished Arm choice analysis for session %d\n',i);
end

end
