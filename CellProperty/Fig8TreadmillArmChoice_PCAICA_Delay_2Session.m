%  Li Yuan, UCSD, 10-May-2022
%  Calculate arm selectivity based on average rate
% -------------------------------------------------------------------------
function Fig8TreadmillArmChoice_PCAICA_Delay_2Session(inFile,AnalyzeSes,shuffleTimes)

% set parameters for analysis

p.savePlot = 1;
p.writeToFile = 1;

% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    close all
    % selectivity for the whole map
    AssemblyDelay_ArmChoice_SameDelay.session = sessInfo(i).mainDir;
    AssemblyDelay_ArmChoice_SameDelay.rat = sessInfo(i).animal;
    AssemblyDelay_ArmChoice_SameDelay.day = sessInfo(i).day;
    
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\Assembly_Delay ArmChoice');
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
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load assembly file
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_Delay-25ms.mat');
    load(assemblyFile);
    
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;    
    tList = SpikeProp.tList(rateLabel);
    
    patNum = CellAssembly_Delay.patNum ;
    AssmblPtrnCellIDs = CellAssembly_Delay.AssmblPtrnCellIDs;
    AssmblWght = CellAssembly_Delay.AssmblWght;
    AssmblStrength = CellAssembly_Delay.AssmblStrength;
    event_Time = CellAssembly_Delay.event_Time;
    event_strength = CellAssembly_Delay.event_strength;
    event_Num = CellAssembly_Delay.event_Num;
    AssemblyDelay_ArmChoice_SameDelay.TList = tList;
    
    if prod(~isempty(char(sessDirs{1})))
        
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
                if strcmp(sessDirs{j}(end-3:end-2),'10')
                    maxT = 10;
                    delayTend1_2 = delayTstart1+maxT;
                elseif strcmp(sessDirs{j}(end-3:end-2),'30')
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
                        % choice and reward
                        timeTemp = PathZone.posEndT.Choice(m)-PathZone.posStartT.Choice(m);
                        time_All(sesNumSum+m,4) = timeTemp;
                    end
                    sesNumSum = sesNumSum + sum(map_1D{j}.ratesByECLR.valid);
                end               
                
            end
            
            
            % find Left / Right turn trials
            leftTurnIdx = find(turnLabel==1); % idx on all trials
            rightTurnIdx = find(turnLabel==0); % idx on all trials
            % find left / right return arm
            leftReturnIdx = find(posLabel==2); % idx on all trials
            rightReturnIdx = find(posLabel==1); % idx on all trials
            
            for k = 1:patNum                                   
                sessionTemp = 1;
                % initiate variables
                % return, delay, stem
                spike_All = nan(sesNumSum,4);
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
                     if strcmp(sessDirs{j}(end-3:end-2),'10')
                         maxT = 10;
                         delayTend1_2 = delayTstart1+maxT;
                     elseif strcmp(sessDirs{j}(end-3:end-2),'30')
                         maxT = 30;
                         delayTend1_2 = delayTstart1+maxT;
                     else
                         error('Delay time is wrong')
                     end
                
                    % get each spike time
                    tSp = event_Time{k};
                    strength = event_strength{k};
                    
                    for m = 1:sum(map_1D{j}.ratesByECLR.valid)
                        % get spike count
                        % return
                        spikeInd = (tSp>=PathZone.posStartT.Return(m) & tSp<PathZone.posEndT.Return(m));
                        spike_All(sessionTemp,1) = sum(strength(spikeInd));
                        % delay
                        spikeInd = (tSp>=delayTstart1(m) & tSp<delayTend1_2(m));
                        spike_All(sessionTemp,2) = sum(strength(spikeInd));
                        % stem
                        spikeInd = (tSp>=(Fig8DelayZonePos.delayPos1.endT(m)+1/30) & tSp<PathZone.posEndT.Center(m));
                        spike_All(sessionTemp,3) = sum(strength(spikeInd));                      
                        % choice
                        spikeInd = (tSp>=(PathZone.posStartT.Choice(m)) & tSp<PathZone.posEndT.Choice(m));
                        spike_All(sessionTemp,4) = sum(strength(spikeInd));
                        
                        sessionTemp = sessionTemp+1;
                    end
                end
                
%                 % LR comparison
%                 rateL = nansum(nansum(spike_All(leftTurnIdx,1)))./sum(sum(time_All(leftTurnIdx,1)));
%                 rateR= nansum(nansum(spike_All(rightTurnIdx,1)))./sum(sum(time_All(rightTurnIdx,1)));
%                 LRrateDiff = (rateL-rateR);
%                 
                % return zone (LR turning comparison)
                returnRateLAll = nansum(nansum(spike_All(leftReturnIdx,1)))./sum(sum(time_All(leftReturnIdx,1)));
                returnRateRAll = nansum(nansum(spike_All(rightReturnIdx,1)))./sum(sum(time_All(rightReturnIdx,1)));
                returnRateLRDiffAll = (returnRateLAll-returnRateRAll);
                
                % delay zone (LR turning comparison)
                delayRateLAll = nansum(nansum(spike_All(leftTurnIdx,2)))./sum(sum(time_All(leftTurnIdx,2)));
                delayRateRAll = nansum(nansum(spike_All(rightTurnIdx,2)))./sum(sum(time_All(rightTurnIdx,2)));
                delayRateLRDiffAll = (delayRateLAll-delayRateRAll);
                
                % stem zone (LR turning comparison)
                stemRateLAll = nansum(nansum(spike_All(leftTurnIdx,3)))./sum(sum(time_All(leftTurnIdx,3)));
                stemRateRAll = nansum(nansum(spike_All(rightTurnIdx,3)))./sum(sum(time_All(rightTurnIdx,3)));
                stemRateLRDiffAll = (stemRateLAll-stemRateRAll);
                
                % choice
                stemChoRateLAll = nansum(nansum(spike_All(leftTurnIdx,4)))./sum(sum(time_All(leftTurnIdx,4)));
                stemChoRateRAll = nansum(nansum(spike_All(rightTurnIdx,4)))./sum(sum(time_All(rightTurnIdx,4)));
                ChoRateLRDiffAll = (stemChoRateLAll-stemChoRateRAll);

                AssemblyDelay_ArmChoice_SameDelay.(behaveType).returnRateLRDiffAll(k) = returnRateLRDiffAll;
                AssemblyDelay_ArmChoice_SameDelay.(behaveType).delayRateLRDiffAll(k) = delayRateLRDiffAll;
                AssemblyDelay_ArmChoice_SameDelay.(behaveType).stemRateLRDiffAll(k) = stemRateLRDiffAll;
                AssemblyDelay_ArmChoice_SameDelay.(behaveType).ChoRateLRDiffAll(k) = ChoRateLRDiffAll;
                
                
                % initiate shuffle quantification               
                shf_returnRateLRDiffAll = nan(shuffleTimes,1);
                shf_delayRateLRDiffAll = nan(shuffleTimes,1);
                shf_stemRateLRDiffAll = nan(shuffleTimes,1);
                shf_ChoRateLRDiffAll = nan(shuffleTimes,1);
                
                for n = 1:shuffleTimes
                    
                    % all the followings are switch left & right
                    % all trial shulffle
                    % randi 0 = left, 1 = right
                    shfOrder = logical(randi([0,1],[sesNumSum,1]));
                    
%                     shf_LrateAll = nansum(nansum(spike_All(~shfOrder,1)))./sum(sum(time_All(~shfOrder,1)));
%                     shf_RrateAll = nansum(nansum(spike_All(shfOrder,1)))./sum(sum(time_All(shfOrder,1)));
%                     shf_LRrateDiffAll(n) = (shf_LrateAll-shf_RrateAll);
% %                     
                    shf_returnRateLall = nansum(nansum(spike_All(~shfOrder,1)))./nansum(nansum(time_All(~shfOrder,1)));
                    shf_returnRateRall = nansum(nansum(spike_All(shfOrder,1)))./nansum(nansum(time_All(shfOrder,1)));
                    shf_returnRateLRDiffAll(n) = (shf_returnRateLall-shf_returnRateRall);
                    
                    shf_delayRateLall = nansum(nansum(spike_All(~shfOrder,2)))./nansum(nansum(time_All(~shfOrder,2)));
                    shf_delayRateRall = nansum(nansum(spike_All(shfOrder,2)))./nansum(nansum(time_All(shfOrder,2)));
                    shf_delayRateLRDiffAll(n) = (shf_delayRateLall-shf_delayRateRall);
                    
                    shf_stemRateLall = nansum(nansum(spike_All(~shfOrder,3)))./nansum(nansum(time_All(~shfOrder,3)));
                    shf_stemRateRall = nansum(nansum(spike_All(shfOrder,3)))./nansum(nansum(time_All(shfOrder,3)));
                    shf_stemRateLRDiffAll(n) = (shf_stemRateLall-shf_stemRateRall);
                    
                    shf_ChoRateLall = nansum(nansum(spike_All(~shfOrder,4)))./nansum(nansum(time_All(~shfOrder,4)));
                    shf_ChoRateRall = nansum(nansum(spike_All(shfOrder,4)))./nansum(nansum(time_All(shfOrder,4)));
                    shf_ChoRateLRDiffAll(n) = (shf_ChoRateLall-shf_ChoRateRall);
                    
                end 
                
                h = figure(k);
                if sesGroup == 1
                    h.Position = [100 100 1600 9000];
                end
                
                % D(L) vs D(R)
                subplot(4,4,sesGroup)
                tempSort = sort(shf_delayRateLRDiffAll(~isnan(shf_delayRateLRDiffAll)));
                if length(tempSort) > shuffleTimes/10
                    if delayRateLRDiffAll >= 0
                        shf_delayRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.95));
                        shf_delayRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.99));
                    else
                        shf_delayRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.05));
                        shf_delayRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.01));
                    end
                    Violin(shf_delayRateLRDiffAll,1);
                    hold on
                    plot([0.7,1.3],ones(1,2).*shf_delayRateLRDiffAll95,'r-')
                    plot(1,delayRateLRDiffAll,'r*')
                    hold off
                    AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_delayRateLRDiffAll(k) = 1;
                else
                    shf_delayRateLRDiffAll95 = NaN;
                    shf_delayRateLRDiffAll99 = NaN;
                    text(0.1,0.5,'Comparison not exist','FontSize',8);
                    AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_delayRateLRDiffAll(k) = 0;
                end
                AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_delayRateLRDiffAll95(k) = shf_delayRateLRDiffAll95;
                AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_delayRateLRDiffAll99(k) = shf_delayRateLRDiffAll99;
                TITLE1 = sprintf('%s%d%s%d%s%d','Arm Choice-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-Assembly-',k);
                TITLE2 = behaveType;
                title({TITLE1;TITLE2},'Interpreter','None')
                ylabel('Delay (L vs R)')
%                 xticklabels({''});
%                 yticklabels({''});
                
                % Return(L) vs Return(R)
                subplot(4,4,4+sesGroup)
                tempSort = sort(shf_returnRateLRDiffAll(~isnan(shf_returnRateLRDiffAll)));
                if length(tempSort) > shuffleTimes/10
                    if stemRateLRDiffAll >= 0
                        shf_returnRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.95));
                        shf_returnRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.99));
                    else
                        shf_returnRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.05));
                        shf_returnRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.01));
                    end
                    Violin(shf_returnRateLRDiffAll,1);
                    hold on
                    plot([0.7,1.3],ones(1,2).*shf_returnRateLRDiffAll95,'r-')
                    plot(1,stemRateLRDiffAll,'r*')
                    hold off
                    AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_returnRateLRDiffAll(k) = 1;
                else
                    shf_returnRateLRDiffAll95 = NaN;
                    shf_returnRateLRDiffAll99 = NaN;
                    text(0.1,0.5,'Comparison not exist','FontSize',8);
                    AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_returnRateLRDiffAll(k) = 0;
                end
                AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_returnRateLRDiffAll95(k) = shf_returnRateLRDiffAll95;
                AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_returnRateLRDiffAll99(k) = shf_returnRateLRDiffAll99;
                ylabel('Return (L vs R)')
%                 xticklabels({''});
%                 yticklabels({''});
                
                % stem (L) vs (R)
                subplot(4,4,8+sesGroup)
                tempSort = sort(shf_stemRateLRDiffAll(~isnan(shf_stemRateLRDiffAll)));
                if length(tempSort) > shuffleTimes/10
                    if stemRateLRDiffAll >= 0
                        shf_stemRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.95));
                        shf_stemRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.99));
                    else
                        shf_stemRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.05));
                        shf_stemRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.01));
                    end
                    Violin(shf_stemRateLRDiffAll,1);
                    hold on
                    plot([0.7,1.3],ones(1,2).*shf_stemRateLRDiffAll95,'r-')
                    plot(1,stemRateLRDiffAll,'r*')
                    hold off
                    AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_stemRateLRDiffAll(k) = 1;
                else
                    shf_stemRateLRDiffAll95 = NaN;
                    shf_stemRateLRDiffAll99 = NaN;
                    text(0.1,0.5,'Comparison not exist','FontSize',8);
                    AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_ChoRateLRDiffAll(k) = 0;
                end
                AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_stemRateLRDiffAll95(k) = shf_stemRateLRDiffAll95;
                AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_stemRateLRDiffAll99(k) = shf_stemRateLRDiffAll99;
                ylabel('Stem (L vs R)')
%                 xticklabels({''});
%                 yticklabels({''});

                %  Choice (L) vs (R)
                subplot(4,4,12+sesGroup)
                tempSort = sort(shf_ChoRateLRDiffAll(~isnan(shf_ChoRateLRDiffAll)));
                if length(tempSort) > shuffleTimes/10
                    if stemRateLRDiffAll >= 0
                        shf_ChoRateLRDiffAll95 = tempSort(ceil(length(tempSort)*0.95));
                        shf_ChoRateLRDiffAll99 = tempSort(ceil(length(tempSort)*0.99));
                    else
                        shf_ChoRateLRDiffAll95 = tempSort(floor(length(tempSort)*0.05));
                        shf_ChoRateLRDiffAll99 = tempSort(floor(length(tempSort)*0.01));
                    end
                    Violin(shf_ChoRateLRDiffAll,1);
                    hold on
                    plot([0.7,1.3],ones(1,2).*shf_ChoRateLRDiffAll95,'r-')
                    plot(1,stemRateLRDiffAll,'r*')
                    hold off
                    AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_stemChoRateLRDiffAll(k) = 1;
                else
                    shf_ChoRateLRDiffAll95 = NaN;
                    shf_ChoRateLRDiffAll99 = NaN;
                    text(0.1,0.5,'Comparison not exist','FontSize',8);
                    AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_stemChoRateLRDiffAll(k) = 0;
                end
                AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_stemChoRateLRDiffAll95(k) = shf_ChoRateLRDiffAll95;
                AssemblyDelay_ArmChoice_SameDelay.(behaveType).shf_stemChoRateLRDiffAll99(k) = shf_ChoRateLRDiffAll99;
                ylabel('Choice (L vs R)')
%                 xticklabels({''});
%                 yticklabels({''});
            end       
        end
    end
    
    % save figure and save .mat file
    if p.savePlot
        for k = 1:patNum
            figure(k)
            set(gcf, 'PaperPositionMode', 'auto')
            figName = sprintf('%s%s%d%s%d%s%d%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-Delay_Assembly-',k,'-ArmChoice');
            print(figName,'-dpng','-r300');
        end
    end
    
    if p.writeToFile
        save(fullfile(savedir2,'AssemblyDelay_ArmChoice_SameDelay.mat'), 'AssemblyDelay_ArmChoice_SameDelay');
    end
    clear AssemblyDelay_ArmChoice_SameDelay map_1D
    close all
    fprintf('Finished Arm choice analysis for session %d\n',i);
end

end
