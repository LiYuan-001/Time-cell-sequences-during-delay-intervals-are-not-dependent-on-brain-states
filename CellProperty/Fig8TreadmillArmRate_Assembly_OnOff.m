function Fig8TreadmillArmRate_Assembly_OnOff(inFile,AnalyzeSes)

%  Li Yuan, UCSD, 11-May-2022
%  Calculate firing rate of each arm
%  for each assembly
% -------------------------------------------------------------------------
% set parameters for analysis

p.savePlot = 0;
p.writeToFile = 1;

% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    close all
    % selectivity for the whole map
    ArmRate_SameDelay_Assembly_OnOff.session = sessInfo(i).mainDir;
    ArmRate_SameDelay_Assembly_OnOff.rat = sessInfo(i).animal;
    ArmRate_SameDelay_Assembly_OnOff.day = sessInfo(i).day;
    
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
    
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load assembly file
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_DelayOnOff_LR-25ms.mat');
    load(assemblyFile);
    
    % get time stamps for the assembly
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;    
    tList = SpikeProp.tList(rateLabel);
    
    patNum_on = CellAssembly_DelayLR.DelayOn.patNum; 
    AssmblPtrnCellIDs_on = CellAssembly_DelayLR.DelayOn.AssmblPtrnCellIDs;
    AssmblWght_on = CellAssembly_DelayLR.DelayOn.AssmblWght;
    AssmblStrength_on = CellAssembly_DelayLR.DelayOn.AssmblStrength;
    event_Time_on = CellAssembly_DelayLR.DelayOn.event_Time;
    event_strength_on = CellAssembly_DelayLR.DelayOn.event_strength;
    event_Num_on = CellAssembly_DelayLR.DelayOn.event_Num;
    
    patNum_off = CellAssembly_DelayLR.DelayOff.patNum; 
    AssmblPtrnCellIDs_off = CellAssembly_DelayLR.DelayOff.AssmblPtrnCellIDs;
    AssmblWght_onff = CellAssembly_DelayLR.DelayOff.AssmblWght;
    AssmblStrength_off = CellAssembly_DelayLR.DelayOff.AssmblStrength;
    event_Time_off = CellAssembly_DelayLR.DelayOff.event_Time;
    event_strength_off = CellAssembly_DelayLR.DelayOff.event_strength;
    event_Num_off = CellAssembly_DelayLR.DelayOff.event_Num;
 
    
    if prod(~isempty(char(sessDirs{1})))
        
        if length(sessDirs) == 8
            sesGroupTemp = 1:4;
        else
            sesGroupTemp = 1:2;
        end
        
        for sesGroup = sesGroupTemp
            sesNumSum = 0;
            posLabel = [];
            turnLabel = [];
            time_All = [];
            
            behaveType = strsplit(sessDirs{sesGroup},'_');
            behaveType = behaveType{1};
            
            if length(sessDirs) == 8
                sesTemp = [sesGroup,sesGroup+4];
            elseif length(sessDirs) == 4
                sesTemp = [sesGroup,sesGroup+2];
            else
                error('session number is wrong')
            end
            for j = sesTemp
                
                % load maps for each cluster
                pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'PathZone.mat');
                load(pathZoneFile);
                delayFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
                load(delayFile);
                
                if exist(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'), 'file')==2
                    map_1D{j} = load(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'));
                    
                    for m = 1:sum(map_1D{j}.ratesByECLR.valid)
                        if map_1D{j}.ratesByECLR.ECLR(m) == 1
                            posLabel(sesNumSum+m,1) = 2;
                            turnLabel(sesNumSum+m) = 1;
                        elseif map_1D{j}.ratesByECLR.ECLR(m) == 2
                            posLabel(sesNumSum+m,1) = 1;
                            turnLabel(sesNumSum+m) = 0;
                        elseif map_1D{j}.ratesByECLR.ECLR(m) == 3
                            posLabel(sesNumSum+m,1) = 1;
                            turnLabel(sesNumSum+m) = 1;
                        elseif map_1D{j}.ratesByECLR.ECLR(m) == 4
                            posLabel(sesNumSum+m,1) = 2;
                            turnLabel(sesNumSum+m) = 0;
                        else
                            error('Turning label ERROR');
                        end
                        % return
                        timeTemp = PathZone.posEndT.Return(m)-PathZone.posStartT.Return(m);
                        time_All(sesNumSum+m,1) = timeTemp;
                        % base
                        timeTemp = PathZone.posEndT.Base(m)-PathZone.posStartT.Base(m);
                        time_All(sesNumSum+m,2) = timeTemp;
                        % delay
                        timeTemp = Fig8DelayZonePos.delayPos2.endT(m)-Fig8DelayZonePos.delayPos2.startT(m);
                        time_All(sesNumSum+m,3) = timeTemp;
                        % stem
                        timeTemp = PathZone.posEndT.Center(m)-(Fig8DelayZonePos.delayPos2.endT(m)+1/30);
                        time_All(sesNumSum+m,4) = timeTemp;
                        % choice
                        timeTemp = PathZone.posEndT.Choice(m)-PathZone.posStartT.Choice(m);
                        time_All(sesNumSum+m,5) = timeTemp;
                        % reward
                        timeTemp = PathZone.posEndT.Reward(m)-PathZone.posStartT.Reward(m);
                        time_All(sesNumSum+m,6) = timeTemp;
                        
                    end
                    sesNumSum = sesNumSum + sum(map_1D{j}.ratesByECLR.valid);
                end
                
            end
            
            
            % find Left / Right turn trials
            leftIdx = find(posLabel(:,1)==1); % idx on all trials
            rightIdx = find(posLabel(:,1)==2); % idx on all trials
            
            LcorrectIdx = find(turnLabel' == 1 & posLabel(:,1)==2); % idx on all trials
            LerrorIdx = find(turnLabel' ==0 & posLabel(:,1)==1); % idx on all trials
            RcorrectIdx = find(turnLabel' == 1 & posLabel(:,1)==1); % idx on all trials
            RerrorIdx = find(turnLabel' == 0 & posLabel(:,1)==2); % idx on all trials
            
            
            for k = 1:patNum_on
                sessionTemp = 1;
                % initiate variables
                % return, delay, stem
                spike_All = nan(sesNumSum,6);
                strength_All = nan(sesNumSum,6);
                
                if length(sessDirs) == 8
                    sesTemp = [sesGroup,sesGroup+4];
                elseif length(sessDirs) == 4
                    sesTemp = [sesGroup,sesGroup+2];
                else
                    error('session number is wrong')
                end
                for j = sesTemp
                    % load path
                    pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'PathZone.mat');
                    load(pathZoneFile);
                    delayFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
                    load(delayFile);
                    
                    % get spikes
                    tSp = event_Time_on{k};
                    strength = event_strength_on{k};
                    
                    for m = 1:sum(map_1D{j}.ratesByECLR.valid)
                        % get spike count
                        % return
                        spikeInd = (tSp>=PathZone.posStartT.Return(m) & tSp<PathZone.posEndT.Return(m));
                        spike_All(sessionTemp,1) = sum(spikeInd);
                        strength_All(sessionTemp,1) = sum(strength(spikeInd));
                        % base
                        spikeInd = (tSp>=PathZone.posStartT.Base(m) & tSp<PathZone.posEndT.Base(m));
                        spike_All(sessionTemp,2) = sum(spikeInd);
                        strength_All(sessionTemp,2) = sum(strength(spikeInd));
                        % delay
                        spikeInd = (tSp>=Fig8DelayZonePos.delayPos2.startT(m) & tSp<Fig8DelayZonePos.delayPos2.endT(m));
                        spike_All(sessionTemp,3) = sum(spikeInd);
                        strength_All(sessionTemp,3) = sum(strength(spikeInd));
                        % stem
                        spikeInd = (tSp>=(Fig8DelayZonePos.delayPos2.endT(m)+1/30) & tSp<PathZone.posEndT.Center(m));
                        spike_All(sessionTemp,4) = sum(spikeInd);
                        strength_All(sessionTemp,4) = sum(strength(spikeInd));
                        % choice
                        spikeInd = (tSp>=PathZone.posStartT.Choice(m) & tSp<PathZone.posEndT.Choice(m));
                        spike_All(sessionTemp,5) = sum(spikeInd);
                        strength_All(sessionTemp,5) = sum(strength(spikeInd));
                        % reward
                        spikeInd = (tSp>=PathZone.posStartT.Reward(m) & tSp<PathZone.posEndT.Reward(m));
                        spike_All(sessionTemp,6) = sum(spikeInd);
                        strength_All(sessionTemp,6) = sum(strength(spikeInd));
                        
                        sessionTemp = sessionTemp+1;
                    end
                end
                
                % firing rate in different arms
                rateReturn_L = nansum(nansum(spike_All(leftIdx,1)))./sum(sum(time_All(leftIdx,1)));
                rateReturn_R= nansum(nansum(spike_All(rightIdx,1)))./sum(sum(time_All(rightIdx,1)));
                
                rateReturn = nansum(nansum(spike_All(:,1)))./sum(sum(time_All(:,1)));
                eventStrengthReturn = nansum(nansum(strength_All(:,1)))./sum(sum(spike_All(:,1)));
                % base
                rateBase = nansum(nansum(spike_All(:,2)))./sum(sum(time_All(:,2)));
                eventStrengthBase = nansum(nansum(strength_All(:,2)))./sum(sum(spike_All(:,2)));
                % delay zone
                rateDelay = nansum(nansum(spike_All(:,3)))./sum(sum(time_All(:,3)));
                eventStrengthDelay = nansum(nansum(strength_All(:,3)))./sum(sum(spike_All(:,3)));
                % stem zone
                rateStem = nansum(nansum(spike_All(:,4)))./sum(sum(time_All(:,4)));
                eventStrengthStem = nansum(nansum(strength_All(:,4)))./sum(sum(spike_All(:,4)));
                % choice zone
                rateChoice = nansum(nansum(spike_All(:,5)))./sum(sum(time_All(:,5)));
                eventStrengthChoice = nansum(nansum(strength_All(:,5)))./sum(sum(spike_All(:,5)));
                % reward zone
                rateReward = nansum(nansum(spike_All(:,6)))./sum(sum(time_All(:,6)));
                eventStrengthReward = nansum(nansum(strength_All(:,6)))./sum(sum(spike_All(:,6)));

                ArmRate_SameDelay_Assembly_OnOff.on.(behaveType).rateReturn(k) = rateReturn;
                ArmRate_SameDelay_Assembly_OnOff.on.(behaveType).rateReturn_L(k) = rateReturn_L;
                ArmRate_SameDelay_Assembly_OnOff.on.(behaveType).rateReturn_R(k) = rateReturn_R;
                ArmRate_SameDelay_Assembly_OnOff.on.(behaveType).rateBase(k) = rateBase;
                ArmRate_SameDelay_Assembly_OnOff.on.(behaveType).rateDelay(k) = rateDelay;
                ArmRate_SameDelay_Assembly_OnOff.on.(behaveType).rateStem(k) = rateStem;
                ArmRate_SameDelay_Assembly_OnOff.on.(behaveType).rateChoice(k) = rateChoice;
                ArmRate_SameDelay_Assembly_OnOff.on.(behaveType).rateReward(k) = rateReward;
                
                ArmRate_SameDelay_Assembly_OnOff.on.(behaveType).eventStrengthReturn(k) = eventStrengthReturn;
                ArmRate_SameDelay_Assembly_OnOff.on.(behaveType).eventStrengthBase(k) = eventStrengthBase;
                ArmRate_SameDelay_Assembly_OnOff.on.(behaveType).eventStrengthDelay(k) = eventStrengthDelay;
                ArmRate_SameDelay_Assembly_OnOff.on.(behaveType).eventStrengthStem(k) = eventStrengthStem;
                ArmRate_SameDelay_Assembly_OnOff.on.(behaveType).eventStrengthChoice(k) = eventStrengthChoice;
                ArmRate_SameDelay_Assembly_OnOff.on.(behaveType).eventStrengthReward(k) = eventStrengthReward;
                
            end
            
            for k = 1:patNum_off
                sessionTemp = 1;
                % initiate variables
                % return, delay, stem
                spike_All = nan(sesNumSum,6);
                strength_All = nan(sesNumSum,6);
                
                if length(sessDirs) == 8
                    sesTemp = [sesGroup,sesGroup+4];
                elseif length(sessDirs) == 4
                    sesTemp = [sesGroup,sesGroup+2];
                else
                    error('session number is wrong')
                end
                for j = sesTemp
                    % load path
                    pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'PathZone.mat');
                    load(pathZoneFile);
                    delayFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
                    load(delayFile);
                    
                    % get spikes
                    tSp = event_Time_off{k};
                    strength = event_strength_off{k};
                    
                    for m = 1:sum(map_1D{j}.ratesByECLR.valid)
                        % get spike count
                        % return
                        spikeInd = (tSp>=PathZone.posStartT.Return(m) & tSp<PathZone.posEndT.Return(m));
                        spike_All(sessionTemp,1) = sum(spikeInd);
                        strength_All(sessionTemp,1) = sum(strength(spikeInd));
                        % base
                        spikeInd = (tSp>=PathZone.posStartT.Base(m) & tSp<PathZone.posEndT.Base(m));
                        spike_All(sessionTemp,2) = sum(spikeInd);
                        strength_All(sessionTemp,2) = sum(strength(spikeInd));
                        % delay
                        spikeInd = (tSp>=Fig8DelayZonePos.delayPos2.startT(m) & tSp<Fig8DelayZonePos.delayPos2.endT(m));
                        spike_All(sessionTemp,3) = sum(spikeInd);
                        strength_All(sessionTemp,3) = sum(strength(spikeInd));
                        % stem
                        spikeInd = (tSp>=(Fig8DelayZonePos.delayPos2.endT(m)+1/30) & tSp<PathZone.posEndT.Center(m));
                        spike_All(sessionTemp,4) = sum(spikeInd);
                        strength_All(sessionTemp,4) = sum(strength(spikeInd));
                        % choice
                        spikeInd = (tSp>=PathZone.posStartT.Choice(m) & tSp<PathZone.posEndT.Choice(m));
                        spike_All(sessionTemp,5) = sum(spikeInd);
                        strength_All(sessionTemp,5) = sum(strength(spikeInd));
                        % reward
                        spikeInd = (tSp>=PathZone.posStartT.Reward(m) & tSp<PathZone.posEndT.Reward(m));
                        spike_All(sessionTemp,6) = sum(spikeInd);
                        strength_All(sessionTemp,6) = sum(strength(spikeInd));
                        
                        sessionTemp = sessionTemp+1;
                    end
                end
                
                % firing rate in different arms
                rateReturn_L = nansum(nansum(spike_All(leftIdx,1)))./sum(sum(time_All(leftIdx,1)));
                rateReturn_R= nansum(nansum(spike_All(rightIdx,1)))./sum(sum(time_All(rightIdx,1)));
                
                rateReturn = nansum(nansum(spike_All(:,1)))./sum(sum(time_All(:,1)));
                eventStrengthReturn = nansum(nansum(strength_All(:,1)))./sum(sum(spike_All(:,1)));
                % base
                rateBase = nansum(nansum(spike_All(:,2)))./sum(sum(time_All(:,2)));
                eventStrengthBase = nansum(nansum(strength_All(:,2)))./sum(sum(spike_All(:,2)));
                % delay zone
                rateDelay = nansum(nansum(spike_All(:,3)))./sum(sum(time_All(:,3)));
                eventStrengthDelay = nansum(nansum(strength_All(:,3)))./sum(sum(spike_All(:,3)));
                % stem zone
                rateStem = nansum(nansum(spike_All(:,4)))./sum(sum(time_All(:,4)));
                eventStrengthStem = nansum(nansum(strength_All(:,4)))./sum(sum(spike_All(:,4)));
                % choice zone
                rateChoice = nansum(nansum(spike_All(:,5)))./sum(sum(time_All(:,5)));
                eventStrengthChoice = nansum(nansum(strength_All(:,5)))./sum(sum(spike_All(:,5)));
                % reward zone
                rateReward = nansum(nansum(spike_All(:,6)))./sum(sum(time_All(:,6)));
                eventStrengthReward = nansum(nansum(strength_All(:,4)))./sum(sum(spike_All(:,4)));

                ArmRate_SameDelay_Assembly_OnOff.off.(behaveType).rateReturn(k) = rateReturn;
                ArmRate_SameDelay_Assembly_OnOff.off.(behaveType).rateReturn_L(k) = rateReturn_L;
                ArmRate_SameDelay_Assembly_OnOff.off.(behaveType).rateReturn_R(k) = rateReturn_R;
                ArmRate_SameDelay_Assembly_OnOff.off.(behaveType).rateBase(k) = rateBase;
                ArmRate_SameDelay_Assembly_OnOff.off.(behaveType).rateDelay(k) = rateDelay;
                ArmRate_SameDelay_Assembly_OnOff.off.(behaveType).rateStem(k) = rateStem;
                ArmRate_SameDelay_Assembly_OnOff.off.(behaveType).rateChoice(k) = rateChoice;
                ArmRate_SameDelay_Assembly_OnOff.off.(behaveType).rateReward(k) = rateReward;
                
                ArmRate_SameDelay_Assembly_OnOff.off.(behaveType).eventStrengthReturn(k) = eventStrengthReturn;
                ArmRate_SameDelay_Assembly_OnOff.off.(behaveType).eventStrengthBase(k) = eventStrengthBase;
                ArmRate_SameDelay_Assembly_OnOff.off.(behaveType).eventStrengthDelay(k) = eventStrengthDelay;
                ArmRate_SameDelay_Assembly_OnOff.off.(behaveType).eventStrengthStem(k) = eventStrengthStem;
                ArmRate_SameDelay_Assembly_OnOff.off.(behaveType).eventStrengthChoice(k) = eventStrengthChoice;
                ArmRate_SameDelay_Assembly_OnOff.off.(behaveType).eventStrengthReward(k) = eventStrengthReward;
                
            end
        end
        
        % sleep sessions
        sleepDirs = sessInfo(i).sleepDirs;
        posFile = fullfile(sessInfo(i).mainDir,'processedData','indataS.mat');
        load(posFile);
        for j =  1:length(sleepDirs)
            pos = indata(j);
            startT = pos.t(1);
            endT = pos.t(end);
            
            for k = 1:patNum_on
                % get spikes
                tSp = event_Time_on{k};
                strength = event_strength_on{k};
                spikeInd = (tSp>=startT & tSp<endT);
                rateSleep = sum(spikeInd)/(endT-startT);
                eventStrengthSleep = sum(strength(spikeInd))/length(spikeInd);
                ArmRate_SameDelay_Assembly_OnOff.on.(sleepDirs{j}).rateSleep(k) = rateSleep;
                ArmRate_SameDelay_Assembly_OnOff.on.(sleepDirs{j}).eventStrengthSleep(k) = eventStrengthSleep;
            end
            
            for k = 1:patNum_off
                % get spikes
                tSp = event_Time_off{k};
                strength = event_strength_off{k};
                spikeInd = (tSp>=startT & tSp<endT);
                rateSleep = sum(spikeInd)/(endT-startT);
                eventStrengthSleep = sum(strength(spikeInd))/length(spikeInd);
                ArmRate_SameDelay_Assembly_OnOff.off.(sleepDirs{j}).rateSleep(k) = rateSleep;
                ArmRate_SameDelay_Assembly_OnOff.off.(sleepDirs{j}).eventStrengthSleep(k) = eventStrengthSleep;
            end
        end
    
        
        if p.writeToFile
            save(fullfile(savedir2,'ArmRate_SameDelay_Assembly_OnOff-25ms.mat'), 'ArmRate_SameDelay_Assembly_OnOff');
        end
        clear ArmRate_SameDelay_Assembly_OnOff map_1D
        close all
        fprintf('Finished Arm analysis for session %d\n',i);
    end

end
