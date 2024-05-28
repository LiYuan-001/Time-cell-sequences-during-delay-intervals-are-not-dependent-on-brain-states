function Fig8TreadmillArmRate_DelayAssembly_Identity_Quant(inFile,AnalyzeSes)
% close all
close all

p.avgRateThres = 0.5;
sessDirs = {'on10','off10','on30','off30'};
for j = 1:length(sessDirs)
    eventTime.(sessDirs{j}) = [];
    eventRate.(sessDirs{j}) = [];
end


assemblyID = [];
assemblyDelayRate = [];
assemblyCenterRate = [];
assemblyRewardRate = [];
timeCell = [];
lateCell = [];
DelayActiveCell = [];
% % Read in input information
sessInfo = SessInfoImport(inFile);

patNum_All = 0;
patNum_All2 = 0;
delayActiveCellNum = 0;
delayInActiveCellNum = 0;

for i = AnalyzeSes(1:end)
    
    sessDirs2 = sessInfo(i).sessDirs;
    % load assembly file
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_Delay.mat');
    load(assemblyFile);
    assemblyRateFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8TreadmillArmRate_SameDelay_DelayAssembly');
    load(assemblyRateFile);
     
    if CellAssembly_Delay.patNum > 0
        
        for j = 1:length(sessDirs)            
            leftSpike = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).delay_L(:,1);
            rightSpike = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).delay_R(:,1);
            leftTime = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).delay_L(:,3);
            rightTime = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).delay_R(:,3);
            spikeNum = leftSpike+rightSpike;
            timeAll = leftTime+rightTime;         
            
            leftSpike = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).stem_L(:,1) + Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).choice_L(:,1);
            rightSpike = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).stem_R(:,1) + Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).choice_R(:,1);
            leftTime = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).stem_L(:,3) + Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).choice_L(:,3);
            rightTime = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).stem_R(:,3) + Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).choice_R(:,3);            
            spikeNum_stemChoice = leftSpike+rightSpike;
            time_stemChoice = leftTime+rightTime;
            
            leftSpike = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).reward_L(:,1);
            rightSpike = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).reward_R(:,1);
            leftTime = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).reward_L(:,3);
            rightTime = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).reward_R(:,3);
            spike_reward = leftSpike+rightSpike;
            time_reward = leftTime+rightTime;
        end
        assemblyDelayRate = [assemblyDelayRate;spikeNum./timeAll];
        assemblyCenterRate = [assemblyCenterRate;spikeNum_stemChoice./time_stemChoice];
        assemblyRewardRate = [assemblyRewardRate;spike_reward./time_reward];
        
        % Get cell identity, early/late cell, delay active cell or non-delay
        % pyr cell
        % load average firing rate file
        SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
        load(SpikeShapeFile);
        % load time field file
        timeFieldFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayField_TrialbyTrial_2Session.mat');
        load(timeFieldFile);
        % load delayFire map
        delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap_2Session.mat');
        load(delayFile);
        % load arm rate file
        armRateFile = fullfile(sessInfo(i).mainDir,'Cell Property', 'Fig8TreadmillArmRate.mat');
        load(armRateFile);
        
        clusterNum = length(SpikeProp.max_AvgRate);
        idx = find(SpikeProp.max_AvgRate>0.1 & SpikeProp.max_AvgRate<5);
        
        rate_Delay = [Fig8TreadmillArmRate.on10.rateDelay(idx,1)';Fig8TreadmillArmRate.off10.rateDelay(idx,1)';...
            Fig8TreadmillArmRate.on30.rateDelay(idx,1)';Fig8TreadmillArmRate.off30.rateDelay(idx,1)'];
        % delay active cell
        delay_onIdx = (sum(rate_Delay>p.avgRateThres,1))>0;
        delayActiveCellNum = delayActiveCellNum + sum(delay_onIdx);

        % other pyr cell
        delay_offIdx = (sum(rate_Delay>p.avgRateThres,1))==0;
        delayInActiveCellNum = delayInActiveCellNum + sum(delay_offIdx);
                
        for j = 1:length(sessDirs)
            
            fieldLabel_Def1.(sessDirs{j}) = [];
            timeMap_Def1.(sessDirs{j}) = [];
            endField.(sessDirs{j}) = [];
            
            fieldLabel_Def1.(sessDirs{j}) = DelayField_TrialbyTrial_2Session.(sessDirs{j}).timeField_1(idx);
            
            sess1 = strcat(sessDirs{j},'_1');
            sess2 = strcat(sessDirs{j},'_2');
            
            for k = idx
                timeMap_Def1.(sessDirs{j}) = [timeMap_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate1_Combined_Smooth{k}];
                if sum(DelayField_TrialbyTrial_2Session.(sessDirs{j}).fieldLabel_1{k}(end-6:end))>0
                    endField.(sessDirs{j}) = [endField.(sessDirs{j}),1];
                else
                    endField.(sessDirs{j}) = [endField.(sessDirs{j}),0];
                end
            end
        end
        
        % early time cell
        on10Idx_Time = (fieldLabel_Def1.on10>0 & delay_onIdx & ~endField.on10);
        off10Idx_Time = (fieldLabel_Def1.off10>0 & delay_onIdx & ~endField.off10);
        on30Idx_Time = (fieldLabel_Def1.on30>0 & delay_onIdx & ~endField.on30);
        off30Idx_Time = (fieldLabel_Def1.off30>0 & delay_onIdx & ~endField.off30);
        Idx_Time = (on10Idx_Time+off10Idx_Time+on30Idx_Time+off30Idx_Time) > 0;
        % late time active cell
        on10Idx_ramp = (fieldLabel_Def1.on10>0 & delay_onIdx & endField.on10);
        off10Idx_ramp = (fieldLabel_Def1.off10>0 & delay_onIdx & endField.off10);
        on30Idx_ramp = (fieldLabel_Def1.on30>0 & delay_onIdx & endField.on30);
        off30Idx_ramp = (fieldLabel_Def1.off30>0 & delay_onIdx & endField.off30);
        Idx_Ramp = (on10Idx_ramp+off10Idx_ramp+on30Idx_ramp+off30Idx_ramp) > 0;        
        Idx_Both = (Idx_Time + Idx_Ramp) > 1;
        
        timeCell = [timeCell,Idx_Time];
        lateCell = [lateCell,Idx_Ramp];
        DelayActiveCell = [DelayActiveCell,delay_onIdx];
        % get assembly cell identity
        % delay_offIdx = 1, delay_onIdx = 2, early time cell = 3, late time
        % cell = 4; early+late = 5;
        %
        patNum = CellAssembly_Delay.patNum ;
        AssmblPtrnCellIDs = CellAssembly_Delay.AssmblPtrnCellIDs;
        AssmblWght = CellAssembly_Delay.AssmblWght;
        AssmblStrength = CellAssembly_Delay.AssmblStrength;
        AssmblStrength(AssmblStrength<0)  = 0;
        event_Time = CellAssembly_Delay.event_Time;
        event_strength = CellAssembly_Delay.event_strength;
        event_Num = CellAssembly_Delay.event_Num;
        binTime = CellAssembly_Delay.binTime;
        
        for k = 1:patNum
            patNum_All = patNum_All + 1;
            cellIDTemp = AssmblPtrnCellIDs{k};
            cellLabel = zeros(1,length(cellIDTemp));
            cellLabel(delay_offIdx(cellIDTemp) == 1) = 1;
            cellLabel(delay_onIdx(cellIDTemp) == 1) = 2;
            cellLabel(Idx_Time(cellIDTemp) == 1) = 3;
            cellLabel(Idx_Ramp(cellIDTemp) == 1) = 4;
            cellLabel(Idx_Both(cellIDTemp) == 1) = 5;
            assemblyID{patNum_All} = cellLabel;
        end
        
        patNum_All2 = patNum_All-patNum;
        sessDirs2 = sessInfo(i).sessDirs;
        if length(sessDirs2) == 8
            sesGroupTemp = 1:4;
        else
            sesGroupTemp = 1:2;
        end
        
        for sesGroup = sesGroupTemp
            
            behaveType = strsplit(sessDirs2{sesGroup},'_');
            behaveType = behaveType{1};
            
            if length(sessDirs2) == 8
                sesTemp = [sesGroup,sesGroup+4];
            elseif length(sessDirs2) == 4
                sesTemp = [sesGroup,sesGroup+2];
            else
                error('session number is wrong')
            end
            for j = sesTemp
                
                % load maps for each cluster
                pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs2{j}, 'PathZone.mat');
                load(pathZoneFile);
                delayFile = fullfile(sessInfo(i).mainDir,sessDirs2{j}, 'Fig8DelayZonePos.mat');
                load(delayFile);
                
                delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
                trialNum = size(delayTstart1,2);
                
                if contains(sessDirs2{j},'10')
                    maxT = 10;
                    delayTend1_2 = delayTstart1+maxT;
                elseif contains(sessDirs2{j},'30')
                    maxT = 30;
                    delayTend1_2 = delayTstart1+maxT;
                else
                    error('Delay time is wrong')
                end
                
                for k = 1:patNum
                    % get spikes
                    tSp = event_Time{k};
                    for m = 1:length(delayTstart1)
                        spikeInd = (tSp>=delayTstart1(m) & tSp<delayTend1_2(m));
                        delayEventTemp = tSp(spikeInd) - delayTstart1(m);
                        eventTime.(behaveType){patNum_All2+k} = delayEventTemp;
                    end
                end
            end
        end
        
    end
end

% label assembly identity based on the cell contents
earlyAssembly = zeros(1,patNum_All);
lateAssembly = zeros(1,patNum_All);
partEarlyAssembly = zeros(1,patNum_All);
partLateAssembly = zeros(1,patNum_All);
otherDelayAssembly = zeros(1,patNum_All);

delayOnAssembly = zeros(1,patNum_All);
delayOffAssembly = zeros(1,patNum_All);
delayMixAssembly = zeros(1,patNum_All);

allCellNum = 0;
assemblyActiveCellNum = 0;
assemblyInActiveCellNum = 0;

for k = 1:patNum_All
    cellLabel = assemblyID{k};
    allCellNum = allCellNum + length(cellLabel);
    assemblyActiveCellNum = assemblyActiveCellNum + sum(cellLabel>1);
    assemblyInActiveCellNum = assemblyInActiveCellNum + sum(cellLabel==1);
    % early time cell
    if sum((cellLabel==5) + (cellLabel==3)) == length(cellLabel)
        earlyAssembly(k) = 1;
    elseif sum((cellLabel==5) + (cellLabel==3)) < length(cellLabel) && sum((cellLabel==5) + (cellLabel==3)) > 0
        partEarlyAssembly(k) = 1;
    end
    % late time cell
    if sum((cellLabel==5) + (cellLabel==4)) == length(cellLabel)
        lateAssembly(k) = 1;
    elseif sum((cellLabel==5) + (cellLabel==4)) < length(cellLabel) && sum((cellLabel==5) + (cellLabel==4)) > 0
        partLateAssembly(k) = 1;
    end
    if sum(cellLabel==2)>0
        otherDelayAssembly(k) = 1;
    end
    
    
    % delayoff 
     if sum((cellLabel==1)) == length(cellLabel)
        delayOffAssembly(k) = 1;
    elseif sum((cellLabel==1)) < length(cellLabel) && sum((cellLabel==1)) > 0
        delayMixAssembly(k) = 1;
     end
    % delayon
     if sum((cellLabel>1)) == length(cellLabel)
        delayOnAssembly(k) = 1;
%     elseif sum((cellLabel>1)) < length(cellLabel) && sum((cellLabel>1)) > 0
%         delayOnAssembly(k) = 1;
    end
end
figure
bar([[assemblyActiveCellNum,assemblyInActiveCellNum]/allCellNum*100,[delayActiveCellNum,delayInActiveCellNum]/sum(delayActiveCellNum+delayInActiveCellNum)*100])
set(gca,'XTick',[1,2,3,4],'XTickLabel',{'Assembly active cell%','inactive','Delay active cell%','inactive'})

figure
bar([sum(delayOnAssembly),sum(delayOffAssembly),sum(delayMixAssembly)]/patNum_All*100)
set(gca,'XTick',[1,2,3],'XTickLabel',{'Delay active','Delay inactive','Mix'})

figure
boxplot(assemblyDelayRate(delayOnAssembly==1),'Position',1)
hold on
boxplot(assemblyDelayRate(delayOffAssembly==1),'Position',2)
boxplot(assemblyDelayRate(delayMixAssembly==1),'Position',3)
xlim([0 4])
set(gca,'XTick',[1,2,3],'XTickLabel',{'Delay active','Delay inactive','Mix'})
title('Delay rate')

figure
boxplot(assemblyCenterRate(delayOnAssembly==1),'Position',1)
hold on
boxplot(assemblyCenterRate(delayOffAssembly==1),'Position',2)
boxplot(assemblyCenterRate(delayMixAssembly==1),'Position',3)
xlim([0 4])
set(gca,'XTick',[1,2,3],'XTickLabel',{'Delay active','Delay inactive','Mix'})
title('Stem + Choice rate')

figure
boxplot(assemblyRewardRate(delayOnAssembly==1),'Position',1)
hold on
boxplot(assemblyRewardRate(delayOffAssembly==1),'Position',2)
boxplot(assemblyRewardRate(delayMixAssembly==1),'Position',3)
xlim([0 4])
set(gca,'XTick',[1,2,3],'XTickLabel',{'Delay active','Delay inactive','Mix'})
title('Reward rate')

% plot wenn of the delay on cell distribution
figure
setListData = {1:length(earlyAssembly); find((earlyAssembly+partEarlyAssembly)>0); find(partLateAssembly)};
setLabels = ["All aseembly"; "Part early"; "Part Late"];
h = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);
title('Cell identity')

figure
setListData = {find(otherDelayAssembly); find((earlyAssembly+partEarlyAssembly)>0); find(partLateAssembly)};
setLabels = ["Other delay active"; "Part early"; "Part Late"];
h = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);
title('Cell identity')

timeAssembly = (earlyAssembly + partEarlyAssembly + partLateAssembly)& ~otherDelayAssembly;
mixType = (earlyAssembly + partEarlyAssembly + partLateAssembly) & otherDelayAssembly;
otherDelayOnly = otherDelayAssembly & ~mixType;

%%
timeAssemblyEventTime_10 = [];
for k = 1:length(timeAssembly)
    if timeAssembly(k) == 1
       timeAssemblyEventTime_10 = [timeAssemblyEventTime_10,eventTime.on10{k},eventTime.off10{k}];
    end
end
mixAssemblyEventTime_10 = [];
for k = 1:length(timeAssembly)
    if mixType(k) == 1
       mixAssemblyEventTime_10 = [mixAssemblyEventTime_10,eventTime.on10{k},eventTime.off10{k}];
    end
end
otherAssemblyEventTime_10 = [];
for k = 1:length(timeAssembly)
    if otherDelayOnly(k) == 1
       otherAssemblyEventTime_10 = [otherAssemblyEventTime_10,eventTime.on10{k},eventTime.off10{k}];
    end
end

figure
lim_10 = 0:0.5:10;
count_timeAssembly_10 = histcounts(timeAssemblyEventTime_10,lim_10);
count_mixAssembly_10 = histcounts(mixAssemblyEventTime_10,lim_10);
count_otherAssembly_10 = histcounts(otherAssemblyEventTime_10,lim_10);

plot(count_timeAssembly_10)
hold on
plot(count_mixAssembly_10)
plot(count_otherAssembly_10)

legend({'Stable cell assembly','Mix assembly','Non-stable cell assembly'},'Location','EastOutside')

timeAssemblyEventTime_30 = [];
for k = 1:length(timeAssembly)
    if timeAssembly(k) == 1
       timeAssemblyEventTime_30 = [timeAssemblyEventTime_30,eventTime.on30{k},eventTime.off30{k}];
    end
end
mixAssemblyEventTime_30 = [];
for k = 1:length(timeAssembly)
    if mixType(k) == 1
       mixAssemblyEventTime_30 = [mixAssemblyEventTime_30,eventTime.on30{k},eventTime.off30{k}];
    end
end
otherAssemblyEventTime_30 = [];
for k = 1:length(timeAssembly)
    if otherDelayOnly(k) == 1
       otherAssemblyEventTime_30 = [otherAssemblyEventTime_30,eventTime.on30{k},eventTime.off30{k}];
    end
end

figure
lim_30 = 0:0.5:30;
count_timeAssembly_30 = histcounts(timeAssemblyEventTime_30,lim_30);
count_mixAssembly_30 = histcounts(mixAssemblyEventTime_30,lim_30);
count_otherAssembly_30 = histcounts(otherAssemblyEventTime_30,lim_30);
plot(count_timeAssembly_30)
hold on
plot(count_mixAssembly_30)
plot(count_otherAssembly_30)

figure
boxplot(assemblyDelayRate(timeAssembly==1),'Position',1)
hold on
boxplot(assemblyDelayRate(mixType==1),'Position',2)
boxplot(assemblyDelayRate(otherDelayOnly==1),'Position',3)
xlim([0 4])
set(gca,'XTick',[1,2,3],'XTickLabel',{'Stable cell','mix','Other delay'})
title('Delay rate')

figure
boxplot(assemblyCenterRate(timeAssembly==1),'Position',1)
hold on
boxplot(assemblyCenterRate(otherDelayAssembly==1),'Position',2)
boxplot(assemblyDelayRate(otherDelayAssembly==1),'Position',3)
xlim([0 4])
set(gca,'XTick',[1,2,3],'XTickLabel',{'Stable cell','mix','Other delay'})
title('Stem + Choice rate')

figure
boxplot(assemblyRewardRate(timeAssembly==1),'Position',1)
hold on
boxplot(assemblyRewardRate(otherDelayAssembly==1),'Position',2)
boxplot(assemblyRewardRate(otherDelayAssembly==1),'Position',3)
xlim([0 4])
set(gca,'XTick',[1,2,3],'XTickLabel',{'Stable cell','mix','Other delay'})
title('Reward rate')


figure
setListData = {1:length(DelayActiveCell); find(timeCell); find(lateCell)};
setLabels = ["Delay active"; "Time cell"; "Constent on"];
h = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);
title('Cell identity')

figure
bar([sum(earlyAssembly),sum(partEarlyAssembly),sum(lateAssembly),sum(partLateAssembly)]/patNum_All*100)
set(gca,'XTick',[1,2,3,4],'XTickLabel',{'Time cell','Part time cell','Late cell','Part late cell'})

timeMixAssembly = (earlyAssembly+partEarlyAssembly) & partLateAssembly;
earlyOnly = (earlyAssembly+partEarlyAssembly) & ~timeMixAssembly;
lateOnly = partLateAssembly & ~timeMixAssembly;

timeAssembly = earlyAssembly + partEarlyAssembly + partLateAssembly;
nontimeAssembly = [1:length(earlyAssembly)] & ~timeAssembly;

figure
boxplot(assemblyDelayRate(timeAssembly==1),'Position',1)
hold on
boxplot(assemblyDelayRate(nontimeAssembly==1),'Position',2)
xlim([0 4])
set(gca,'XTick',[1,2],'XTickLabel',{'Stable delay','Other delay'})
title('Delay rate')

figure
boxplot(assemblyCenterRate(timeAssembly==1),'Position',1)
hold on
boxplot(assemblyCenterRate(nontimeAssembly==1),'Position',2)
xlim([0 4])
set(gca,'XTick',[1,2],'XTickLabel',{'Stable delay','Other delay'})
title('Stem + Choice rate')

figure
boxplot(assemblyRewardRate(timeAssembly==1),'Position',1)
hold on
boxplot(assemblyRewardRate(nontimeAssembly==1),'Position',2)
xlim([0 4])
set(gca,'XTick',[1,2],'XTickLabel',{'Stable delay','Other delay'})
title('Reward rate')


figure
boxplot(assemblyDelayRate(earlyAssembly==1),'Position',1)
hold on
boxplot(assemblyDelayRate(partEarlyAssembly==1),'Position',2)
boxplot(assemblyDelayRate(partLateAssembly==1),'Position',4)
xlim([0 5])
set(gca,'XTick',[1,2,3,4],'XTickLabel',{'Time cell','Part time cell','Late cell','Part late cell'})

end