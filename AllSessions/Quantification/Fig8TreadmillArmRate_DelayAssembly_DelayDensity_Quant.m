function Fig8TreadmillArmRate_DelayAssembly_DelayDensity_Quant(inFile,AnalyzeSes)
% close all
close all

sessDirs = {'on10','off10','on30','off30'};

for j = 1:length(sessDirs)
    eventTime.(sessDirs{j}) = [];
    eventRate.(sessDirs{j}) = [];
end

% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    sessDirs2 = sessInfo(i).sessDirs;
    % load assembly file
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_Delay.mat');
    load(assemblyFile);
    
    patNum = CellAssembly_Delay.patNum ;
    AssmblPtrnCellIDs = CellAssembly_Delay.AssmblPtrnCellIDs;
    AssmblWght = CellAssembly_Delay.AssmblWght;
    AssmblStrength = CellAssembly_Delay.AssmblStrength;
    AssmblStrength(AssmblStrength<0)  = 0;
    event_Time = CellAssembly_Delay.event_Time;
    event_strength = CellAssembly_Delay.event_strength;
    event_Num = CellAssembly_Delay.event_Num;
    binTime = CellAssembly_Delay.binTime;
    
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
                    eventTime.(behaveType) = [eventTime.(behaveType),delayEventTemp];
                end
            end
        end
    end
end

figure(1)
lim_10 = 0:0.5:10;
lim_30 = 0:0.5:30;

count_on10 = histcounts(eventTime.on10,lim_10);
count_off10 = histcounts(eventTime.off10,lim_10);
count_on30 = histcounts(eventTime.on30,lim_30);
count_off30 = histcounts(eventTime.off30,lim_30);

plot_max = max([count_on10,count_off10,count_on30,count_off30]);
colormap(jet)

returnDiff = (ArmStrength.on30.return - ArmStrength.off30.return)./(ArmStrength.on30.return + ArmStrength.off30.return);
delayDiff = (ArmStrength.on30.delay - ArmStrength.off30.delay)./(ArmStrength.on30.delay + ArmStrength.off30.delay);
stemDiff = (ArmStrength.on30.stem - ArmStrength.off30.stem)./(ArmStrength.on30.stem + ArmStrength.off30.stem);
choiceDiff = (ArmStrength.on30.choice - ArmStrength.off30.choice)./(ArmStrength.on30.choice + ArmStrength.off30.choice);


boxplot(ArmStrength.on30.rateDelay','positions',1,'Colors','m','Symbol','');
hold on
boxplot(AssemblyRate.off10.rateDelay','positions',2,'Colors','m','Symbol','');
boxplot(AssemblyRate.on30.rateDelay','positions',3,'Colors','m','Symbol','');
boxplot(AssemblyRate.off30.rateDelay','positions',4,'Colors','m','Symbol','');

boxplot(AssemblyRate.on10.rateStem','positions',6,'Colors','k','Symbol','');
boxplot(AssemblyRate.off10.rateStem','positions',7,'Colors','k','Symbol','');
boxplot(AssemblyRate.on30.rateStem','positions',8,'Colors','k','Symbol','');
boxplot(AssemblyRate.off30.rateStem','positions',9,'Colors','k','Symbol','');

boxplot(AssemblyRate.on10.rateChoice','positions',11,'Colors','r','Symbol','');
boxplot(AssemblyRate.off10.rateChoice','positions',12,'Colors','r','Symbol','');
boxplot(AssemblyRate.on30.rateChoice','positions',13,'Colors','r','Symbol','');
boxplot(AssemblyRate.off30.rateChoice','positions',14,'Colors','r','Symbol','');

boxplot(AssemblyRate.on10.rateReward','positions',16,'Colors','g','Symbol','');
boxplot(AssemblyRate.off10.rateReward','positions',17,'Colors','g','Symbol','');
boxplot(AssemblyRate.on30.rateReward','positions',18,'Colors','g','Symbol','');
boxplot(AssemblyRate.off30.rateReward','positions',19,'Colors','g','Symbol','');

boxplot(AssemblyRate.on10.rateReturn','positions',21,'Colors','b','Symbol','');
boxplot(AssemblyRate.off10.rateReturn','positions',22,'Colors','b','Symbol','');
boxplot(AssemblyRate.on30.rateReturn','positions',23,'Colors','b','Symbol','');
boxplot(AssemblyRate.off30.rateReturn','positions',24,'Colors','b','Symbol','');

boxplot(AssemblyRate.sleep1','positions',26,'Colors','k','Symbol','');
boxplot(AssemblyRate.sleep2','positions',27,'Colors','k','Symbol','');

xlim([0 28])
set(gca,'XTick',[1,6,11,16,21,26,27],'XTickLabel',{'Delay','Stem','Choice','Reward','Return','Sleep1','Sleep2'});
title('Assembly event average rate')

figure(2)
boxplot(ArmStrength.on10.eventStrengthDelay','positions',1,'Colors','m','Symbol','');
hold on
boxplot(ArmStrength.off10.eventStrengthDelay','positions',2,'Colors','m','Symbol','');
boxplot(ArmStrength.on30.eventStrengthDelay','positions',3,'Colors','m','Symbol','');
boxplot(ArmStrength.off30.eventStrengthDelay','positions',4,'Colors','m','Symbol','');

boxplot(ArmStrength.on10.eventStrengthStem','positions',6,'Colors','k','Symbol','');
boxplot(ArmStrength.off10.eventStrengthStem','positions',7,'Colors','k','Symbol','');
boxplot(ArmStrength.on30.eventStrengthStem','positions',8,'Colors','k','Symbol','');
boxplot(ArmStrength.off30.eventStrengthStem','positions',9,'Colors','k','Symbol','');

boxplot(ArmStrength.on10.eventStrengthChoice','positions',11,'Colors','r','Symbol','');
boxplot(ArmStrength.off10.eventStrengthChoice','positions',12,'Colors','r','Symbol','');
boxplot(ArmStrength.on30.eventStrengthChoice','positions',13,'Colors','r','Symbol','');
boxplot(ArmStrength.off30.eventStrengthChoice','positions',14,'Colors','r','Symbol','');

boxplot(ArmStrength.on10.eventStrengthReward','positions',16,'Colors','g','Symbol','');
boxplot(ArmStrength.off10.eventStrengthReward','positions',17,'Colors','g','Symbol','');
boxplot(ArmStrength.on30.eventStrengthReward','positions',18,'Colors','g','Symbol','');
boxplot(ArmStrength.off30.eventStrengthReward','positions',19,'Colors','g','Symbol','');

boxplot(ArmStrength.on10.eventStrengthReturn','positions',21,'Colors','b','Symbol','');
boxplot(ArmStrength.off10.eventStrengthReturn','positions',22,'Colors','b','Symbol','');
boxplot(ArmStrength.on30.eventStrengthReturn','positions',23,'Colors','b','Symbol','');
boxplot(ArmStrength.off30.eventStrengthReturn','positions',24,'Colors','b','Symbol','');

boxplot(ArmStrength.sleep1','positions',26,'Colors','k','Symbol','');
boxplot(ArmStrength.sleep2','positions',27,'Colors','k','Symbol','');

xlim([0 28])
set(gca,'XTick',[1,6,11,16,21,26,27],'XTickLabel',{'Delay','Stem','Choice','Reward','Return','Sleep1','Sleep2'});
title('Assembly event average strength')

%    title({'Delay','Stem-Delay/Stem+Delay'})
%    set(gca,'XTick',[1,2,5,6],'XTickLabel',{'All trials Ctr','Lesion','Correct Ctr','Lesion'});
%    [H_Stem2Delayall,P_Stem2Delayall]=kstest2(CtrThres.Stem2Delayall',LesionThres.Stem2Delayall');
%    [H_Stem2Delayc,P_Stem2Delayc]=kstest2(CtrThres.Stem2Delayc',LesionThres.Stem2Delayc');

figure(3)
Violin(AssemblyRate.on10.rateDelay',1);
hold on
Violin(AssemblyRate.off10.rateDelay',2);
Violin(AssemblyRate.on30.rateDelay',3);
Violin(AssemblyRate.off30.rateDelay',4);

Violin(AssemblyRate.on10.rateStem',6);
Violin(AssemblyRate.off10.rateStem',7);
Violin(AssemblyRate.on30.rateStem',8);
Violin(AssemblyRate.off30.rateStem',9);

Violin(AssemblyRate.on10.rateChoice',11);
Violin(AssemblyRate.off10.rateChoice',12);
Violin(AssemblyRate.on30.rateChoice',13);
Violin(AssemblyRate.off30.rateChoice',14);

Violin(AssemblyRate.on10.rateReward',16);
Violin(AssemblyRate.off10.rateReward',17);
Violin(AssemblyRate.on30.rateReward',18);
Violin(AssemblyRate.off30.rateReward',19);

Violin(AssemblyRate.on10.rateReturn',21);
Violin(AssemblyRate.off10.rateReturn',22);
Violin(AssemblyRate.on30.rateReturn',23);
Violin(AssemblyRate.off30.rateReturn',24);

Violin(AssemblyRate.sleep1',26);
Violin(AssemblyRate.sleep2',27);

xlim([0 28])
set(gca,'XTick',[1,6,11,16,21,26,27],'XTickLabel',{'Delay','Stem','Choice','Reward','Return','Sleep1','Sleep2'});
title('Assembly event average rate')


figure(2)
boxplot(AssemblyRate.on10.rateStem','positions',1.6,'Colors','m','Symbol','');
hold on
boxplot(AssemblyRate.off10.rateStem','positions',4.6,'Colors','m','Symbol','');
boxplot(AssemblyRate.on30.rateStem','positions',7.6,'Colors','m','Symbol','');
boxplot(AssemblyRate.off30.rateStem','positions',10.6,'Colors','m','Symbol','');

% boxplot(ArmRate.on10.rateStem','positions',2,'Colors','r','Symbol','');
hold on
% boxplot(ArmRate.off10.rateStem','positions',5,'Colors','k','Symbol','');
% boxplot(ArmRate.on30.rateStem','positions',8,'Colors','r','Symbol','');
% boxplot(ArmRate.off10.rateStem','positions',11,'Colors','k','Symbol','');
xlim([0 12])
set(gca,'XTick',[1,2,4,5,7,8,10,11],'XTickLabel',xaxisName);

figure(3)
boxplot(AssemblyRate.on10.rateReturn','positions',1.6,'Colors','m','Symbol','');
hold on
boxplot(AssemblyRate.off10.rateReturn','positions',4.6,'Colors','m','Symbol','');
boxplot(AssemblyRate.on30.rateReturn','positions',7.6,'Colors','m','Symbol','');
boxplot(AssemblyRate.off30.rateReturn','positions',10.6,'Colors','m','Symbol','');

% boxplot(ArmRate.on10.rateReturn','positions',2,'Colors','r','Symbol','');
hold on
% boxplot(ArmRate.off10.rateReturn','positions',5,'Colors','k','Symbol','');
% boxplot(ArmRate.on30.rateReturn','positions',8,'Colors','r','Symbol','');
% boxplot(ArmRate.off10.rateReturn','positions',11,'Colors','k','Symbol','');
xlim([0 12])
set(gca,'XTick',[1,2,4,5,7,8,10,11],'XTickLabel',xaxisName);

figure(4)
boxplot(AssemblyRate.rate','positions',1.6,'Colors','m','Symbol','');
% % boxplot(ArmRate.rate','positions',2,'Colors','r','Symbol','');
% hold on
% % boxplot(ArmRate.rate','positions',3,'Colors','k','Symbol','');
xlim([0 4])
set(gca,'XTick',[1,2,3],'XTickLabel',{'On&Off','On','Off'});
end