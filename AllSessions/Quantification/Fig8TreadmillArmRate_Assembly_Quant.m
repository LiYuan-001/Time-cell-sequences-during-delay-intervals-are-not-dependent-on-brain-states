function Fig8TreadmillArmRate_Assembly_Quant(inFile,AnalyzeSes)
% close all
close all

AssemblyRate.rate = [];
neuron.rate = [];

sessDirs = {'on10','off10','on30','off30'};
% sessDirs = {'off10','off30'};
armNames = {'rateReturn','rateReturn_L','rateReturn_R','rateBase',...
    'rateDelay','rateStem','rateChoice','rateReward'};
armNames2 = {'eventStrengthReturn','eventStrengthBase',...
    'eventStrengthDelay','eventStrengthStem','eventStrengthChoice','eventStrengthReward'};

xaxisName = {'on10 on&off','on10 on','off10 on&off','off10 off',...
    'on30 on&off','on30 on','off30 on&off','off30 off'};
       
for j = 1:length(sessDirs)
    for k = 1:length(armNames)
        % arm compare value
        AssemblyRate.(sessDirs{j}).(armNames{k}) = [];
    end
    
    for k = 1:length(armNames2)
        % arm compare value
        ArmEventStrength.(sessDirs{j}).(armNames2{k}) = [];
    end
    
    rateMatrix.(sessDirs{j}) = [];
    AssemblyRate.sleep1 = [];
    AssemblyRate.sleep2 = [];
    ArmEventStrength.sleep1 = [];
    ArmEventStrength.sleep2 = [];
end

% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    Fig8ArmRateFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8TreadmillArmRate_SameDelay_DelayAssembly.mat');
    load(Fig8ArmRateFile);
    
%     idx = find(SpikeProp.AvgRate.Fig8Rate>0.1 & SpikeProp.AvgRate.Fig8Rate<10);
%     ArmRate.rate = [ArmRate.rate,SpikeProp.AvgRate.Fig8Rate];
    
    for j = 1:length(sessDirs)
        % rate comparison
        for k = 1:length(armNames)
            AssemblyRate.(sessDirs{j}).(armNames{k}) = [AssemblyRate.(sessDirs{j}).(armNames{k}),...
                Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).(armNames{k})];
        end
        
        for k = 1:length(armNames2)
            ArmEventStrength.(sessDirs{j}).(armNames2{k}) = [ArmEventStrength.(sessDirs{j}).(armNames2{k}),...
                Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).(armNames2{k})];
        end
        
        
        rateMatrix.(sessDirs{j}) = [AssemblyRate.(sessDirs{j}).rateReturn',AssemblyRate.(sessDirs{j}).rateDelay',...
            AssemblyRate.(sessDirs{j}).rateStem'];
    end
    
    sleepDirs = sessInfo(i).sleepDirs;
    for j =  1:length(sleepDirs)    
        AssemblyRate.(sleepDirs{j}) = [AssemblyRate.(sleepDirs{j}),Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sleepDirs{j}).rateSleep];
        ArmEventStrength.(sleepDirs{j}) = [ArmEventStrength.(sleepDirs{j}),Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sleepDirs{j}).eventStrengthSleep];
    end
end

figure(1)
boxplot(AssemblyRate.on10.rateDelay','positions',1,'Colors','m','Symbol','');
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
boxplot(ArmEventStrength.on10.eventStrengthDelay','positions',1,'Colors','m','Symbol','');
hold on
boxplot(ArmEventStrength.off10.eventStrengthDelay','positions',2,'Colors','m','Symbol','');
boxplot(ArmEventStrength.on30.eventStrengthDelay','positions',3,'Colors','m','Symbol','');
boxplot(ArmEventStrength.off30.eventStrengthDelay','positions',4,'Colors','m','Symbol','');

boxplot(ArmEventStrength.on10.eventStrengthStem','positions',6,'Colors','k','Symbol','');
boxplot(ArmEventStrength.off10.eventStrengthStem','positions',7,'Colors','k','Symbol','');
boxplot(ArmEventStrength.on30.eventStrengthStem','positions',8,'Colors','k','Symbol','');
boxplot(ArmEventStrength.off30.eventStrengthStem','positions',9,'Colors','k','Symbol','');

boxplot(ArmEventStrength.on10.eventStrengthChoice','positions',11,'Colors','r','Symbol','');
boxplot(ArmEventStrength.off10.eventStrengthChoice','positions',12,'Colors','r','Symbol','');
boxplot(ArmEventStrength.on30.eventStrengthChoice','positions',13,'Colors','r','Symbol','');
boxplot(ArmEventStrength.off30.eventStrengthChoice','positions',14,'Colors','r','Symbol','');

boxplot(ArmEventStrength.on10.eventStrengthReward','positions',16,'Colors','g','Symbol','');
boxplot(ArmEventStrength.off10.eventStrengthReward','positions',17,'Colors','g','Symbol','');
boxplot(ArmEventStrength.on30.eventStrengthReward','positions',18,'Colors','g','Symbol','');
boxplot(ArmEventStrength.off30.eventStrengthReward','positions',19,'Colors','g','Symbol','');

boxplot(ArmEventStrength.on10.eventStrengthReturn','positions',21,'Colors','b','Symbol','');
boxplot(ArmEventStrength.off10.eventStrengthReturn','positions',22,'Colors','b','Symbol','');
boxplot(ArmEventStrength.on30.eventStrengthReturn','positions',23,'Colors','b','Symbol','');
boxplot(ArmEventStrength.off30.eventStrengthReturn','positions',24,'Colors','b','Symbol','');

boxplot(ArmEventStrength.sleep1','positions',26,'Colors','k','Symbol','');
boxplot(ArmEventStrength.sleep2','positions',27,'Colors','k','Symbol','');

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