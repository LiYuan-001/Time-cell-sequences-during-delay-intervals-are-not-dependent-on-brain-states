function Fig8TreadmillArmRate_Assembly_OnOff_Quant(inFile,AnalyzeSes)
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
        AssemblyRate.on.(sessDirs{j}).(armNames{k}) = [];
        AssemblyRate.off.(sessDirs{j}).(armNames{k}) = [];
    end
    
    for k = 1:length(armNames2)
        % arm compare value
        ArmEventStrength.on.(sessDirs{j}).(armNames2{k}) = [];
        ArmEventStrength.off.(sessDirs{j}).(armNames2{k}) = [];
    end
    
    rateMatrix.(sessDirs{j}) = [];
    AssemblyRate.on.sleep1 = [];
    AssemblyRate.on.sleep2 = [];
    ArmEventStrength.on.sleep1 = [];
    ArmEventStrength.on.sleep2 = [];
    
    AssemblyRate.off.sleep1 = [];
    AssemblyRate.off.sleep2 = [];
    ArmEventStrength.off.sleep1 = [];
    ArmEventStrength.off.sleep2 = [];
    
end

% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    Fig8ArmRateFile = fullfile(sessInfo(i).mainDir,'Cell Property','ArmRate_SameDelay_Assembly_OnOff-25ms.mat');
    load(Fig8ArmRateFile);
    
%     idx = find(SpikeProp.AvgRate.Fig8Rate>0.1 & SpikeProp.AvgRate.Fig8Rate<10);
%     ArmRate.rate = [ArmRate.rate,SpikeProp.AvgRate.Fig8Rate];
    
    for j = 1:length(sessDirs)
        % rate comparison
        for k = 1:length(armNames)
            if isfield(ArmRate_SameDelay_Assembly_OnOff,'on')
                AssemblyRate.on.(sessDirs{j}).(armNames{k}) = [AssemblyRate.on.(sessDirs{j}).(armNames{k}),...
                    ArmRate_SameDelay_Assembly_OnOff.on.(sessDirs{j}).(armNames{k})];
            end
            if isfield(ArmRate_SameDelay_Assembly_OnOff,'off')
                AssemblyRate.off.(sessDirs{j}).(armNames{k}) = [AssemblyRate.off.(sessDirs{j}).(armNames{k}),...
                    ArmRate_SameDelay_Assembly_OnOff.off.(sessDirs{j}).(armNames{k})];
            end
        end
        
        for k = 1:length(armNames2)
            if isfield(ArmRate_SameDelay_Assembly_OnOff,'on')
                ArmEventStrength.on.(sessDirs{j}).(armNames2{k}) = [ArmEventStrength.on.(sessDirs{j}).(armNames2{k}),...
                    ArmRate_SameDelay_Assembly_OnOff.on.(sessDirs{j}).(armNames2{k})];
            end
            if isfield(ArmRate_SameDelay_Assembly_OnOff,'off')
                ArmEventStrength.off.(sessDirs{j}).(armNames2{k}) = [ArmEventStrength.off.(sessDirs{j}).(armNames2{k}),...
                    ArmRate_SameDelay_Assembly_OnOff.off.(sessDirs{j}).(armNames2{k})];
            end
        end
        
        
%         rateMatrix.(sessDirs{j}) = [AssemblyRate.(sessDirs{j}).rateReturn',AssemblyRate.(sessDirs{j}).rateDelay',...
%             AssemblyRate.(sessDirs{j}).rateStem'];
    end
    
    sleepDirs = sessInfo(i).sleepDirs;
    for j =  1:length(sleepDirs)
        if isfield(ArmRate_SameDelay_Assembly_OnOff,'on')
            AssemblyRate.on.(sleepDirs{j}) = [AssemblyRate.on.(sleepDirs{j}),ArmRate_SameDelay_Assembly_OnOff.on.(sleepDirs{j}).rateSleep];
            ArmEventStrength.on.(sleepDirs{j}) = [ArmEventStrength.on.(sleepDirs{j}),ArmRate_SameDelay_Assembly_OnOff.on.(sleepDirs{j}).eventStrengthSleep];
        end
        if isfield(ArmRate_SameDelay_Assembly_OnOff,'off')
            AssemblyRate.off.(sleepDirs{j}) = [AssemblyRate.off.(sleepDirs{j}),ArmRate_SameDelay_Assembly_OnOff.off.(sleepDirs{j}).rateSleep];
            ArmEventStrength.off.(sleepDirs{j}) = [ArmEventStrength.off.(sleepDirs{j}),ArmRate_SameDelay_Assembly_OnOff.off.(sleepDirs{j}).eventStrengthSleep];
        end
    end
end

figure(1)
subplot(2,1,1)
boxplot(AssemblyRate.on.on10.rateDelay','positions',1,'Colors','m','Symbol','');
hold on
boxplot(AssemblyRate.on.off10.rateDelay','positions',2,'Colors','m','Symbol','');
boxplot(AssemblyRate.on.on30.rateDelay','positions',3,'Colors','m','Symbol','');
boxplot(AssemblyRate.on.off30.rateDelay','positions',4,'Colors','m','Symbol','');

boxplot(AssemblyRate.on.on10.rateStem','positions',6,'Colors','k','Symbol','');
boxplot(AssemblyRate.on.off10.rateStem','positions',7,'Colors','k','Symbol','');
boxplot(AssemblyRate.on.on30.rateStem','positions',8,'Colors','k','Symbol','');
boxplot(AssemblyRate.on.off30.rateStem','positions',9,'Colors','k','Symbol','');

boxplot(AssemblyRate.on.on10.rateChoice','positions',11,'Colors','r','Symbol','');
boxplot(AssemblyRate.on.off10.rateChoice','positions',12,'Colors','r','Symbol','');
boxplot(AssemblyRate.on.on30.rateChoice','positions',13,'Colors','r','Symbol','');
boxplot(AssemblyRate.on.off30.rateChoice','positions',14,'Colors','r','Symbol','');

boxplot(AssemblyRate.on.on10.rateReward','positions',16,'Colors','g','Symbol','');
boxplot(AssemblyRate.on.off10.rateReward','positions',17,'Colors','g','Symbol','');
boxplot(AssemblyRate.on.on30.rateReward','positions',18,'Colors','g','Symbol','');
boxplot(AssemblyRate.on.off30.rateReward','positions',19,'Colors','g','Symbol','');

boxplot(AssemblyRate.on.on10.rateReturn','positions',21,'Colors','b','Symbol','');
boxplot(AssemblyRate.on.off10.rateReturn','positions',22,'Colors','b','Symbol','');
boxplot(AssemblyRate.on.on30.rateReturn','positions',23,'Colors','b','Symbol','');
boxplot(AssemblyRate.on.off30.rateReturn','positions',24,'Colors','b','Symbol','');

boxplot(AssemblyRate.on.sleep1','positions',26,'Colors','k','Symbol','');
boxplot(AssemblyRate.on.sleep2','positions',27,'Colors','k','Symbol','');

xlim([0 28])
set(gca,'XTick',[1,6,11,16,21,26,27],'XTickLabel',{'Delay','Stem','Choice','Reward','Return','Sleep1','Sleep2'});
title('Assembly event average rate on')

subplot(2,1,2)
boxplot(AssemblyRate.off.on10.rateDelay','positions',1,'Colors','m','Symbol','');
hold on
boxplot(AssemblyRate.off.off10.rateDelay','positions',2,'Colors','m','Symbol','');
boxplot(AssemblyRate.off.on30.rateDelay','positions',3,'Colors','m','Symbol','');
boxplot(AssemblyRate.off.off30.rateDelay','positions',4,'Colors','m','Symbol','');

boxplot(AssemblyRate.off.on10.rateStem','positions',6,'Colors','k','Symbol','');
boxplot(AssemblyRate.off.off10.rateStem','positions',7,'Colors','k','Symbol','');
boxplot(AssemblyRate.off.on30.rateStem','positions',8,'Colors','k','Symbol','');
boxplot(AssemblyRate.off.off30.rateStem','positions',9,'Colors','k','Symbol','');

boxplot(AssemblyRate.off.on10.rateChoice','positions',11,'Colors','r','Symbol','');
boxplot(AssemblyRate.off.off10.rateChoice','positions',12,'Colors','r','Symbol','');
boxplot(AssemblyRate.off.on30.rateChoice','positions',13,'Colors','r','Symbol','');
boxplot(AssemblyRate.off.off30.rateChoice','positions',14,'Colors','r','Symbol','');

boxplot(AssemblyRate.off.on10.rateReward','positions',16,'Colors','g','Symbol','');
boxplot(AssemblyRate.off.off10.rateReward','positions',17,'Colors','g','Symbol','');
boxplot(AssemblyRate.off.on30.rateReward','positions',18,'Colors','g','Symbol','');
boxplot(AssemblyRate.off.off30.rateReward','positions',19,'Colors','g','Symbol','');

boxplot(AssemblyRate.off.on10.rateReturn','positions',21,'Colors','b','Symbol','');
boxplot(AssemblyRate.off.off10.rateReturn','positions',22,'Colors','b','Symbol','');
boxplot(AssemblyRate.off.on30.rateReturn','positions',23,'Colors','b','Symbol','');
boxplot(AssemblyRate.off.off30.rateReturn','positions',24,'Colors','b','Symbol','');

boxplot(AssemblyRate.off.sleep1','positions',26,'Colors','k','Symbol','');
boxplot(AssemblyRate.off.sleep2','positions',27,'Colors','k','Symbol','');

xlim([0 28])
set(gca,'XTick',[1,6,11,16,21,26,27],'XTickLabel',{'Delay','Stem','Choice','Reward','Return','Sleep1','Sleep2'});
title('Assembly event average rate off')

figure(2)
subplot(2,1,1)
boxplot(ArmEventStrength.on.on10.eventStrengthDelay','positions',1,'Colors','m','Symbol','');
hold on
boxplot(ArmEventStrength.on.off10.eventStrengthDelay','positions',2,'Colors','m','Symbol','');
boxplot(ArmEventStrength.on.on30.eventStrengthDelay','positions',3,'Colors','m','Symbol','');
boxplot(ArmEventStrength.on.off30.eventStrengthDelay','positions',4,'Colors','m','Symbol','');

boxplot(ArmEventStrength.on.on10.eventStrengthStem','positions',6,'Colors','k','Symbol','');
boxplot(ArmEventStrength.on.off10.eventStrengthStem','positions',7,'Colors','k','Symbol','');
boxplot(ArmEventStrength.on.on30.eventStrengthStem','positions',8,'Colors','k','Symbol','');
boxplot(ArmEventStrength.on.off30.eventStrengthStem','positions',9,'Colors','k','Symbol','');

boxplot(ArmEventStrength.on.on10.eventStrengthChoice','positions',11,'Colors','r','Symbol','');
boxplot(ArmEventStrength.on.off10.eventStrengthChoice','positions',12,'Colors','r','Symbol','');
boxplot(ArmEventStrength.on.on30.eventStrengthChoice','positions',13,'Colors','r','Symbol','');
boxplot(ArmEventStrength.on.off30.eventStrengthChoice','positions',14,'Colors','r','Symbol','');

boxplot(ArmEventStrength.on.on10.eventStrengthReward','positions',16,'Colors','g','Symbol','');
boxplot(ArmEventStrength.on.off10.eventStrengthReward','positions',17,'Colors','g','Symbol','');
boxplot(ArmEventStrength.on.on30.eventStrengthReward','positions',18,'Colors','g','Symbol','');
boxplot(ArmEventStrength.on.off30.eventStrengthReward','positions',19,'Colors','g','Symbol','');

boxplot(ArmEventStrength.on.on10.eventStrengthReturn','positions',21,'Colors','b','Symbol','');
boxplot(ArmEventStrength.on.off10.eventStrengthReturn','positions',22,'Colors','b','Symbol','');
boxplot(ArmEventStrength.on.on30.eventStrengthReturn','positions',23,'Colors','b','Symbol','');
boxplot(ArmEventStrength.on.off30.eventStrengthReturn','positions',24,'Colors','b','Symbol','');

boxplot(ArmEventStrength.on.sleep1','positions',26,'Colors','k','Symbol','');
boxplot(ArmEventStrength.on.sleep2','positions',27,'Colors','k','Symbol','');

xlim([0 28])
set(gca,'XTick',[1,6,11,16,21,26,27],'XTickLabel',{'Delay','Stem','Choice','Reward','Return','Sleep1','Sleep2'});
title('Assembly event average strength on')

subplot(2,1,2)
boxplot(ArmEventStrength.off.on10.eventStrengthDelay','positions',1,'Colors','m','Symbol','');
hold on
boxplot(ArmEventStrength.off.off10.eventStrengthDelay','positions',2,'Colors','m','Symbol','');
boxplot(ArmEventStrength.off.on30.eventStrengthDelay','positions',3,'Colors','m','Symbol','');
boxplot(ArmEventStrength.off.off30.eventStrengthDelay','positions',4,'Colors','m','Symbol','');

boxplot(ArmEventStrength.off.on10.eventStrengthStem','positions',6,'Colors','k','Symbol','');
boxplot(ArmEventStrength.off.off10.eventStrengthStem','positions',7,'Colors','k','Symbol','');
boxplot(ArmEventStrength.off.on30.eventStrengthStem','positions',8,'Colors','k','Symbol','');
boxplot(ArmEventStrength.off.off30.eventStrengthStem','positions',9,'Colors','k','Symbol','');

boxplot(ArmEventStrength.off.on10.eventStrengthChoice','positions',11,'Colors','r','Symbol','');
boxplot(ArmEventStrength.off.off10.eventStrengthChoice','positions',12,'Colors','r','Symbol','');
boxplot(ArmEventStrength.off.on30.eventStrengthChoice','positions',13,'Colors','r','Symbol','');
boxplot(ArmEventStrength.off.off30.eventStrengthChoice','positions',14,'Colors','r','Symbol','');

boxplot(ArmEventStrength.off.on10.eventStrengthReward','positions',16,'Colors','g','Symbol','');
boxplot(ArmEventStrength.off.off10.eventStrengthReward','positions',17,'Colors','g','Symbol','');
boxplot(ArmEventStrength.off.on30.eventStrengthReward','positions',18,'Colors','g','Symbol','');
boxplot(ArmEventStrength.off.off30.eventStrengthReward','positions',19,'Colors','g','Symbol','');

boxplot(ArmEventStrength.off.on10.eventStrengthReturn','positions',21,'Colors','b','Symbol','');
boxplot(ArmEventStrength.off.off10.eventStrengthReturn','positions',22,'Colors','b','Symbol','');
boxplot(ArmEventStrength.off.on30.eventStrengthReturn','positions',23,'Colors','b','Symbol','');
boxplot(ArmEventStrength.off.off30.eventStrengthReturn','positions',24,'Colors','b','Symbol','');

boxplot(ArmEventStrength.off.sleep1','positions',26,'Colors','k','Symbol','');
boxplot(ArmEventStrength.off.sleep2','positions',27,'Colors','k','Symbol','');

xlim([0 28])
set(gca,'XTick',[1,6,11,16,21,26,27],'XTickLabel',{'Delay','Stem','Choice','Reward','Return','Sleep1','Sleep2'});
title('Assembly event average strength off')

%    title({'Delay','Stem-Delay/Stem+Delay'})
%    set(gca,'XTick',[1,2,5,6],'XTickLabel',{'All trials Ctr','Lesion','Correct Ctr','Lesion'});
%    [H_Stem2Delayall,P_Stem2Delayall]=kstest2(CtrThres.Stem2Delayall',LesionThres.Stem2Delayall');
%    [H_Stem2Delayc,P_Stem2Delayc]=kstest2(CtrThres.Stem2Delayc',LesionThres.Stem2Delayc');

end