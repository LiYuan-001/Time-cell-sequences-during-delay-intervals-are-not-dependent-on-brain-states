function Fig8TreadmillArmRate_DelayAssembly_Quant(inFile,AnalyzeSes)
% close all
close all

AssemblyRate.rate = [];
neuron.rate = [];

sessDirs = {'on10','off10','on30','off30'};
sessDirs2 = {'on10_1','off10_1','on30_1','off30_1','on10_2','off10_2','on30_2','off30_2'};
% sessDirs = {'off10','off30'};
                
armNames = {'return_L','return_R','base_L','base_R','delay_L','delay_R','stem_L','stem_R','choice_L','choice_R','reward_L','reward_R'};
armNames2 = {'return','base','delay','stem','choice','reward'};

for j = 1:length(sessDirs)
%     for k = 1:length(armNames)
%         % arm compare value
%         AssemblyRate.(sessDirs{j}).(armNames{k}) = [];
%     end
    
    for k = 1:length(armNames2)
        % arm compare value
        ArmStrength.(sessDirs{j}).(armNames2{k}) = [];
        ArmRate.(sessDirs{j}).(armNames2{k}) = [];
    end
    
    rateMatrix.(sessDirs{j}) = [];
end
for j = 1:length(sessDirs2)
    
    for k = 1:length(armNames2)
        % arm compare value
        ArmStrength.(sessDirs2{j}).(armNames2{k}) = [];
        ArmRate.(sessDirs2{j}).(armNames2{k}) = [];
    end
    
    rateMatrix.(sessDirs2{j}) = [];
end
ArmRate.sleep1 = [];
ArmRate.sleep2 = [];
ArmStrength.sleep1 = [];
ArmStrength.sleep2 = [];

% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    Fig8ArmRateFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8TreadmillArmRate_SameDelay_DelayAssembly.mat');
    load(Fig8ArmRateFile);
    
    %     idx = find(SpikeProp.AvgRate.Fig8Rate>0.1 & SpikeProp.AvgRate.Fig8Rate<10);
    %     ArmRate.rate = [ArmRate.rate,SpikeProp.AvgRate.Fig8Rate];
    if isfield(Fig8TreadmillArmRate_SameDelay_DelayAssembly,'on10')
        for j = 1:length(sessDirs)
            %         % rate comparison
            %         for k = 1:length(armNames)
            %             AssemblyRate.(sessDirs{j}).(armNames{k}) = [AssemblyRate.(sessDirs{j}).(armNames{k}),...
            %                 Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).(armNames{k})];
            %         end
            
            for k = 1:length(armNames2)
                leftName = strcat(armNames2{k},'_L');
                rightName = strcat(armNames2{k},'_R');
                leftStrength = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).(leftName)(:,2);
                rightStrength = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).(rightName)(:,2);
                leftTime = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).(leftName)(:,3);
                rightTime = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).(rightName)(:,3);
                strengthTemp = (leftStrength + rightStrength)./(leftTime + rightTime);
                ArmStrength.(sessDirs{j}).(armNames2{k}) = [ArmStrength.(sessDirs{j}).(armNames2{k});...
                    strengthTemp];
            end
            
            %         rateMatrix.(sessDirs{j}) = [AssemblyRate.(sessDirs{j}).rateReturn',AssemblyRate.(sessDirs{j}).rateDelay',...
            %             AssemblyRate.(sessDirs{j}).rateStem'];
        end
        
        for j = 1:length(sessDirs2)
            
            for k = 1:length(armNames2)
                leftName = strcat(armNames2{k},'_L');
                rightName = strcat(armNames2{k},'_R');
                leftSpike = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs2{j}).(leftName)(:,1);
                rightSpike = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs2{j}).(rightName)(:,1);
                leftStrength = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs2{j}).(leftName)(:,2);
                rightStrength = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs2{j}).(rightName)(:,2);
                leftTime = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs2{j}).(leftName)(:,3);
                rightTime = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs2{j}).(rightName)(:,3);
                strengthTemp = (leftStrength + rightStrength)./(leftTime + rightTime);
                rateTemp = (leftSpike + rightSpike)./(leftTime + rightTime);
                ArmStrength.(sessDirs2{j}).(armNames2{k}) = [ArmStrength.(sessDirs2{j}).(armNames2{k});...
                    strengthTemp];
                ArmRate.(sessDirs2{j}).(armNames2{k}) = [ArmRate.(sessDirs2{j}).(armNames2{k});...
                    rateTemp];
            end
            
            %         rateMatrix.(sessDirs{j}) = [AssemblyRate.(sessDirs{j}).rateReturn',AssemblyRate.(sessDirs{j}).rateDelay',...
            %             AssemblyRate.(sessDirs{j}).rateStem'];
        end
        
        sleepDirs = sessInfo(i).sleepDirs;
        for j =  1:length(sleepDirs)
            %         AssemblyRate.(sleepDirs{j}) = [AssemblyRate.(sleepDirs{j}),Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sleepDirs{j}).rateSleep];
            rateTemp = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sleepDirs{j}).sleep(:,1)./Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sleepDirs{j}).sleep(:,3);
            strengthTemp = Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sleepDirs{j}).sleep(:,2)./Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sleepDirs{j}).sleep(:,3);
            ArmStrength.(sleepDirs{j}) = [ArmStrength.(sleepDirs{j});strengthTemp];
            ArmRate.(sleepDirs{j}) = [ArmRate.(sleepDirs{j});rateTemp];
        end
    end
end

figure
Violin(ArmRate.on10_1.delay,1)
Violin(ArmRate.on10_2.delay,2)
Violin(ArmRate.on30_1.delay,4)
Violin(ArmRate.on30_2.delay,5)

Violin(ArmRate.off10_1.delay,7)
Violin(ArmRate.off10_2.delay,8)
Violin(ArmRate.off30_1.delay,10)
Violin(ArmRate.off30_2.delay,11)

figure(1)

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