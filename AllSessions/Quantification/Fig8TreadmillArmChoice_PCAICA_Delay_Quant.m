
function Fig8TreadmillArmChoice_PCAICA_Delay_Quant(inFile,AnalyzeSes)
close all
figureCount = 0;
ArmPCAChoice.rate = [];

sessDirs = {'on10','off10','on30','off30'};
armCompareNames = {'delayRateLRDiffAll','returnRateLRDiffAll','stemRateLRDiffAll','ChoRateLRDiffAll'};
shuffle_armCompareNames = {'shf_delayRateLRDiffAll95','shf_returnRateLRDiffAll95','shf_stemRateLRDiffAll95','shf_ChoRateLRDiffAll95'};
assemblyNum_All = 0;

for j = 1:length(sessDirs)
    for k = 1:length(armCompareNames)
        % arm compare value
        ArmPCAChoice.(sessDirs{j}).(armCompareNames{k}) = [];
        ArmPCA_Strength_Choice.(sessDirs{j}).(armCompareNames{k}) = [];
    end
    for n = 1:length(shuffle_armCompareNames)
        % 95% value
        ArmPCAChoice.(sessDirs{j}).(shuffle_armCompareNames{n}) = [];
        ArmPCA_Strength_Choice.(sessDirs{j}).(shuffle_armCompareNames{n}) = [];
    end
    Assembly_Delay_Rate.(sessDirs{j}) = [];
    Assembly_Stem_Rate.(sessDirs{j}) = [];
    Assembly_Reward_Rate.(sessDirs{j}) = [];
    Assembly_Choice_Rate.(sessDirs{j}) = [];
    AssemblyCellNum = []; 
end

delayActive = [];
returnActive = [];
stemActive = [];
choiceActive = [];

% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    
    
    Fig8ArmChoice = fullfile(sessInfo(i).mainDir,'Cell Property','AssemblyDelay_ArmChoice_SameDelay.mat');
    load(Fig8ArmChoice);
    
    if isfield(AssemblyDelay_ArmChoice_SameDelay,'on10')
        
        % load assembly file
        assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_Delay-25ms.mat');
        load(assemblyFile);
    
        assemblyRateFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8TreadmillArmRate_SameDelay_DelayAssembly.mat');
        load(assemblyRateFile);
          
%         assemblyNum = length(AssemblyDelay_ArmChoice_SameDelay.on10.delayRateLRDiffAll);
        
%         for m = 1:CellAssembly_WholeSes.patNum
%             AssemblyCellNum = [AssemblyCellNum,length(CellAssembly_WholeSes.AssmblPtrnCellIDs{m})];
%             assemblyNum_All = assemblyNum_All +1;
%             nameTemp = sprintf('%s%d%s%d%s%d','Rat-',CellAssembly_WholeSes.rat,'Day-',CellAssembly_WholeSes.day,'Assembly-',m);
%             AssemblyName{assemblyNum_All} = nameTemp;
%         end
%   
        for j = 1:length(sessDirs)
            % delay Rate
            Assembly_Delay_Rate.(sessDirs{j}) = [Assembly_Delay_Rate.(sessDirs{j}),...
                Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).rateDelay];
            Assembly_Stem_Rate.(sessDirs{j}) = [Assembly_Stem_Rate.(sessDirs{j}),...
                Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).rateStem];
            Assembly_Reward_Rate.(sessDirs{j}) = [Assembly_Reward_Rate.(sessDirs{j}),...
                Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).rateReward];
            Assembly_Choice_Rate.(sessDirs{j}) = [Assembly_Choice_Rate.(sessDirs{j}),...
                Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).rateChoice];

            % rate comparison
            ArmPCAChoice.(sessDirs{j}).delayRateLRDiffAll = [ArmPCAChoice.(sessDirs{j}).delayRateLRDiffAll,...
                AssemblyDelay_ArmChoice_SameDelay.(sessDirs{j}).delayRateLRDiffAll];
            
            ArmPCAChoice.(sessDirs{j}).stemRateLRDiffAll = [ArmPCAChoice.(sessDirs{j}).stemRateLRDiffAll,...
                AssemblyDelay_ArmChoice_SameDelay.(sessDirs{j}).stemRateLRDiffAll];
            
            ArmPCAChoice.(sessDirs{j}).ChoRateLRDiffAll = [ArmPCAChoice.(sessDirs{j}).ChoRateLRDiffAll,...
                AssemblyDelay_ArmChoice_SameDelay.(sessDirs{j}).ChoRateLRDiffAll];
            
            ArmPCAChoice.(sessDirs{j}).returnRateLRDiffAll = [ArmPCAChoice.(sessDirs{j}).returnRateLRDiffAll,...
                AssemblyDelay_ArmChoice_SameDelay.(sessDirs{j}).returnRateLRDiffAll];
            
            % 95% value
            
            ArmPCAChoice.(sessDirs{j}).shf_delayRateLRDiffAll95 = [ArmPCAChoice.(sessDirs{j}).shf_delayRateLRDiffAll95,...
                AssemblyDelay_ArmChoice_SameDelay.(sessDirs{j}).shf_delayRateLRDiffAll95];
            
            ArmPCAChoice.(sessDirs{j}).shf_stemRateLRDiffAll95 = [ArmPCAChoice.(sessDirs{j}).shf_stemRateLRDiffAll95,...
                AssemblyDelay_ArmChoice_SameDelay.(sessDirs{j}).shf_stemRateLRDiffAll95];
            
            ArmPCAChoice.(sessDirs{j}).shf_ChoRateLRDiffAll95 = [ArmPCAChoice.(sessDirs{j}).shf_ChoRateLRDiffAll95,...
                AssemblyDelay_ArmChoice_SameDelay.(sessDirs{j}).shf_stemChoRateLRDiffAll95];
            
            ArmPCAChoice.(sessDirs{j}).shf_returnRateLRDiffAll95 = [ArmPCAChoice.(sessDirs{j}).shf_returnRateLRDiffAll95,...
                AssemblyDelay_ArmChoice_SameDelay.(sessDirs{j}).shf_returnRateLRDiffAll95];
        end
        
        % delay active assemblies
        delayActive_on10 = Fig8TreadmillArmRate_SameDelay_DelayAssembly.on10.rateDelay>0.1;
        delayActive_off10 = Fig8TreadmillArmRate_SameDelay_DelayAssembly.on30.rateDelay>0.1;
        delayActive_on30 = Fig8TreadmillArmRate_SameDelay_DelayAssembly.off10.rateDelay>0.1;
        delayActive_off30 = Fig8TreadmillArmRate_SameDelay_DelayAssembly.off30.rateDelay>0.1;
        delayActiveTemp = (delayActive_on10 + delayActive_off10 + delayActive_on30 + delayActive_off30)>0;
        
        % return active assemblies
        returnActive_on10 = (Fig8TreadmillArmRate_SameDelay_DelayAssembly.on10.rateReturn)>0.1;
        returnActive_off10 = (Fig8TreadmillArmRate_SameDelay_DelayAssembly.off10.rateReturn)>0.1;
        returnActive_on30 = (Fig8TreadmillArmRate_SameDelay_DelayAssembly.on30.rateReturn)>0.1;
        returnActive_off30 = (Fig8TreadmillArmRate_SameDelay_DelayAssembly.off30.rateReturn)>0.1;
        returnActiveTemp  = (returnActive_on10 + returnActive_off10 + returnActive_on30 + returnActive_off30)>0;
        
        % stem active assemblies
        stemActive_on10 = (Fig8TreadmillArmRate_SameDelay_DelayAssembly.on10.rateStem )>0.1;
        stemActive_off10 = (Fig8TreadmillArmRate_SameDelay_DelayAssembly.off10.rateStem)>0.1;
        stemActive_on30 = (Fig8TreadmillArmRate_SameDelay_DelayAssembly.on30.rateStem)>0.1;
        stemActive_off30 = (Fig8TreadmillArmRate_SameDelay_DelayAssembly.off30.rateStem)>0.1;
        stemActiveTemp  = (stemActive_on10 + stemActive_off10 + stemActive_on30 + stemActive_off30)>0;
        
        % choice active assemblies
        choiceActive_on10 = Fig8TreadmillArmRate_SameDelay_DelayAssembly.on10.rateChoice>0.1;
        choiceActive_off10 = Fig8TreadmillArmRate_SameDelay_DelayAssembly.off10.rateChoice>0.1;
        choiceActive_on30 = Fig8TreadmillArmRate_SameDelay_DelayAssembly.on30.rateChoice>0.1;
        choiceActive_off30 = Fig8TreadmillArmRate_SameDelay_DelayAssembly.off30.rateChoice>0.1;
        choiceActiveTemp  = (choiceActive_on10 + choiceActive_off10 + choiceActive_on30 + choiceActive_off30)>0;
        
        delayActive = [delayActive,delayActiveTemp];
        returnActive = [returnActive,returnActiveTemp];
        stemActive = [stemActive,stemActiveTemp];
        choiceActive = [choiceActive,choiceActiveTemp];
    end
        
end



for k = 1:length(armCompareNames)
    for j = 1:length(sessDirs)
        % All trials
        rateComp = ArmPCAChoice.(sessDirs{j}).(armCompareNames{k});
        shfValue = ArmPCAChoice.(sessDirs{j}).(shuffle_armCompareNames{k});
        SigArmChoice.(sessDirs{j}).(armCompareNames{k}) = (abs(rateComp)-abs(shfValue))>0;
    end
end

 % figure
    figure
    Violin(abs(ArmPCAChoice.on10.returnRateLRDiffAll(returnActive==1))',1,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.off10.returnRateLRDiffAll(returnActive==1))',3,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.on30.returnRateLRDiffAll(returnActive==1))',5,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.off30.returnRateLRDiffAll(returnActive==1))',7,'ViolinColor',[0.2,0.2,0.2]);
    
    Violin(abs(ArmPCAChoice.on10.returnRateLRDiffAll(returnActive & SigArmChoice.on10.returnRateLRDiffAll))',2);
    Violin(abs(ArmPCAChoice.off10.returnRateLRDiffAll(returnActive & SigArmChoice.off10.returnRateLRDiffAll))',4);
    Violin(abs(ArmPCAChoice.on30.returnRateLRDiffAll(returnActive & SigArmChoice.on30.returnRateLRDiffAll))',6);
    Violin(abs(ArmPCAChoice.off30.returnRateLRDiffAll(returnActive & SigArmChoice.off30.returnRateLRDiffAll))',8);
    
    title({strcat('Assembly cell Arm Choice-Return')},'Interpreter','None')
    set(gca,'XTick',[1,3,5,7],'XTickLabel',{'on10','off10','on30','off30'});

    figure
    Violin(abs(ArmPCAChoice.on10.delayRateLRDiffAll(delayActive==1))',1,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.off10.delayRateLRDiffAll(delayActive==1))',3,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.on30.delayRateLRDiffAll(delayActive==1))',5,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.off30.delayRateLRDiffAll(delayActive==1))',7,'ViolinColor',[0.2,0.2,0.2]);
    
    Violin(abs(ArmPCAChoice.on10.delayRateLRDiffAll(returnActive & delayActive & SigArmChoice.on10.returnRateLRDiffAll))',2);
    Violin(abs(ArmPCAChoice.off10.delayRateLRDiffAll(returnActive & delayActive & SigArmChoice.off10.returnRateLRDiffAll))',4);
    Violin(abs(ArmPCAChoice.on30.delayRateLRDiffAll(returnActive & delayActive & SigArmChoice.on30.returnRateLRDiffAll))',6);
    Violin(abs(ArmPCAChoice.off30.delayRateLRDiffAll(returnActive & delayActive & SigArmChoice.off30.returnRateLRDiffAll))',8);
    
    title({strcat('Assembly cell Arm Choice-Delay (Return sig)')},'Interpreter','None')
    set(gca,'XTick',[1,3,5,7],'XTickLabel',{'on10','off10','on30','off30'});
    
    figure
    Violin(abs(ArmPCAChoice.on10.stemRateLRDiffAll(stemActive==1))',1,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.off10.stemRateLRDiffAll(stemActive==1))',3,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.on30.stemRateLRDiffAll(stemActive==1))',5,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.off30.stemRateLRDiffAll(stemActive==1))',7,'ViolinColor',[0.2,0.2,0.2]);
    
    Violin(abs(ArmPCAChoice.on10.stemRateLRDiffAll(stemActive & delayActive & SigArmChoice.on10.stemRateLRDiffAll))',2);
    Violin(abs(ArmPCAChoice.off10.stemRateLRDiffAll(stemActive & delayActive & SigArmChoice.off10.stemRateLRDiffAll))',4);
    Violin(abs(ArmPCAChoice.on30.stemRateLRDiffAll(stemActive & delayActive & SigArmChoice.on30.stemRateLRDiffAll))',6);
    Violin(abs(ArmPCAChoice.off30.stemRateLRDiffAll(stemActive & delayActive & SigArmChoice.off30.stemRateLRDiffAll))',8);
    
    title({strcat('Assembly cell Arm Choice-Return')},'Interpreter','None')
    set(gca,'XTick',[1,3,5,7],'XTickLabel',{'on10','off10','on30','off30'});

    figure
    Violin(abs(ArmPCAChoice.on10.delayRateLRDiffAll(delayActive==1))',1,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.off10.delayRateLRDiffAll(delayActive==1))',3,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.on30.delayRateLRDiffAll(delayActive==1))',5,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.off30.delayRateLRDiffAll(delayActive==1))',7,'ViolinColor',[0.2,0.2,0.2]);
    
    Violin(abs(ArmPCAChoice.on10.delayRateLRDiffAll(stemActive & delayActive & SigArmChoice.on10.stemRateLRDiffAll))',2);
    Violin(abs(ArmPCAChoice.off10.delayRateLRDiffAll(stemActive & delayActive & SigArmChoice.off10.stemRateLRDiffAll))',4);
    Violin(abs(ArmPCAChoice.on30.delayRateLRDiffAll(stemActive & delayActive & SigArmChoice.on30.stemRateLRDiffAll))',6);
    Violin(abs(ArmPCAChoice.off30.delayRateLRDiffAll(stemActive & delayActive & SigArmChoice.off30.stemRateLRDiffAll))',8);
    
    title({strcat('Assembly cell Arm Choice-Delay (Stem+Choice sig)')},'Interpreter','None')
    set(gca,'XTick',[1,3,5,7],'XTickLabel',{'on10','off10','on30','off30'});
    
    figure
    Violin(abs(ArmPCAChoice.on10.delayRateLRDiffAll(delayActive==1))',1,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.off10.delayRateLRDiffAll(delayActive==1))',3,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.on30.delayRateLRDiffAll(delayActive==1))',5,'ViolinColor',[0.2,0.2,0.2]);
    Violin(abs(ArmPCAChoice.off30.delayRateLRDiffAll(delayActive==1))',7,'ViolinColor',[0.2,0.2,0.2]);
    
    Violin(abs(ArmPCAChoice.on10.delayRateLRDiffAll(delayActive & SigArmChoice.on10.delayRateLRDiffAll))',2);
    Violin(abs(ArmPCAChoice.off10.delayRateLRDiffAll(delayActive & SigArmChoice.off10.delayRateLRDiffAll))',4);
    Violin(abs(ArmPCAChoice.on30.delayRateLRDiffAll(delayActive & SigArmChoice.on30.delayRateLRDiffAll))',6);
    Violin(abs(ArmPCAChoice.off30.delayRateLRDiffAll(delayActive & SigArmChoice.off30.delayRateLRDiffAll))',8);
    title({strcat('Assembly cell Arm Choice-Delay (Delay sig)')},'Interpreter','None')
    set(gca,'XTick',[1,3,5,7],'XTickLabel',{'on10','off10','on30','off30'});
    
    
    figure
    bar(1,100*sum(returnActive & SigArmChoice.on10.returnRateLRDiffAll)/sum(returnActive),'FaceColor',[1,0,0.2]);
    hold on
    bar(2,100*sum(returnActive & SigArmChoice.off10.returnRateLRDiffAll)/sum(returnActive),'FaceColor',[0.6,0.6,0.6]);
    bar(3,100*sum(returnActive & SigArmChoice.on30.returnRateLRDiffAll)/sum(returnActive),'FaceColor',[1,0,0.2]);
    bar(4,100*sum(returnActive & SigArmChoice.off30.returnRateLRDiffAll)/sum(returnActive),'FaceColor',[0.6,0.6,0.6]);
    
    bar(6,100*sum(delayActive & SigArmChoice.on10.delayRateLRDiffAll)/sum(delayActive),'FaceColor',[1,0,0.2]);
    hold on
    bar(7,100*sum(delayActive & SigArmChoice.off10.delayRateLRDiffAll)/sum(delayActive),'FaceColor',[0.6,0.6,0.6]);
    bar(8,100*sum(delayActive & SigArmChoice.on30.delayRateLRDiffAll)/sum(delayActive),'FaceColor',[1,0,0.2]);
    bar(9,100*sum(delayActive & SigArmChoice.off30.delayRateLRDiffAll)/sum(delayActive),'FaceColor',[0.6,0.6,0.6]);
    
    bar(11,100*sum(stemActive & SigArmChoice.on10.stemRateLRDiffAll)/sum(stemActive),'FaceColor',[1,0,0.2]);
    hold on
    bar(12,100*sum(stemActive & SigArmChoice.off10.stemRateLRDiffAll)/sum(stemActive),'FaceColor',[0.6,0.6,0.6]);
    bar(13,100*sum(stemActive & SigArmChoice.on30.stemRateLRDiffAll)/sum(stemActive),'FaceColor',[1,0,0.2]);
    bar(14,100*sum(stemActive & SigArmChoice.off30.stemRateLRDiffAll)/sum(stemActive),'FaceColor',[0.6,0.6,0.6]);
    
    bar(16,100*sum(choiceActive & SigArmChoice.on10.ChoRateLRDiffAll)/sum(choiceActive),'FaceColor',[1,0,0.2]);
    hold on
    bar(17,100*sum(choiceActive & SigArmChoice.off10.ChoRateLRDiffAll)/sum(choiceActive),'FaceColor',[0.6,0.6,0.6]);
    bar(18,100*sum(choiceActive & SigArmChoice.on30.ChoRateLRDiffAll)/sum(choiceActive),'FaceColor',[1,0,0.2]);
    bar(19,100*sum(choiceActive & SigArmChoice.off30.ChoRateLRDiffAll)/sum(choiceActive),'FaceColor',[0.6,0.6,0.6]);
    
    
    set(gca,'XTick',[1,2,3,4,6,7,8,9,11,12,13,14,16:19],'XTickLabel',{'on10','off10','on30','off30'});


% for k = 1:length(armCompareNames)
%     for j = 1:length(sessDirs)
%         % All trials
%         rateComp = ArmPCAChoice.(sessDirs{j}).(armCompareNames{k});
%         shfValue = ArmPCAChoice.(sessDirs{j}).(shuffle_armCompareNames{k});
%         SigArmChoice.(sessDirs{j}).(armCompareNames{k}) = (abs(rateComp)-abs(shfValue))>0;
%     end
%     
%     figureCount = figureCount+1;
%     figure(figureCount)
%     Violin(ArmPCAChoice.on10.(armCompareNames{k})',1);
%     Violin(ArmPCAChoice.off10.(armCompareNames{k})',3);
%     Violin(ArmPCAChoice.on30.(armCompareNames{k})',5);
%     Violin(ArmPCAChoice.off30.(armCompareNames{k})',7);
%     
%     % siginificant cells
%     figureCount = figureCount+1;
%     figure(figureCount)
%     hold on
%     if ~isempty(SigArmChoice.on10.(armCompareNames{k}))
%         bar(1,100*sum(SigArmChoice.on10.(armCompareNames{k}))/length(SigArmChoice.on10.(armCompareNames{k})));
%     end
%     if ~isempty(SigArmChoice.off10.(armCompareNames{k}))
%         bar(3,100*sum(SigArmChoice.off10.(armCompareNames{k}))/length(SigArmChoice.on10.(armCompareNames{k})));
%     end
%     if ~isempty(SigArmChoice.on30.(armCompareNames{k}))
%         bar(5,100*sum(SigArmChoice.on30.(armCompareNames{k}))/length(SigArmChoice.on10.(armCompareNames{k})));
%     end
%     if ~isempty(SigArmChoice.off30.(armCompareNames{k}))
%         bar(7,100*sum(SigArmChoice.off30.(armCompareNames{k}))/length(SigArmChoice.on10.(armCompareNames{k})));
%     end
%     ylabel('%')
%     ylim([0 50])
%     title({strcat('Assembly Significant Arm Choice-',armCompareNames{k})},'Interpreter','None')
%     set(gca,'XTick',[1,3,5,7],'XTickLabel',{'on10','off10','on30','off30'});
%     
%     
%     
%     %     if ~isempty(SigArmChoice.nodelay.(armCompareNames{k}).Correct) && ~isempty(SigLayer2.Lesion.nodelay.(armCompareNames{k}).Correct)
%     %     [SigLayer2_HCorrect.(armCompareNames{k}).nodelay,SigLayer2_PCorrect.(armCompareNames{k}).nodelay]=kstest2(SigArmChoice.nodelay.(armCompareNames{k}).Correct',SigLayer2.Lesion.nodelay.(armCompareNames{k}).Correct');
%     %     end
%     %     if ~isempty(SigArmChoice.delay10.(armCompareNames{k}).Correct) && ~isempty(SigLayer2.Lesion.delay10.(armCompareNames{k}).Correct)
%     %         [SigLayer2_HCorrect.(armCompareNames{k}).delay10,SigLayer2_PCorrect.(armCompareNames{k}).delay10]=kstest2(SigArmChoice.delay10.(armCompareNames{k}).Correct',SigLayer2.Lesion.delay10.(armCompareNames{k}).Correct');
%     %     end
%     %     if ~isempty(SigArmChoice.delay60.(armCompareNames{k}).Correct) && ~isempty(SigLayer2.Lesion.delay60.(armCompareNames{k}).Correct)
%     %         [SigLayer2_HCorrect.(armCompareNames{k}).delay60,SigLayer2_PCorrect.(armCompareNames{k}).delay60]=kstest2(SigArmChoice.delay60.(armCompareNames{k}).Correct',SigLayer2.Lesion.delay60.(armCompareNames{k}).Correct');
%     %     end
%     
% end
% 
% figure
% Violin(abs(Assembly_ArmChoice_SameDelay.on10.returnRateLRDiffAll),1);
% Violin(abs(Assembly_ArmChoice_SameDelay.off10.returnRateLRDiffAll),2);
% Violin(abs(Assembly_ArmChoice_SameDelay.on10.returnRateLRDiffAll),3);
% Violin(abs(Assembly_ArmChoice_SameDelay.off10.returnRateLRDiffAll),4);


end