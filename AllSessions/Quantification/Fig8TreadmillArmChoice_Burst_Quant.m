
function Fig8TreadmillArmChoice_Burst_Quant(inFile,AnalyzeSes)
close all
figureCount = 0;
BurstChoice.rate = [];

sessDirs = {'on10','off10','on30','off30'};
armCompareNames = {'delayRateLRDiffAll','stemRateLRDiffAll'};
shuffle_armCompareNames = {'shf_delayRateLRDiffAll95','shf_stemRateLRDiffAll95'};
allClusterNum = 0;
avgRate = [];

for j = 1:length(sessDirs)
    for k = 1:length(armCompareNames)
        % arm compare value
        BurstChoice.(sessDirs{j}).(armCompareNames{k}) = [];
        SingleChoice.(sessDirs{j}).(armCompareNames{k}) = [];
    end
    for n = 1:length(shuffle_armCompareNames)
        % 95% value
        BurstChoice.(sessDirs{j}).(shuffle_armCompareNames{n}) = [];
        SingleChoice.(sessDirs{j}).(shuffle_armCompareNames{n}) = [];
    end
end
% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    Fig8TreadmillArmChoice_BurstFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8TreadmillArmChoice_Burst.mat');
    load(Fig8TreadmillArmChoice_BurstFile);    
%     % load time field file
%     timeFieldFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayTimeField_Trial_Shuffle.mat');
%     load(timeFieldFile);
    
    clusterNum = length(SpikeProp.AvgRate.Fig8Rate);
    allClusterNum = allClusterNum + clusterNum;
    
    % get each phase names (no delay etc)
    sessDirs2 = sessInfo(i).sessDirs;
    avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
 
    
    for j = 1:length(sessDirs)
        Fig8TreadmillArmChoice_Burst.on10.Burst.delayRateLRDiffAll
        % rate comparison
        BurstChoice.(sessDirs{j}).delayRateLRDiffAll = [BurstChoice.(sessDirs{j}).delayRateLRDiffAll,...
                Fig8TreadmillArmChoice_Burst.(sessDirs{j}).Burst.delayRateLRDiffAll];
            
        BurstChoice.(sessDirs{j}).stemRateLRDiffAll = [BurstChoice.(sessDirs{j}).stemRateLRDiffAll,...
                Fig8TreadmillArmChoice_Burst.(sessDirs{j}).Burst.stemRateLRDiffAll];
            
            SingleChoice.(sessDirs{j}).delayRateLRDiffAll = [SingleChoice.(sessDirs{j}).delayRateLRDiffAll,...
                Fig8TreadmillArmChoice_Burst.(sessDirs{j}).SingleSpk.delayRateLRDiffAll];
            
        SingleChoice.(sessDirs{j}).stemRateLRDiffAll = [SingleChoice.(sessDirs{j}).stemRateLRDiffAll,...
                Fig8TreadmillArmChoice_Burst.(sessDirs{j}).SingleSpk.stemRateLRDiffAll];
            
        % 95% value
        
        BurstChoice.(sessDirs{j}).shf_delayRateLRDiffAll95 = [BurstChoice.(sessDirs{j}).shf_delayRateLRDiffAll95,...
                Fig8TreadmillArmChoice_Burst.(sessDirs{j}).Burst.shf_delayRateLRDiffAll95];
        BurstChoice.(sessDirs{j}).shf_stemRateLRDiffAll95 = [BurstChoice.(sessDirs{j}).shf_stemRateLRDiffAll95,...
                Fig8TreadmillArmChoice_Burst.(sessDirs{j}).Burst.shf_stemRateLRDiffAll95];      
            
            SingleChoice.(sessDirs{j}).shf_delayRateLRDiffAll95 = [SingleChoice.(sessDirs{j}).shf_delayRateLRDiffAll95,...
                Fig8TreadmillArmChoice_Burst.(sessDirs{j}).SingleSpk.shf_delayRateLRDiffAll95];
        SingleChoice.(sessDirs{j}).shf_stemRateLRDiffAll95 = [SingleChoice.(sessDirs{j}).shf_stemRateLRDiffAll95,...
                Fig8TreadmillArmChoice_Burst.(sessDirs{j}).SingleSpk.shf_stemRateLRDiffAll95]; 
    end
    
%     for j = 1:length(sessDirs2)
%         if i == AnalyzeSes(1)
%             fieldLabel_Def1.(sessDirs2{j}) = [];
%         end
%         
%         fieldLabel_Def1.(sessDirs2{j}) = [fieldLabel_Def1.(sessDirs2{j}),DelayField_TrialbyTrial_Shuffle.(sessDirs2{j}).timeField_1];
%     end
end

% timeCellIdxTemp = (fieldLabel_Def1.on10_1 + fieldLabel_Def1.on10_2 + fieldLabel_Def1.off10_1 + fieldLabel_Def1.off10_2 ...
%     + fieldLabel_Def1.on30_1 + fieldLabel_Def1.on30_2 + fieldLabel_Def1.off30_1 + fieldLabel_Def1.off30_2);
pyrIdxTemp = (avgRate>0.1 & avgRate<5);
    
% timeCellIdx = timeCellIdxTemp>0 & pyrIdxTemp>0;
% non_timeCellIdx = timeCellIdxTemp==0 & pyrIdxTemp>0;

for k = 1:length(armCompareNames)
    for j = 1:length(sessDirs)
        % All trials
        rateComp = BurstChoice.(sessDirs{j}).(armCompareNames{k});
        shfValue = BurstChoice.(sessDirs{j}).(shuffle_armCompareNames{k});
        SigArmChoice_Burst.(sessDirs{j}).(armCompareNames{k}) = (abs(rateComp)-abs(shfValue))>0;
        
        rateComp = SingleChoice.(sessDirs{j}).(armCompareNames{k});
        shfValue = SingleChoice.(sessDirs{j}).(shuffle_armCompareNames{k});
        SigArmChoice_SingleSpk.(sessDirs{j}).(armCompareNames{k}) = (abs(rateComp)-abs(shfValue))>0;
    end
    
    % figure
    figureCount = figureCount+1;
    figure(figureCount)
    Violin(BurstChoice.on10.(armCompareNames{k})',1);
    Violin(BurstChoice.off10.(armCompareNames{k})',3);
    Violin(BurstChoice.on30.(armCompareNames{k})',5);
    Violin(BurstChoice.off30.(armCompareNames{k})',7);
    
    title({strcat('All cell Arm Choice-Burst-',armCompareNames{k})},'Interpreter','None')
    set(gca,'XTick',[1,3,5,7],'XTickLabel',{'on10','off10','on30','off30'});
    
    % siginificant cells
    figureCount = figureCount+1;
    figure(figureCount)
    hold on
    if ~isempty(SigArmChoice_Burst.on10.(armCompareNames{k}))
        bar(1,100*sum(SigArmChoice_Burst.on10.(armCompareNames{k})(pyrIdxTemp))/sum(pyrIdxTemp));
    end
    if ~isempty(SigArmChoice_Burst.off10.(armCompareNames{k})(pyrIdxTemp))
        bar(3,100*sum(SigArmChoice_Burst.off10.(armCompareNames{k})(pyrIdxTemp))/sum(pyrIdxTemp));
    end
    if ~isempty(SigArmChoice_Burst.on30.(armCompareNames{k})(pyrIdxTemp))
        bar(5,100*sum(SigArmChoice_Burst.on30.(armCompareNames{k})(pyrIdxTemp))/sum(pyrIdxTemp));
    end
    if ~isempty(SigArmChoice_Burst.off30.(armCompareNames{k})(pyrIdxTemp))
        bar(7,100*sum(SigArmChoice_Burst.off30.(armCompareNames{k})(pyrIdxTemp))/sum(pyrIdxTemp));
    end
    ylabel('%')
    ylim([0 20])
    title({strcat('All cell Significant Arm Choice-Burst-',armCompareNames{k})},'Interpreter','None')
    set(gca,'XTick',[1,3,5,7],'XTickLabel',{'on10','off10','on30','off30'});
    
    % single spike
    figureCount = figureCount+1;
    figure(figureCount)
    Violin(SingleChoice.on10.(armCompareNames{k})',1);
    Violin(SingleChoice.off10.(armCompareNames{k})',3);
    Violin(SingleChoice.on30.(armCompareNames{k})',5);
    Violin(SingleChoice.off30.(armCompareNames{k})',7);
    
    title({strcat('All cell Arm Choice-',armCompareNames{k})},'Interpreter','None')
    set(gca,'XTick',[1,3,5,7],'XTickLabel',{'on10','off10','on30','off30'});
    
    % siginificant cells
    figureCount = figureCount+1;
    figure(figureCount)
    hold on
    if ~isempty(SigArmChoice_SingleSpk.on10.(armCompareNames{k}))
        bar(1,100*sum(SigArmChoice_SingleSpk.on10.(armCompareNames{k})(pyrIdxTemp))/sum(pyrIdxTemp));
    end
    if ~isempty(SigArmChoice_SingleSpk.off10.(armCompareNames{k})(pyrIdxTemp))
        bar(3,100*sum(SigArmChoice_SingleSpk.off10.(armCompareNames{k})(pyrIdxTemp))/sum(pyrIdxTemp));
    end
    if ~isempty(SigArmChoice_SingleSpk.on30.(armCompareNames{k})(pyrIdxTemp))
        bar(5,100*sum(SigArmChoice_SingleSpk.on30.(armCompareNames{k})(pyrIdxTemp))/sum(pyrIdxTemp));
    end
    if ~isempty(SigArmChoice_SingleSpk.off30.(armCompareNames{k})(pyrIdxTemp))
        bar(7,100*sum(SigArmChoice_SingleSpk.off30.(armCompareNames{k})(pyrIdxTemp))/sum(pyrIdxTemp));
    end
    ylabel('%')
    ylim([0 20])
    title({strcat('All cell Significant Arm Choice-SingleSpk-',armCompareNames{k})},'Interpreter','None')
    set(gca,'XTick',[1,3,5,7],'XTickLabel',{'on10','off10','on30','off30'});
    
    %     if ~isempty(SigArmChoice.nodelay.(armCompareNames{k}).Correct) && ~isempty(SigLayer2.Lesion.nodelay.(armCompareNames{k}).Correct)
    %     [SigLayer2_HCorrect.(armCompareNames{k}).nodelay,SigLayer2_PCorrect.(armCompareNames{k}).nodelay]=kstest2(SigArmChoice.nodelay.(armCompareNames{k}).Correct',SigLayer2.Lesion.nodelay.(armCompareNames{k}).Correct');
    %     end
    %     if ~isempty(SigArmChoice.delay10.(armCompareNames{k}).Correct) && ~isempty(SigLayer2.Lesion.delay10.(armCompareNames{k}).Correct)
    %         [SigLayer2_HCorrect.(armCompareNames{k}).delay10,SigLayer2_PCorrect.(armCompareNames{k}).delay10]=kstest2(SigArmChoice.delay10.(armCompareNames{k}).Correct',SigLayer2.Lesion.delay10.(armCompareNames{k}).Correct');
    %     end
    %     if ~isempty(SigArmChoice.delay60.(armCompareNames{k}).Correct) && ~isempty(SigLayer2.Lesion.delay60.(armCompareNames{k}).Correct)
    %         [SigLayer2_HCorrect.(armCompareNames{k}).delay60,SigLayer2_PCorrect.(armCompareNames{k}).delay60]=kstest2(SigArmChoice.delay60.(armCompareNames{k}).Correct',SigLayer2.Lesion.delay60.(armCompareNames{k}).Correct');
    %     end
    
end

end