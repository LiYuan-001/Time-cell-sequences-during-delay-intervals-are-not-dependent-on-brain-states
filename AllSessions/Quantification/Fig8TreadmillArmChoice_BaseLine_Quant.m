
function Fig8TreadmillArmChoice_BaseLine_Quant(inFile,AnalyzeSes)
close all
figureCount = 0;
ArmRateChoice.rate = [];
p.avgRateThres = 0.5;

sessDirs = {'on10','off10','on30','off30'};
armCompareNames = {'delayRateLRDiffAll','returnRateLRDiffAll','stemRateLRDiffAll','ChoRateLRDiffAll'};
shuffle_armCompareNames = {'shf_delayRateLRDiffAll95','shf_returnRateLRDiffAll95','shf_stemRateLRDiffAll95','shf_ChoRateLRDiffAll95'};
allClusterNum = 0;
avgRate = [];
rate_Delay = [];
delay_onIdx = [];
rate_Return = [];
return_onIdx = [];  
rate_Stem = [];
stem_onIdx = [];
rate_Choice = [];
choice_onIdx = [];


for j = 1:length(sessDirs)
    for k = 1:length(armCompareNames)
        % arm compare value
        ArmRateChoice.(sessDirs{j}).(armCompareNames{k}) = [];
    end
    for n = 1:length(shuffle_armCompareNames)
        % 95% value
        ArmRateChoice.(sessDirs{j}).(shuffle_armCompareNames{n}) = [];
    end
end
% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    Fig8ArmChoice = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8TreadmillArmChoice_Baseline.mat');
    load(Fig8ArmChoice);  
    % load arm rate file
    armRateFile = fullfile(sessInfo(i).mainDir,'Cell Property', 'Fig8TreadmillArmRate.mat');
    load(armRateFile); 
    % load time field file
    timeFieldFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayTimeField_Trial_Shuffle.mat');
    load(timeFieldFile);
    
    clusterNum = length(SpikeProp.AvgRate.Fig8Rate);
    idx = find(SpikeProp.AvgRate.Fig8Rate>0.1 & SpikeProp.AvgRate.Fig8Rate<5);
    avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
    allClusterNum = allClusterNum + sum(idx);
    
    % get each phase names (no delay etc)
    sessDirs2 = sessInfo(i).sessDirs;
    
    for j = 1:length(sessDirs)
        Fig8TreadmillArmChoice_Baseline.on10.delayRateLRDiffAll
        % rate comparison
        ArmRateChoice.(sessDirs{j}).delayRateLRDiffAll = [ArmRateChoice.(sessDirs{j}).delayRateLRDiffAll,...
                Fig8TreadmillArmChoice_Baseline.(sessDirs{j}).delayRateLRDiffAll(idx)];            
        ArmRateChoice.(sessDirs{j}).stemRateLRDiffAll = [ArmRateChoice.(sessDirs{j}).stemRateLRDiffAll,...
                Fig8TreadmillArmChoice_Baseline.(sessDirs{j}).stemRateLRDiffAll(idx)];
        ArmRateChoice.(sessDirs{j}).returnRateLRDiffAll = [ArmRateChoice.(sessDirs{j}).returnRateLRDiffAll,...
                Fig8TreadmillArmChoice_Baseline.(sessDirs{j}).returnRateLRDiffAll(idx)];
        ArmRateChoice.(sessDirs{j}).ChoRateLRDiffAll = [ArmRateChoice.(sessDirs{j}).ChoRateLRDiffAll,...
                Fig8TreadmillArmChoice_Baseline.(sessDirs{j}).ChoRateLRDiffAll(idx)];
            
            
        % 95% value        
        ArmRateChoice.(sessDirs{j}).shf_delayRateLRDiffAll95 = [ArmRateChoice.(sessDirs{j}).shf_delayRateLRDiffAll95,...
                Fig8TreadmillArmChoice_Baseline.(sessDirs{j}).shf_delayRateLRDiffAll95(idx)];
        ArmRateChoice.(sessDirs{j}).shf_stemRateLRDiffAll95 = [ArmRateChoice.(sessDirs{j}).shf_stemRateLRDiffAll95,...
                Fig8TreadmillArmChoice_Baseline.(sessDirs{j}).shf_stemRateLRDiffAll95(idx)];        
        ArmRateChoice.(sessDirs{j}).shf_returnRateLRDiffAll95 = [ArmRateChoice.(sessDirs{j}).shf_returnRateLRDiffAll95,...
            Fig8TreadmillArmChoice_Baseline.(sessDirs{j}).shf_returnRateLRDiffAll95(idx)];
        ArmRateChoice.(sessDirs{j}).shf_ChoRateLRDiffAll95 = [ArmRateChoice.(sessDirs{j}).shf_ChoRateLRDiffAll95,...
            Fig8TreadmillArmChoice_Baseline.(sessDirs{j}).shf_ChoRateLRDiffAll95(idx)];
        
    end
    
    rate_Temp = [Fig8TreadmillArmRate.on10.rateDelay(idx,1)';Fig8TreadmillArmRate.off10.rateDelay(idx,1)';...
        Fig8TreadmillArmRate.on30.rateDelay(idx,1)';Fig8TreadmillArmRate.off30.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_Temp];
    delay_onIdxTemp = (sum(rate_Temp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
    
    rate_Temp = [Fig8TreadmillArmRate.on10.rateReturn(idx,1)';Fig8TreadmillArmRate.off10.rateReturn(idx,1)';...
        Fig8TreadmillArmRate.on30.rateReturn(idx,1)';Fig8TreadmillArmRate.off30.rateReturn(idx,1)'];
    rate_Return = [rate_Return,rate_Temp];
    return_onIdxTemp = (sum(rate_Temp>p.avgRateThres,1))>0;
    return_onIdx = [return_onIdx,return_onIdxTemp];
    
    rate_Temp = [Fig8TreadmillArmRate.on10.rateStem(idx,1)';Fig8TreadmillArmRate.off10.rateStem(idx,1)';...
        Fig8TreadmillArmRate.on30.rateStem(idx,1)';Fig8TreadmillArmRate.off30.rateStem(idx,1)'];
    rate_Stem = [rate_Stem,rate_Temp];
    stem_onIdxTemp = (sum(rate_Temp>p.avgRateThres,1))>0;
    stem_onIdx = [stem_onIdx,stem_onIdxTemp];
    
    rate_Temp = [Fig8TreadmillArmRate.on10.rateChoice(idx,1)';Fig8TreadmillArmRate.off10.rateChoice(idx,1)';...
        Fig8TreadmillArmRate.on30.rateChoice(idx,1)';Fig8TreadmillArmRate.off30.rateChoice(idx,1)'];
    rate_Choice = [rate_Choice,rate_Temp];
    choice_onIdxTemp = (sum(rate_Temp>p.avgRateThres,1))>0;
    choice_onIdx = [choice_onIdx,choice_onIdxTemp];
    
    
end

% timeCellIdxTemp = (fieldLabel_Def1.on10_1 + fieldLabel_Def1.on10_2 + fieldLabel_Def1.off10_1 + fieldLabel_Def1.off10_2 ...
%     + fieldLabel_Def1.on30_1 + fieldLabel_Def1.on30_2 + fieldLabel_Def1.off30_1 + fieldLabel_Def1.off30_2);
% pyrIdxTemp = (avgRate>0.1 & avgRate<5);
%     
% timeCellIdx = timeCellIdxTemp>0 & pyrIdxTemp>0;
% non_timeCellIdx = timeCellIdxTemp==0 & pyrIdxTemp>0;

for k = 1:length(armCompareNames)
    for j = 1:length(sessDirs)
        % All trials
        rateComp = ArmRateChoice.(sessDirs{j}).(armCompareNames{k});
        shfValue = ArmRateChoice.(sessDirs{j}).(shuffle_armCompareNames{k});
        SigArmChoice.(sessDirs{j}).(armCompareNames{k}) = ((abs(rateComp)-abs(shfValue))>0 & (abs(rateComp)~=0));
    end
end

 figure
    bar(1,100*sum(return_onIdx & SigArmChoice.on10.returnRateLRDiffAll)/length(return_onIdx),'FaceColor',[1,0,0.2]);
    hold on
    bar(2,100*sum(return_onIdx & SigArmChoice.off10.returnRateLRDiffAll)/length(return_onIdx),'FaceColor',[0.6,0.6,0.6]);
    bar(3,100*sum(return_onIdx & SigArmChoice.on30.returnRateLRDiffAll)/length(return_onIdx),'FaceColor',[1,0,0.2]);
    bar(4,100*sum(return_onIdx & SigArmChoice.off30.returnRateLRDiffAll)/length(return_onIdx),'FaceColor',[0.6,0.6,0.6]);
    
    bar(6,100*sum(delay_onIdx & SigArmChoice.on10.delayRateLRDiffAll)/length(return_onIdx),'FaceColor',[1,0,0.2]);
    hold on
    bar(7,100*sum(delay_onIdx & SigArmChoice.off10.delayRateLRDiffAll)/length(return_onIdx),'FaceColor',[0.6,0.6,0.6]);
    bar(8,100*sum(delay_onIdx & SigArmChoice.on30.delayRateLRDiffAll)/length(return_onIdx),'FaceColor',[1,0,0.2]);
    bar(9,100*sum(delay_onIdx & SigArmChoice.off30.delayRateLRDiffAll)/length(return_onIdx),'FaceColor',[0.6,0.6,0.6]);
    
    bar(11,100*sum(stem_onIdx & SigArmChoice.on10.stemRateLRDiffAll)/length(stem_onIdx),'FaceColor',[1,0,0.2]);
    hold on
    bar(12,100*sum(stem_onIdx & SigArmChoice.off10.stemRateLRDiffAll)/length(stem_onIdx),'FaceColor',[0.6,0.6,0.6]);
    bar(13,100*sum(stem_onIdx & SigArmChoice.on30.stemRateLRDiffAll)/length(stem_onIdx),'FaceColor',[1,0,0.2]);
    bar(14,100*sum(stem_onIdx & SigArmChoice.off30.stemRateLRDiffAll)/length(stem_onIdx),'FaceColor',[0.6,0.6,0.6]);
    
    bar(16,100*sum(choice_onIdx & SigArmChoice.on10.ChoRateLRDiffAll)/length(stem_onIdx),'FaceColor',[1,0,0.2]);
    hold on
    bar(17,100*sum(choice_onIdx & SigArmChoice.off10.ChoRateLRDiffAll)/length(stem_onIdx),'FaceColor',[0.6,0.6,0.6]);
    bar(18,100*sum(choice_onIdx & SigArmChoice.on30.ChoRateLRDiffAll)/length(stem_onIdx),'FaceColor',[1,0,0.2]);
    bar(19,100*sum(choice_onIdx & SigArmChoice.off30.ChoRateLRDiffAll)/length(stem_onIdx),'FaceColor',[0.6,0.6,0.6]);
    
    
    set(gca,'XTick',[1,2,3,4,6,7,8,9,11,12,13,14,16:19],'XTickLabel',{'on10','off10','on30','off30'});
    
    
    
    
    % figure
    figureCount = figureCount+1;
    figure(figureCount)
    Violin(ArmRateChoice.on10.(armCompareNames{k})',1);
    Violin(ArmRateChoice.off10.(armCompareNames{k})',3);
    Violin(ArmRateChoice.on30.(armCompareNames{k})',5);
    Violin(ArmRateChoice.off30.(armCompareNames{k})',7);
    
    title({strcat('All cell Arm Choice-',armCompareNames{k})},'Interpreter','None')
    set(gca,'XTick',[1,3,5,7],'XTickLabel',{'on10','off10','on30','off30'});
    
    % siginificant cells
    figureCount = figureCount+1;
    figure(figureCount)
    hold on
    if ~isempty(SigArmChoice.on10.(armCompareNames{k}))
        bar(1,100*sum(SigArmChoice.on10.(armCompareNames{k})(pyrIdxTemp))/sum(pyrIdxTemp));
    end
    if ~isempty(SigArmChoice.off10.(armCompareNames{k})(pyrIdxTemp))
        bar(3,100*sum(SigArmChoice.off10.(armCompareNames{k})(pyrIdxTemp))/sum(pyrIdxTemp));
    end
    if ~isempty(SigArmChoice.on30.(armCompareNames{k})(pyrIdxTemp))
        bar(5,100*sum(SigArmChoice.on30.(armCompareNames{k})(pyrIdxTemp))/sum(pyrIdxTemp));
    end
    if ~isempty(SigArmChoice.off30.(armCompareNames{k})(pyrIdxTemp))
        bar(7,100*sum(SigArmChoice.off30.(armCompareNames{k})(pyrIdxTemp))/sum(pyrIdxTemp));
    end
    ylabel('%')
    ylim([0 20])
    title({strcat('All cell Significant Arm Choice-',armCompareNames{k})},'Interpreter','None')
    set(gca,'XTick',[1,3,5,7],'XTickLabel',{'on10','off10','on30','off30'});
    
    % non-time cell
    % figure
    figureCount = figureCount+1;
    figure(figureCount)
    Violin(ArmRateChoice.on10.(armCompareNames{k})(non_timeCellIdx)',1);
    Violin(ArmRateChoice.off10.(armCompareNames{k})(non_timeCellIdx)',3);
    Violin(ArmRateChoice.on30.(armCompareNames{k})(non_timeCellIdx)',5);
    Violin(ArmRateChoice.off30.(armCompareNames{k})(non_timeCellIdx)',7);
    
    title({strcat('Non-Time cell Arm Choice-',armCompareNames{k})},'Interpreter','None')
    set(gca,'XTick',[1,3,5,7],'XTickLabel',{'on10','off10','on30','off30'});
    
    % siginificant cells
    figureCount = figureCount+1;
    figure(figureCount)
    hold on
    if ~isempty(SigArmChoice.on10.(armCompareNames{k})(non_timeCellIdx))
        bar(1,100*sum(SigArmChoice.on10.(armCompareNames{k})(non_timeCellIdx))/sum(non_timeCellIdx));
    end
    if ~isempty(SigArmChoice.off10.(armCompareNames{k})(non_timeCellIdx))
        bar(3,100*sum(SigArmChoice.off10.(armCompareNames{k})(non_timeCellIdx))/sum(non_timeCellIdx));
    end
    if ~isempty(SigArmChoice.on30.(armCompareNames{k})(non_timeCellIdx))
        bar(5,100*sum(SigArmChoice.on30.(armCompareNames{k})(non_timeCellIdx))/sum(non_timeCellIdx));
    end
    if ~isempty(SigArmChoice.off30.(armCompareNames{k})(non_timeCellIdx))
        bar(7,100*sum(SigArmChoice.off30.(armCompareNames{k})(non_timeCellIdx))/sum(non_timeCellIdx));
    end
    ylabel('%')
    ylim([0 100])
    title({strcat('Non-time cell Significant Arm Choice-',armCompareNames{k})},'Interpreter','None')
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