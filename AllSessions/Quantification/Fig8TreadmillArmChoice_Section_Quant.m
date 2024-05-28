
function Fig8TreadmillArmChoice_Section_Quant(inFile,AnalyzeSes)
close all
figureCount = 0;
ArmRateChoice.rate = [];
p.avgRateThres = 0.5;

sessDirs = {'on10','off10','on30','off30'};
allClusterNum = 0;
avgRate = [];
rate_Delay = [];
delay_onIdx = [];


for j = 1:length(sessDirs)
        % arm compare value
        ArmRateChoice.(sessDirs{j}).delayRateLRDiffAll = [];
        % 95% value
        ArmRateChoice.(sessDirs{j}).shf_delayRateLRDiffAll95 = [];
end
% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    Fig8ArmChoice = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8TreadmillArmChoice_Section.mat');
    load(Fig8ArmChoice);  
    % load arm rate file
    armRateFile = fullfile(sessInfo(i).mainDir,'Cell Property', 'Fig8TreadmillArmRate.mat');
    load(armRateFile); 
    % load time field file
    timeFieldFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayTimeField_Trial_Shuffle.mat');
    load(timeFieldFile);
    
    clusterNum = length(SpikeProp.AvgRate.Fig8Rate);
    idx = find(SpikeProp.AvgRate.Fig8Rate>0 & SpikeProp.AvgRate.Fig8Rate<5);
    avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
    allClusterNum = allClusterNum + sum(idx);
    
    % get each phase names (no delay etc)
    sessDirs2 = sessInfo(i).sessDirs;
    
    for j = 1:length(sessDirs)
        % rate comparison
        ArmRateChoice.(sessDirs{j}).delayRateLRDiffAll = [ArmRateChoice.(sessDirs{j}).delayRateLRDiffAll;...
                Fig8TreadmillArmChoice_Section.(sessDirs{j}).delayRateLRDiffAll(idx,:)];                   
            
        % 95% value        
        ArmRateChoice.(sessDirs{j}).shf_delayRateLRDiffAll95 = [ArmRateChoice.(sessDirs{j}).shf_delayRateLRDiffAll95;...
                Fig8TreadmillArmChoice_Section.(sessDirs{j}).shf_delayRateLRDiffAll95(idx,:)];

    end
    
    rate_Temp = [Fig8TreadmillArmRate.on10.rateDelay(idx,1)';Fig8TreadmillArmRate.off10.rateDelay(idx,1)';...
        Fig8TreadmillArmRate.on30.rateDelay(idx,1)';Fig8TreadmillArmRate.off30.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_Temp];
    delay_onIdxTemp = (sum(rate_Temp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
     
    
end

% timeCellIdxTemp = (fieldLabel_Def1.on10_1 + fieldLabel_Def1.on10_2 + fieldLabel_Def1.off10_1 + fieldLabel_Def1.off10_2 ...
%     + fieldLabel_Def1.on30_1 + fieldLabel_Def1.on30_2 + fieldLabel_Def1.off30_1 + fieldLabel_Def1.off30_2);
% pyrIdxTemp = (avgRate>0 & avgRate<5);
%     
% timeCellIdx = timeCellIdxTemp>0 & pyrIdxTemp>0;
% non_timeCellIdx = timeCellIdxTemp==0 & pyrIdxTemp>0;

for j = 1:length(sessDirs)
    % All trials
    rateComp = ArmRateChoice.(sessDirs{j}).delayRateLRDiffAll;
    shfValue = ArmRateChoice.(sessDirs{j}).shf_delayRateLRDiffAll95;
    
    SigArmChoice.(sessDirs{j}).delayRateLRDiffAll = ((abs(rateComp)-abs(shfValue))>0);
    delay_sec_sig.(sessDirs{j}) = delay_onIdx' & SigArmChoice.(sessDirs{j}).delayRateLRDiffAll;
end
    


    figure
    plot(1:2,100*sum(delay_sec_sig.on10)/sum(delay_onIdx),'r^-');
    hold on
    plot(1:2,100*sum(delay_sec_sig.off10)/sum(delay_onIdx),'ko-');
    plot(1:6,100*sum(delay_sec_sig.on30)/sum(delay_onIdx),'r^--');
    plot(1:6,100*sum(delay_sec_sig.off30)/sum(delay_onIdx),'ko--');

    
    figure
    subplot(1,4,1)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.on10.return/100,'FaceColor',[1,0,0.2]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('Return arm random dist-on 10')
    subplot(1,4,2)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.off10.return/100,'FaceColor',[0.6,0.6,0.6]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('off 10')
    subplot(1,4,3)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.on30.return/100,'FaceColor',[1,0,0.2]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('on 30')
    subplot(1,4,4)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.off30.return/100,'FaceColor',[0.6,0.6,0.6]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('off 30')
    
    figure
    subplot(1,4,1)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.on10.delay/100,'FaceColor',[1,0,0.2]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('Delay arm random dist-on 10')
    subplot(1,4,2)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.off10.delay/100,'FaceColor',[0.6,0.6,0.6]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('off 10')
    subplot(1,4,3)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.on30.delay/100,'FaceColor',[1,0,0.2]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('on 30')
    subplot(1,4,4)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.off30.delay/100,'FaceColor',[0.6,0.6,0.6]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('off 30')
    
    
    figure
    subplot(1,4,1)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.on10.stem/100,'FaceColor',[1,0,0.2]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('Stem arm random dist-on 10')
    subplot(1,4,2)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.off10.stem/100,'FaceColor',[0.6,0.6,0.6]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('off 10')
    subplot(1,4,3)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.on30.stem/100,'FaceColor',[1,0,0.2]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('on 30')
    subplot(1,4,4)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.off30.stem/100,'FaceColor',[0.6,0.6,0.6]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('off 30')
    
    figure
    subplot(1,4,1)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.on10.choice/100,'FaceColor',[1,0,0.2]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('Stem arm random dist-on 10')
    subplot(1,4,2)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.off10.choice/100,'FaceColor',[0.6,0.6,0.6]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('off 10')
    subplot(1,4,3)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.on30.choice/100,'FaceColor',[1,0,0.2]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('on 30')
    subplot(1,4,4)
    bar(1:1:20,SigArmChoice_BaseLine_Dist.off30.choice/100,'FaceColor',[0.6,0.6,0.6]);
    set(gca,'view',[90 -90])
    xlim([0 20])
    ylim([0 1])
    title('off 30')
    
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
    ylim([0 1])
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