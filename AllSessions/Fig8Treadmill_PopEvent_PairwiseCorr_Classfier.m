function Fig8Treadmill_PopEvent_PairwiseCorr_Classfier(inFile,AnalyzeSes)

close all
p.savePlot = 0;
p.writeToFile = 0;

% % Read in input information
sessInfo = SessInfoImport(inFile);
dayCount = 0;

for i = AnalyzeSes(1:end)
    
    sessDirs = sessInfo(i).sessDirs;
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);

%     % load population file
%     reactFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_WholeSes.mat');
%     load(reactFile);
%     % load reactivation delay specific file
%     reactDelayFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_Delay.mat');
%     load(reactDelayFile);
    
    % initiate the data
    PopEvent_CellPairwiseCorr_Classifier.rat = sessInfo(i).animal;
    PopEvent_CellPairwiseCorr_Classifier.day = sessInfo(i).day;
    
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    if sum(rateLabel) >= 20
        reactFile = fullfile(sessInfo(i).mainDir,'Cell Property','PopEvent_CellPairwiseCorr.mat');
        load(reactFile);
        
        dayCount = dayCount + 1;
        
        delayDirection = [PopEvent_CellPairwiseCorr.off10_1.delay_DirectionLabel;PopEvent_CellPairwiseCorr.off30_1.delay_DirectionLabel;...
            PopEvent_CellPairwiseCorr.off10_2.delay_DirectionLabel;PopEvent_CellPairwiseCorr.off30_2.delay_DirectionLabel];
        eventNum = length(delayDirection);
        
        if eventNum >= 20
            delayCellPairCorrTemp = [PopEvent_CellPairwiseCorr.off10_1.delay_pop_CellPairCorr,PopEvent_CellPairwiseCorr.off30_1.delay_pop_CellPairCorr,...
                PopEvent_CellPairwiseCorr.off10_2.delay_pop_CellPairCorr,PopEvent_CellPairwiseCorr.off30_2.delay_pop_CellPairCorr];
            
            nanIdx = sum(~isnan(delayCellPairCorrTemp),2)>(eventNum/10);
            delayCellPairCorr = delayCellPairCorrTemp(nanIdx,:);
            %         delayCellPairCorr(isnan(delayCellPairCorr)) = 0;
            
            
            corr_L = nanmean(delayCellPairCorr(:,delayDirection==1),2);
            corr_R = nanmean(delayCellPairCorr(:,delayDirection==0),2);
            %         figure
            %         plot([corr_L(corr_L>corr_R)';corr_R(corr_L>corr_R)'])
            %         figure
            %         plot([corr_L(corr_L<corr_R)';corr_R(corr_L<corr_R)'])
            
            figure
            Violin(corr_L,1);
            Violin(corr_R,2);
            Violin(corr_L-corr_R,3);
            
            figure
            plot([corr_L((corr_L-corr_R)>0.2)';corr_R((corr_L-corr_R)>0.2)'],'o-')
            xlim([0 3])
            Text = sprintf('%2.2f%s',100*sum((corr_L-corr_R)>0.2)/length(corr_R),'%');
            text(1.3,0.6,Text,'FontSize',10);
            set(gca, 'XTick', [1,2], 'XTickLabel', {'Left','Right'});
            %         ylim([0 1])
            
            figure
            plot([corr_L((corr_R-corr_L)>0.2)';corr_R((corr_R-corr_L)>0.2)'],'o-')
            xlim([0 3])
            Text = sprintf('%2.2f%s',100*sum((corr_R-corr_L)>0.2)/length(corr_R),'%');
            text(1.3,0.6,Text,'FontSize',10);
            set(gca, 'XTick', [1,2], 'XTickLabel', {'Left','Right'});
            %         ylim([0 1])
            
            % start svm fitting
            %         SVMModel = fitcsvm(delayCellPairCorr',delayDirection,'KernelFunction','rbf','Standardize',true,'ClassNames',{'1','0'});
        end
        
        if p.writeToFile == 1
            save(fullfile(savedir2,'PopEvent_CellPairwiseCorr.mat'), 'PopEvent_CellPairwiseCorr');
        end
        clear PopEvent_CellPairwiseCorr
        close all
    end
    fprintf('Finished analysis for session %d\n',i)
end

end
