function Fig8Treadmill_PairwiseCorr_Quant(inFile,AnalyzeSes)

close all

p.savePlot = 0;
p.writeToFile = 0;

% % Read in input information
sessInfo = SessInfoImport(inFile);

if p.savePlot
    % directory for plot figures
    % generate a folder for each rat eah day under the current folder
    savedir = sprintf('%s%s',cd,'\Figures\Quantification Figures\Cell paiwise corr');
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end
end

for i = AnalyzeSes(1:end)
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load pairwise correlations
    pairCorrFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellPairwiseCorr.mat');
    load(pairCorrFile);
    
    clusterNum = length(SpikeProp.AvgRate.Fig8Rate);
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    % take out putative pyr cell pair index
    pairSize = nchoosek(clusterNum,2);
    pyrPairLabel = zeros(pairSize,1);
    pairCount = 0;
    for kk = 1:clusterNum-1
        for mm = kk+1:clusterNum
            pairCount = pairCount + 1;
            if rateLabel(kk)==1 && rateLabel(mm)==1              
                pyrPairLabel(pairCount) = 1;
            end
        end
    end
    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    for j = 1:length(sessDirs)
        if i == AnalyzeSes(1)
            delay_pairWiseCorr.(sessDirs{j}) = [];
            reward_pairWiseCorr.(sessDirs{j}) = [];
        end
        delay_pairWiseCorr.(sessDirs{j}) = [delay_pairWiseCorr.(sessDirs{j});nanmean(CellPairwiseCorr.(sessDirs{j}).delay_pairWiseCorr(pyrPairLabel==1,:),2)];
        reward_pairWiseCorr.(sessDirs{j}) = [reward_pairWiseCorr.(sessDirs{j});nanmean(CellPairwiseCorr.(sessDirs{j}).reward_pairWiseCorr(pyrPairLabel==1,:),2)];
        
    end  
end

figure(1)


% TITLE1 = 'on10 delay-on cells Normalized';
% TITLE2 = 'Delay aligned by barrier';
% title({TITLE1;TITLE2},'Interpreter','None')
% axis on
% set(gca, 'xtick', [0.5 preBinLength-0.5 size(cellMapTemp_Sort_Norm,2)-0.5]);
% set(gca, 'xticklabels', [-3 0 10]);
% caxis(gca,[0 1])


%%
if p.writeToFile == 1
    save(fullfile(savedir,'timecellQuant_Timebin.mat'), 'timecellQuant_Timebin');
end
clear timecellQuant

end