function Fig8Treadmill_DelayActive_Quant(inFile,AnalyzeSes)

close all

p.savePlot = 0;
p.writeToFile = 0;

p.avgRateThres = 0.5;

avgRate = [];
allClusterNum = 0;

% % Read in input information
sessInfo = SessInfoImport(inFile);
% get each phase names (no delay etc)
sessDirs = {'on10','off10','on30','off30'};
    
if p.savePlot
    % directory for plot figures
    % generate a folder for each rat eah day under the current folder
    savedir = sprintf('%s%s',cd,'\Figures\Quantification Figures\Time cell num-dist quant');
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end
end

for i = AnalyzeSes(1:end)
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load arm rate file
    armRateFile = fullfile(sessInfo(i).mainDir,'Cell Property', 'Fig8TreadmillArmRate.mat');
    load(armRateFile);
    % load delayFire map
    delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap_2Session.mat');
    load(delayFile);
    
    clusterNum = length(SpikeProp.AvgRate.Fig8Rate);
    allClusterNum = allClusterNum + clusterNum;
    
    avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
    delayMeanTemp_Def1 = zeros(clusterNum,length(sessDirs));
    delayMeanTemp_Def2 = zeros(clusterNum,length(sessDirs));
    
    for j = 1:length(sessDirs)
        if i == AnalyzeSes(1)
            timeMap_Def1.(sessDirs{j}) = [];
            timeMap_Def2.(sessDirs{j}) = [];
            pre_timeMap_Def1.(sessDirs{j}) = [];
            pre_timeMap_Def2.(sessDirs{j}) = [];
        end
        
        cellCount = 0;
        for k = 1:clusterNum
            cellCount = cellCount+1;
            
            timeMap_Def1.(sessDirs{j}) = [timeMap_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate1_Combined_Smooth{k}];
            timeMap_Def2.(sessDirs{j}) = [timeMap_Def2.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate2_Combined_Smooth{k}];
            
            pre_timeMap_Def1.(sessDirs{j}) = [pre_timeMap_Def1.(sessDirs{j});DelayFire.(sessDirs{j}).pre_delay_spikeRate1_Combined_Smooth{k}];
            pre_timeMap_Def2.(sessDirs{j}) = [pre_timeMap_Def2.(sessDirs{j});DelayFire.(sessDirs{j}).pre_delay_spikeRate2_Combined_Smooth{k}];
            
            delayMeanTemp_Def2(cellCount,j) = mean(DelayFire.(sessDirs{j}).spikeRate2_Combined{k});
        end
        delayMeanTemp_Def1(:,j) = ThetaCellRate.(sessDirs{j}).rate_All;
    end  
    
    rate_DelayTemp = [Fig8TreadmillArmRate.on10.rateDelay(idx,1)';Fig8TreadmillArmRate.off10.rateDelay(idx,1)';...
        Fig8TreadmillArmRate.on30.rateDelay(idx,1)';Fig8TreadmillArmRate.off30.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];

end
pyrIdx = avgRate>0.1 & avgRate<5;

%% delay definition 1: delay aligned by barrier
% as long as it is over average rate threshold in one session
delay_onIdx = sum(rate_Delay>p.avgRateThres,2) & pyrIdx';

delay_onIdx_1 = sum(rate_Delay(:,1)>p.avgRateThres,2) & pyrIdx';
delay_onIdx_2 = sum(rate_Delay(:,2)>p.avgRateThres,2) & pyrIdx';
delay_onIdx_3 = sum(rate_Delay(:,3)>p.avgRateThres,2) & pyrIdx';
delay_onIdx_4 = sum(rate_Delay(:,4)>p.avgRateThres,2) & pyrIdx';

% plot wenn of the cell distribution

cellMapTemp1 = timeMap_Def1.on10_1(delay_onIdx_1,:) + timeMap_Def1.on10_2(delay_onIdx_1,:);
cellMapTemp2 = timeMap_Def1.off10_1(delay_onIdx_2,:) + timeMap_Def1.off10_2(delay_onIdx_2,:);
cellMapTemp3 = timeMap_Def1.on30_1(delay_onIdx_3,:) + timeMap_Def1.on30_2(delay_onIdx_3,:);
cellMapTemp4 = timeMap_Def1.off30_1(delay_onIdx_4,:) + timeMap_Def1.off30_2(delay_onIdx_4,:);

pre_cellMapTemp1 = pre_timeMap_Def1.on10_1(delay_onIdx_1,:) + pre_timeMap_Def1.on10_2(delay_onIdx_1,:);
pre_cellMapTemp2 = pre_timeMap_Def1.off10_1(delay_onIdx_2,:) + pre_timeMap_Def1.off10_2(delay_onIdx_2,:);
pre_cellMapTemp3 = pre_timeMap_Def1.on30_1(delay_onIdx_3,:) + pre_timeMap_Def1.on30_2(delay_onIdx_3,:);
pre_cellMapTemp4 = pre_timeMap_Def1.off30_1(delay_onIdx_4,:) + pre_timeMap_Def1.off30_2(delay_onIdx_4,:);

preBinLength = size(pre_cellMapTemp1,2) - size(cellMapTemp1,2);
maxRate = max([max(max(cellMapTemp1)),max(max(cellMapTemp2)),max(max(cellMapTemp3)),max(max(cellMapTemp4))]);


figure(1)
% find peak and sort
subplot(1,4,1)
[~,peakIdx] = max(cellMapTemp1,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp_Sort = pre_cellMapTemp1(peak_Sort,:);
cellMapTemp_Sort_Norm = cellMapTemp_Sort./max(cellMapTemp1(peak_Sort,:),[],2);

imagesc(cellMapTemp_Sort_Norm);hold on
colormap(jet)
vertplot(preBinLength-0.5, 0, size(cellMapTemp_Sort_Norm,1), 'k--');
TITLE1 = 'on10 delay-on cells Normalized';
TITLE2 = 'Delay aligned by barrier';
title({TITLE1;TITLE2},'Interpreter','None')
axis on
set(gca, 'xtick', [0.5 preBinLength-0.5 size(cellMapTemp_Sort_Norm,2)-0.5]);
set(gca, 'xticklabels', [-3 0 10]);
caxis(gca,[0 1])


subplot(1,4,2)
% find peak and sort
[~,peakIdx] = max(cellMapTemp2,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp_Sort = pre_cellMapTemp2(peak_Sort,:);
cellMapTemp_Sort_Norm = cellMapTemp_Sort./max(cellMapTemp2(peak_Sort,:),[],2);

imagesc(cellMapTemp_Sort_Norm);hold on
colormap(jet)
vertplot(preBinLength-0.5, 0, size(cellMapTemp_Sort_Norm,1), 'k--');
TITLE1 = 'off10 delay-on cells Normalized';
TITLE2 = 'Delay aligned by barrier';
title({TITLE1;TITLE2},'Interpreter','None')
axis on
set(gca, 'xtick', [0.5 preBinLength-0.5 size(cellMapTemp_Sort_Norm,2)-0.5]);
set(gca, 'xticklabels', [-3 0 10]);
caxis(gca,[0 1])


subplot(1,4,3)
% find peak and sort
[~,peakIdx] = max(cellMapTemp3,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp_Sort = pre_cellMapTemp3(peak_Sort,:);
cellMapTemp_Sort_Norm = cellMapTemp_Sort./max(cellMapTemp3(peak_Sort,:),[],2);

imagesc(cellMapTemp_Sort_Norm);hold on
colormap(jet)
vertplot(preBinLength-0.5, 0, size(cellMapTemp_Sort_Norm,1), 'k--');
TITLE1 = 'on30 delay-on cells Normalized';
TITLE2 = 'Delay aligned by barrier';
title({TITLE1;TITLE2},'Interpreter','None')
axis on
set(gca, 'xtick', [0.5 preBinLength-0.5 size(cellMapTemp_Sort_Norm,2)-0.5]);
set(gca, 'xticklabels', [-3 0 30]);
caxis([0 1])

subplot(1,4,4)
% find peak and sort
[~,peakIdx] = max(cellMapTemp4,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp_Sort = pre_cellMapTemp4(peak_Sort,:);
cellMapTemp_Sort_Norm = cellMapTemp_Sort./max(cellMapTemp4(peak_Sort,:),[],2);

imagesc(cellMapTemp_Sort_Norm);hold on
colormap(jet)
vertplot(preBinLength-0.5, 0, size(cellMapTemp_Sort_Norm,1), 'k--');
TITLE1 = 'off30 delay-on cells Normalized';
TITLE2 = 'Delay aligned by barrier';
title({TITLE1;TITLE2},'Interpreter','None')
axis on
set(gca, 'xtick', [0.5 preBinLength-0.5 size(cellMapTemp_Sort_Norm,2)-0.5]);
set(gca, 'xticklabels', [-3 0 30]);
caxis([0 1])

figure(2)
% find peak and sort
subplot(1,4,1)
[~,peakIdx] = max(cellMapTemp1,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp_Sort = pre_cellMapTemp1(peak_Sort,:);
% cellMapTemp_Sort_Norm = cellMapTemp_Sort./max(cellMapTemp1(peak_Sort,:),[],2);

imagesc(cellMapTemp_Sort);hold on
colormap(jet)
vertplot(preBinLength-0.5, 0, size(cellMapTemp_Sort_Norm,1), 'k--');
TITLE1 = 'on10 delay-on cells';
TITLE2 = 'Delay aligned by barrier';
title({TITLE1;TITLE2},'Interpreter','None')
axis on
set(gca, 'xtick', [0.5 preBinLength-0.5 size(cellMapTemp_Sort_Norm,2)-0.5]);
set(gca, 'xticklabels', [-3 0 10]);
caxis(gca,[0 maxRate/3])


subplot(1,4,2)
% find peak and sort
[~,peakIdx] = max(cellMapTemp2,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp_Sort = pre_cellMapTemp2(peak_Sort,:);
% cellMapTemp_Sort_Norm = cellMapTemp_Sort./max(cellMapTemp_Sort,[],2);

imagesc(cellMapTemp_Sort);hold on
colormap(jet)
vertplot(preBinLength-0.5, 0, size(cellMapTemp_Sort_Norm,1), 'k--');
TITLE1 = 'off10 delay-on cells';
TITLE2 = 'Delay aligned by barrier';
title({TITLE1;TITLE2},'Interpreter','None')
axis on
set(gca, 'xtick', [0.5 preBinLength-0.5 size(cellMapTemp_Sort_Norm,2)-0.5]);
set(gca, 'xticklabels', [-3 0 10]);
caxis([0 maxRate/3])


subplot(1,4,3)
% find peak and sort
[~,peakIdx] = max(cellMapTemp3,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp_Sort = pre_cellMapTemp3(peak_Sort,:);
% cellMapTemp_Sort_Norm = cellMapTemp_Sort./max(cellMapTemp_Sort,[],2);

imagesc(cellMapTemp_Sort);hold on
colormap(jet)
vertplot(preBinLength-0.5, 0, size(cellMapTemp_Sort_Norm,1), 'k--');
TITLE1 = 'on30 delay-on cells';
TITLE2 = 'Delay aligned by barrier';
title({TITLE1;TITLE2},'Interpreter','None')
axis on
set(gca, 'xtick', [0.5 preBinLength-0.5 size(cellMapTemp_Sort_Norm,2)-0.5]);
set(gca, 'xticklabels', [-3 0 30]);
caxis([0 maxRate/3])

subplot(1,4,4)
% find peak and sort
[~,peakIdx] = max(cellMapTemp4,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp_Sort = pre_cellMapTemp4(peak_Sort,:);
% cellMapTemp_Sort_Norm = cellMapTemp_Sort./max(cellMapTemp_Sort,[],2);

imagesc(cellMapTemp_Sort);hold on
colormap(jet)
vertplot(preBinLength-0.5, 0, size(cellMapTemp_Sort_Norm,1), 'k--');
TITLE1 = 'off30 delay-on cells';
TITLE2 = 'Delay aligned by barrier';
title({TITLE1;TITLE2},'Interpreter','None')
axis on
set(gca, 'xtick', [0.5 preBinLength-0.5 size(cellMapTemp_Sort_Norm,2)-0.5]);
set(gca, 'xticklabels', [-3 0 30]);
caxis([0 maxRate/3])

% delayon_Rate.on10 = mean(cellMapTemp1,2);
% delayon_Rate.off10 = mean(cellMapTemp2,2);
% delayon_Rate.on30 = mean(cellMapTemp3,2);
% delayon_Rate.off30 = mean(cellMapTemp4,2);
% 
% delayon_Rate_Early.on10 = mean(cellMapTemp1(:,1:33),2);
% delayon_Rate_Early.off10 = mean(cellMapTemp2(:,1:33),2);
% delayon_Rate_Early.on30 = mean(cellMapTemp3(:,1:33),2);
% delayon_Rate_Early.off30 = mean(cellMapTemp4(:,1:33),2);
% 
% delayon_Rate_Late.on10 = mean(cellMapTemp1(:,34:end),2);
% delayon_Rate_Late.off10 = mean(cellMapTemp2(:,34:end),2);
% delayon_Rate_Late.on30 = mean(cellMapTemp3(:,34:end),2);
% delayon_Rate_Late.off30 = mean(cellMapTemp4(:,34:end),2);


%%
if p.writeToFile == 1
    save(fullfile(savedir,'timecellQuant_Timebin.mat'), 'timecellQuant_Timebin');
end
clear timecellQuant

end