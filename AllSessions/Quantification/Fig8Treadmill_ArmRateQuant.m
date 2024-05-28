function Fig8Treadmill_ArmRateQuant(inFile,AnalyzeSes)

close all

p.savePlot = 0;
p.writeToFile = 0;

p.avgRateThres = 0.5;

avgRate = [];
rate_Delay = [];
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

armNames = {'rateReturn','rateDelay','rateStem','rateChoice','rateReward'};
session = {'on','off'};
timeNames = {'timeReturn','timeDelay','timeStem','timeChoice','timeReward'};
                
for j = 1:length(session)
    for k = 1:length(armNames)
        % arm compare value
        neuronRate.(session{j}).(armNames{k}) = [];
        burstRate.(session{j}).(armNames{k}) = [];
        singleSpkRate.(session{j}).(armNames{k}) = [];
        AssemblyRate.(session{j}).(armNames{k}) = [];
    end
end


for i = AnalyzeSes(1:end)
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load arm rate file
    armRateFile = fullfile(sessInfo(i).mainDir,'Cell Property', 'Fig8TreadmillArmRate.mat');
    load(armRateFile);
    % load arm rate burst file
    armRateFile = fullfile(sessInfo(i).mainDir,'Cell Property', 'Fig8TreadmillArmRate_Burst.mat');
    load(armRateFile);
    Fig8ArmRateFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8TreadmillArmRate_SameDelay_DelayAssembly.mat');
    load(Fig8ArmRateFile);
    
    clusterNum = length(SpikeProp.AvgRate.Fig8Rate);
    allClusterNum = allClusterNum + clusterNum;
    
    avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
    
    rate_DelayTemp = [Fig8TreadmillArmRate.on10.rateDelay(:,1)';Fig8TreadmillArmRate.off10.rateDelay(:,1)';...
        Fig8TreadmillArmRate.on30.rateDelay(:,1)';Fig8TreadmillArmRate.off30.rateDelay(:,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];
    
    for k = 1:length(armNames)
        spikeNumTemp = (Fig8TreadmillArmRate.on10.(armNames{k})(:,3)+ Fig8TreadmillArmRate.on30.(armNames{k})(:,3));
        timeTemp =(Fig8TreadmillArmRate.on10.(armNames{k})(:,2)+ Fig8TreadmillArmRate.on30.(armNames{k})(:,2));
        neuronRate.on.(armNames{k}) = [neuronRate.on.(armNames{k});spikeNumTemp./timeTemp];
        
        spikeNumTemp = (Fig8TreadmillArmRate.off10.(armNames{k})(:,3)+ Fig8TreadmillArmRate.off30.(armNames{k})(:,3));
        timeTemp =(Fig8TreadmillArmRate.off10.(armNames{k})(:,2)+ Fig8TreadmillArmRate.off30.(armNames{k})(:,2));
        neuronRate.off.(armNames{k}) = [neuronRate.off.(armNames{k});spikeNumTemp./timeTemp];
        
        spikeNumTemp = (Fig8TreadmillArmRate_Burst.on10.Burst.(armNames{k})(:,3)+ Fig8TreadmillArmRate_Burst.on30.Burst.(armNames{k})(:,3));
        timeTemp =(Fig8TreadmillArmRate_Burst.on10.Burst.(armNames{k})(:,2)+ Fig8TreadmillArmRate_Burst.on30.Burst.(armNames{k})(:,2));
        burstRate.on.(armNames{k}) = [burstRate.on.(armNames{k});spikeNumTemp./timeTemp];
        
        spikeNumTemp = (Fig8TreadmillArmRate_Burst.off10.Burst.(armNames{k})(:,3)+ Fig8TreadmillArmRate_Burst.off30.Burst.(armNames{k})(:,3));
        timeTemp =(Fig8TreadmillArmRate_Burst.off10.Burst.(armNames{k})(:,2)+ Fig8TreadmillArmRate_Burst.off30.Burst.(armNames{k})(:,2));
        burstRate.off.(armNames{k}) = [burstRate.off.(armNames{k});spikeNumTemp./timeTemp];
        
        spikeNumTemp = (Fig8TreadmillArmRate_Burst.on10.SingleSpk.(armNames{k})(:,3)+ Fig8TreadmillArmRate_Burst.on30.SingleSpk.(armNames{k})(:,3));
        timeTemp =(Fig8TreadmillArmRate_Burst.on10.SingleSpk.(armNames{k})(:,2)+ Fig8TreadmillArmRate_Burst.on30.SingleSpk.(armNames{k})(:,2));
        singleSpkRate.on.(armNames{k}) = [singleSpkRate.on.(armNames{k});spikeNumTemp./timeTemp];
        
        spikeNumTemp = (Fig8TreadmillArmRate_Burst.off10.SingleSpk.(armNames{k})(:,3)+ Fig8TreadmillArmRate_Burst.off30.SingleSpk.(armNames{k})(:,3));
        timeTemp =(Fig8TreadmillArmRate_Burst.off10.SingleSpk.(armNames{k})(:,2)+ Fig8TreadmillArmRate_Burst.off30.SingleSpk.(armNames{k})(:,2));
        singleSpkRate.off.(armNames{k}) = [singleSpkRate.off.(armNames{k});spikeNumTemp./timeTemp];
        
        
%         if isfield(Fig8TreadmillArmRate_SameDelay_DelayAssembly,'on10')
%             spikeNumTemp = Fig8TreadmillArmRate_SameDelay_DelayAssembly.on10.(armNames{k}).*Fig8TreadmillArmRate_SameDelay_DelayAssembly.on10.(timeNames{k}) + ...
%                 Fig8TreadmillArmRate_SameDelay_DelayAssembly.on30.(armNames{k}).*Fig8TreadmillArmRate_SameDelay_DelayAssembly.on30.(timeNames{k});
%             timeTemp = Fig8TreadmillArmRate_SameDelay_DelayAssembly.on10.(timeNames{k}) + Fig8TreadmillArmRate_SameDelay_DelayAssembly.on30.(timeNames{k});
%             AssemblyRate.on.(armNames{k}) = spikeNumTemp./timeTemp;
%             
%             spikeNumTemp = Fig8TreadmillArmRate_SameDelay_DelayAssembly.off10.(armNames{k}).*Fig8TreadmillArmRate_SameDelay_DelayAssembly.off10.(timeNames{k}) + ...
%                 Fig8TreadmillArmRate_SameDelay_DelayAssembly.off30.(armNames{k}).*Fig8TreadmillArmRate_SameDelay_DelayAssembly.off30.(timeNames{k});
%             timeTemp = Fig8TreadmillArmRate_SameDelay_DelayAssembly.off10.(timeNames{k}) + Fig8TreadmillArmRate_SameDelay_DelayAssembly.off30.(timeNames{k});
%             AssemblyRate.off.(armNames{k}) = spikeNumTemp./timeTemp;
%         end
        
    end    
end
pyrIdx = avgRate>0.1 & avgRate<5;

%% delay definition 1: delay aligned by barrier
% as long as it is over average rate threshold in one session
delay_onIdx = sum(rate_Delay>p.avgRateThres,1) & pyrIdx;

% delay active cell rate difference in treadmill on sessions

for k = 1:length(armNames)
    neuronRate_mean.on.(armNames{k}) = mean(neuronRate.on.(armNames{k})(pyrIdx));
    neuronRate_mean.off.(armNames{k}) = mean(neuronRate.off.(armNames{k})(pyrIdx));
    neuronRate_err.on.(armNames{k}) = std(neuronRate.on.(armNames{k})(pyrIdx))./sqrt(sum(pyrIdx));
    neuronRate_err.off.(armNames{k}) = std(neuronRate.off.(armNames{k})(pyrIdx))./sqrt(sum(pyrIdx));
    
    burstRate_mean.on.(armNames{k}) = mean(burstRate.on.(armNames{k})(pyrIdx));
    burstRate_mean.off.(armNames{k}) = mean(burstRate.off.(armNames{k})(pyrIdx));
    burstRate_err.on.(armNames{k}) = std(burstRate.on.(armNames{k})(pyrIdx))./sqrt(sum(pyrIdx));
    burstRate_err.off.(armNames{k}) = std(burstRate.off.(armNames{k})(pyrIdx))./sqrt(sum(pyrIdx));
    
    singleSpkRate_mean.on.(armNames{k}) = mean(singleSpkRate.on.(armNames{k})(pyrIdx));
    singleSpkRate_mean.off.(armNames{k}) = mean(singleSpkRate.off.(armNames{k})(pyrIdx));
    singleSpkRate_err.on.(armNames{k}) = std(singleSpkRate.on.(armNames{k})(pyrIdx))./sqrt(sum(pyrIdx));
    singleSpkRate_err.off.(armNames{k}) = std(singleSpkRate.off.(armNames{k})(pyrIdx))./sqrt(sum(pyrIdx));
%     
%     AssemblyRate_mean.on.(armNames{k}) = mean(AssemblyRate.on.(armNames{k}));
%     AssemblyRate_mean.off.(armNames{k}) = mean(AssemblyRate.off.(armNames{k}));
%     AssemblyRate_err.on.(armNames{k}) = std(AssemblyRate.on.(armNames{k}))./sqrt(length(AssemblyRate.on.(armNames{k})));
%     AssemblyRate_err.off.(armNames{k}) = std(AssemblyRate.off.(armNames{k}))./sqrt(length(AssemblyRate.off.(armNames{k})));
    
%     figure(1)
%     Violin(log2(neuronRate.on.(armNames{k})(pyrIdx)),2*k-1,'ViolinColor',[1,0,0],'ShowData',false)
%     hold on
%     Violin(log2(neuronRate.off.(armNames{k})(pyrIdx)),2*k-0.5,'ViolinColor',[0,0,0],'ShowData',false)
%     set(gca,'XTick',[1:2:10],'XTickLabel',armNames);
%     xlim([0 10])
%     title('Pyr Rate')
%     
    
    figure(1)
    errorbar(k,neuronRate_mean.on.(armNames{k}),neuronRate_err.on.(armNames{k}),'rd');
    hold on
    errorbar(k+0.2,neuronRate_mean.off.(armNames{k}),neuronRate_err.off.(armNames{k}),'kd');
    title('Pyr Rate')
    set(gca,'XTick',[1:5],'XTickLabel',armNames);
    xlim([0 6])
    
    
    figure(2)
    errorbar(k,burstRate_mean.on.(armNames{k}),burstRate_err.on.(armNames{k}),'rd');
    hold on
    errorbar(k+0.2,burstRate_mean.off.(armNames{k}),burstRate_err.off.(armNames{k}),'kd');
    title('Burst rate')
    set(gca,'XTick',[1:5],'XTickLabel',armNames);
    xlim([0 6])
    
    
    figure(3)
    errorbar(k,singleSpkRate_mean.on.(armNames{k}),singleSpkRate_err.on.(armNames{k}),'rd');
    hold on
    errorbar(k+0.2,singleSpkRate_mean.off.(armNames{k}),singleSpkRate_err.off.(armNames{k}),'kd');
    title('Single spike rate')
    set(gca,'XTick',[1:5],'XTickLabel',armNames);
    xlim([0 6])
%     
%     figure(4)
%     errorbar(k,AssemblyRate_mean.on.(armNames{k}),AssemblyRate_err.on.(armNames{k}),'rd');
%     hold on
%     errorbar(k+0.2,AssemblyRate_mean.off.(armNames{k}),AssemblyRate_err.off.(armNames{k}),'kd');
%     title('Assembly rate')
end

% set(gca,'XTick',[1:5],'XTickLabel',armNames);
% xlim([0 6])
% delay active cells rate difference in treadmill off sessions


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