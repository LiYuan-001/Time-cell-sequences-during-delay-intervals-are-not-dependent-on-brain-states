% it needs to be stable in each session
% and stable in 2 session combined
%
function Fig8TreadmillTimeCellQuant_Selectivity_Block(inFile,AnalyzeSes)

close all

p.savePlot = 1;
p.writeToFile = 1;
p.avgRateThres = 0.5;

avgRate = [];
fieldLabel.allSession = [];
allClusterNum = 0;

% % Read in input information
sessInfo = SessInfoImport(inFile);

if p.savePlot
    % directory for plot figures
    % generate a folder for each rat eah day under the current folder
    savedir = sprintf('%s%s',cd,'\Figures\Quantification Figures\Time cell num-dist quant');
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end
end
    
delay_onIdx = [];
rate_Delay = [];
for i = AnalyzeSes(1:end)
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load time field file
    timeFieldFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayField_TrialbyTrial_2Session.mat');
    load(timeFieldFile);
    % load delayFire map
    delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap_2Session.mat');
    load(delayFile);
    % load correlation values
    stabilityFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayFireStability_TrialbyTrial_2Session.mat');
    load(stabilityFile);
    % load arm rate file
    armRateFile = fullfile(sessInfo(i).mainDir,'Cell Property', 'Fig8TreadmillArmRate.mat');
    load(armRateFile);
    
    
    clusterNum = length(SpikeProp.AvgRate.Fig8Rate);
    idx = find(SpikeProp.AvgRate.Fig8Rate>0.1 & SpikeProp.AvgRate.Fig8Rate<5);
    avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
        
    % get each phase names (no delay etc)
    if length(sessInfo(1).sessDirs) == 8
        sessDirs = {'on10','off10','on30','off30'};
    elseif contains(sessInfo(1).sessDirs{1},'on')
        sessDirs = {'on10','on30'};
    else
        sessDirs = {'off10','off30'};
    end
    
 
    for j = 1:length(sessDirs)
        if i == AnalyzeSes(1)
            fieldLabel_Def1.(sessDirs{j}) = [];
            fieldLabel_Def2.(sessDirs{j}) = [];
            timeMap_Def1.(sessDirs{j}) = [];
            timeMap_Def2.(sessDirs{j}) = [];
            pre_Map_Def1.(sessDirs{j}) = [];
            post_Map_Def1.(sessDirs{j}) = [];
            corr_Def1.(sessDirs{j}) = [];
            corr_Def2.(sessDirs{j}) = [];
            endField_1.(sessDirs{j}) = [];
            endField_2.(sessDirs{j}) = [];
            corrLabel_Def1.(sessDirs{j}) = [];
            corrLabel_Def2.(sessDirs{j}) = [];
        end
        
        fieldLabel_Def1.(sessDirs{j}) = [fieldLabel_Def1.(sessDirs{j}),DelayField_TrialbyTrial_2Session.(sessDirs{j}).timeField_1(idx)];
        fieldLabel_Def2.(sessDirs{j}) = [fieldLabel_Def2.(sessDirs{j}),DelayField_TrialbyTrial_2Session.(sessDirs{j}).timeField_2(idx)];
        
        % average correlation value of the 2 blocks
        corr_Def1.(sessDirs{j}) = [corr_Def1.(sessDirs{j}),DelayFireStability_TrialbyTrial_2Session.(sessDirs{j}).delayCorr_Def1_Trial(idx)]; 
        corr_Def2.(sessDirs{j}) = [corr_Def2.(sessDirs{j}),DelayFireStability_TrialbyTrial_2Session.(sessDirs{j}).delayCorr_Def2_Trial(idx)]; 
        
        corrLabel_Def1.(sessDirs{j}) = [corrLabel_Def1.(sessDirs{j}),DelayFireStability_TrialbyTrial_2Session.(sessDirs{j}).delayCorr_Def1_Trial_Shuffle_Label95(idx)];
        corrLabel_Def2.(sessDirs{j}) = [corrLabel_Def2.(sessDirs{j}),DelayFireStability_TrialbyTrial_2Session.(sessDirs{j}).delayCorr_Def2_Trial_Shuffle_Label95(idx)];
        
        sess1 = strcat(sessDirs{j},'_1');
        sess2 = strcat(sessDirs{j},'_2');

%         blockFieldTemp = DelayField_TrialbyTrial_Shuffle.(sess1).timeField_1(idx) + DelayField_TrialbyTrial_Shuffle.(sess2).timeField_1(idx);
%         block_field_Def1.(sessDirs{j}) = [block_field_Def1.(sessDirs{j}),blockFieldTemp];
%         blockFieldTemp = DelayField_TrialbyTrial_Shuffle.(sess1).timeField_2(idx) + DelayField_TrialbyTrial_Shuffle.(sess2).timeField_2(idx);
%         block_field_Def2.(sessDirs{j}) = [block_field_Def2.(sessDirs{j}),blockFieldTemp];
        
        for k = idx
            timeMap_Def1.(sessDirs{j}) = [timeMap_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate1_Combined_Smooth{k}];
            timeMap_Def2.(sessDirs{j}) = [timeMap_Def2.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate2_Combined_Smooth{k}];

            pre_Map_Def1.(sessDirs{j}) = [pre_Map_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).pre_spikeRate1_Combined_Smooth{k}];
            post_Map_Def1.(sessDirs{j}) = [post_Map_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).post_mapMean_Def1_2Session{k}(1:20)];
            
            if sum(DelayField_TrialbyTrial_2Session.(sessDirs{j}).fieldLabel_1{k}(end-6:end))>0
                endField_1.(sessDirs{j}) = [endField_1.(sessDirs{j}),1];
            else
                endField_1.(sessDirs{j}) = [endField_1.(sessDirs{j}),0];
            end
            
            if sum(DelayField_TrialbyTrial_2Session.(sessDirs{j}).fieldLabel_2{k}(end-6:end))>0
                endField_2.(sessDirs{j}) = [endField_2.(sessDirs{j}),1];
            else
                endField_2.(sessDirs{j}) = [endField_2.(sessDirs{j}),0];
            end
            
        end        
    end 
    
    rate_DelayTemp = [Fig8TreadmillArmRate.on10.rateDelay(idx,1)';Fig8TreadmillArmRate.off10.rateDelay(idx,1)';...
        Fig8TreadmillArmRate.on30.rateDelay(idx,1)';Fig8TreadmillArmRate.off30.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];
    delay_onIdxTemp = (sum(rate_DelayTemp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
    
end

close all

%% on10, off10, on30, off30
% all time limited cells
h1 = figure;
h1.Position = [100 100 1200 400];
% on30 off30
on30Idx_Time = (fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.on30(on30Idx_Time,:);

[~,peakIdx] = max(cellMapTemp1,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = cellMapTemp1(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);

subplot(1,4,1)
imagesc(cellMapTemp1_Sort_Norm(:,1:66))
colormap(jet)
title('on 30 0-10 s','Interpreter','None')
axis on
set(gca, 'xtick', [1,66]);
set(gca, 'xticklabels', [0 10]);
clim([0 1])
subplot(1,4,2)
imagesc(cellMapTemp1_Sort_Norm(:,67:end))
colormap(jet)
title('on 30 10-30 s','Interpreter','None')
axis on
set(gca, 'xtick', [1,200-66]);
set(gca, 'xticklabels', [10 30]);
clim([0 1])

on30_block1_selectivity = max(cellMapTemp1_Sort_Norm(:,1:66),[],2)./nanmean(cellMapTemp1_Sort_Norm(:,1:66),2);
on30_block2_selectivity = max(cellMapTemp1_Sort_Norm(:,67:end),[],2)./nanmean(cellMapTemp1_Sort_Norm(:,67:end),2);
figure(2)
Violin(on30_block1_selectivity,1,'ShowData',false)
Violin(on30_block2_selectivity,2,'ShowData',false)

off30Idx_Time = (fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.off30(off30Idx_Time,:);

[~,peakIdx] = max(cellMapTemp1,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = cellMapTemp1(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);

figure(1)
subplot(1,4,3)
imagesc(cellMapTemp1_Sort_Norm(:,1:66))
colormap(jet)
title('off 30 0-10 s','Interpreter','None')
axis on
set(gca, 'xtick', [1,66]);
set(gca, 'xticklabels', [0 10]);
clim([0 1])
subplot(1,4,4)
imagesc(cellMapTemp1_Sort_Norm(:,67:end))
colormap(jet)
title('off 30 10-30 s','Interpreter','None')
axis on
set(gca, 'xtick', [1,200-66]);
set(gca, 'xticklabels', [10 30]);
clim([0 1])

off30_block1_selectivity = max(cellMapTemp1_Sort_Norm(:,1:66),[],2)./nanmean(cellMapTemp1_Sort_Norm(:,1:66),2);
off30_block2_selectivity = max(cellMapTemp1_Sort_Norm(:,67:end),[],2)./nanmean(cellMapTemp1_Sort_Norm(:,67:end),2);
figure(2)
Violin(off30_block1_selectivity,3,'ShowData',false)
Violin(off30_block2_selectivity,4,'ShowData',false)

% time cells noemalized by all max
% h2 = figure(2);
% h2.Position = [100 100 1200 400];
%% on10, off10, on30, off30
% all time cells
h2 = figure;
h2.Position = [100 100 1200 400];
% on10 off10
on10Idx_Time = (fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
on10_mean = mean(cellMapTemp1,2);

off10Idx_Time = (fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx);
cellMapTemp2 = timeMap_Def1.off10(off10Idx_Time,:);
off10_mean = mean(cellMapTemp2,2);

plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),...
    1,[1,4,1],[1,4,2],{'Delay start at barrier','Stable cells on10'},{'-3','0','10','13'},{'Stable cells off10'},{'-3','0','10','13'});

% on30 and off30
on30Idx_Time = (fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.on30(on30Idx_Time,:);
on30_mean = mean(cellMapTemp1,2);

off30Idx_Time = (fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx);
cellMapTemp2 = timeMap_Def1.off30(off30Idx_Time,:);
off30_mean = mean(cellMapTemp2,2);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on30(on30Idx_Time,:),post_Map_Def1.on30(on30Idx_Time,:),pre_Map_Def1.off30(off30Idx_Time,:),post_Map_Def1.off30(off30Idx_Time,:),...
    1,[1,4,3],[1,4,4],{'Delay start at barrier','Stable cells on10'},{'-3','0','30','33'},{'Stable cells off10'},{'-3','0','30','33'});



end

function rgb = pixelcolour(map,maxrate,cmap)
cmap1 = ...
    [    0         0    0.5625; ...
         0         0    0.6875; ...
         0         0    0.8125; ...
         0         0    0.9375; ...
         0    0.0625    1.0000; ...
         0    0.1875    1.0000; ...
         0    0.3125    1.0000; ...
         0    0.4375    1.0000; ...
         0    0.5625    1.0000; ...
         0    0.6875    1.0000; ...
         0    0.8125    1.0000; ...
         0    0.9375    1.0000; ...
    0.0625    1.0000    1.0000; ...
    0.1875    1.0000    0.8750; ...
    0.3125    1.0000    0.7500; ...
    0.4375    1.0000    0.6250; ...
    0.5625    1.0000    0.5000; ...
    0.6875    1.0000    0.3750; ...
    0.8125    1.0000    0.2500; ...
    0.9375    1.0000    0.1250; ...
    1.0000    1.0000         0; ...
    1.0000    0.8750         0; ...
    1.0000    0.7500         0; ...
    1.0000    0.6250         0; ...
    1.0000    0.5000         0; ...
    1.0000    0.3750         0; ...
    1.0000    0.2500         0; ...
    1.0000    0.1250         0; ...
    1.0000         0         0; ...
    0.8750         0         0; ...
    0.7500         0         0; ...
    0.6250         0         0 ];

cmap2 = ...
   [0.0417         0         0; ...
    0.1250         0         0; ...
    0.2083         0         0; ...
    0.2917         0         0; ...
    0.3750         0         0; ...
    0.4583         0         0; ...
    0.5417         0         0; ...
    0.6250         0         0; ...
    0.7083         0         0; ...
    0.7917         0         0; ...
    0.8750         0         0; ...
    0.9583         0         0; ...
    1.0000    0.0417         0; ...
    1.0000    0.1250         0; ...
    1.0000    0.2083         0; ...
    1.0000    0.2917         0; ...
    1.0000    0.3750         0; ...
    1.0000    0.4583         0; ...
    1.0000    0.5417         0; ...
    1.0000    0.6250         0; ...
    1.0000    0.7083         0; ...
    1.0000    0.7917         0; ...
    1.0000    0.8750         0; ...
    1.0000    0.9583         0; ...
    1.0000    1.0000    0.0625; ...
    1.0000    1.0000    0.1875; ...
    1.0000    1.0000    0.3125; ...
    1.0000    1.0000    0.4375; ...
    1.0000    1.0000    0.5625; ...
    1.0000    1.0000    0.6875; ...
    1.0000    1.0000    0.8125; ...
    1.0000    1.0000    0.9375];
if strcmp(cmap,'jet')
   steps = (31*(map/maxrate))+1;
   steps = round(steps);
   if steps>32; steps = 32; end
   if steps<1; steps = 1; end
   rgb = cmap1(steps,:);
else
   steps = (31*(map/maxrate))+1;
   steps = round(steps);
   if steps>32; steps = 32; end
   if steps<1; steps = 1; end
   rgb = cmap2(steps,:);
end    
end