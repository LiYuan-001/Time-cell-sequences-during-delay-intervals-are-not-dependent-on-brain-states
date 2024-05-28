% it needs to be stable in each session
% and stable in 2 session combined
%
function Fig8TreadmillTimeCellQuant_FieldLength(inFile,AnalyzeSes)

close all

p.savePlot = 0;
p.writeToFile = 0;
p.avgRateThres = 1;

avgRate = [];

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
    % average rate file
    ThetaCellRateFile = fullfile(sessInfo(i).mainDir,'Cell Property','ThetaCellRate.mat');
    load(ThetaCellRateFile);
    % load info file
    Fig8TreadmillDelayPropFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayInfoStat_2Session.mat');
    load(Fig8TreadmillDelayPropFile);
    
    
    clusterNum = length(SpikeProp.AvgRate.Fig8Rate);
    idx = find(SpikeProp.AvgRate.Fig8Rate>0.1 & SpikeProp.AvgRate.Fig8Rate<5);
    avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
        
    % get each phase names (no delay etc)
    sessDirs = {'on10','off10','on30','off30'};
    
    for j = 1:length(sessDirs)
        if i == AnalyzeSes(1)
            timefieldLabel_Def1.(sessDirs{j}) = [];
            fielMapdLabel_Def1.(sessDirs{j}) = [];
            timeMap_Def1.(sessDirs{j}) = [];
            pretimeMap_Def1.(sessDirs{j}) = [];
            endField.(sessDirs{j}) = [];
            InfoSpike.(sessDirs{j}) = [];
            InfoSec.(sessDirs{j}) = [];
            sparsity.(sessDirs{j}) = [];
            selectivity.(sessDirs{j}) = [];
            zCoherence.(sessDirs{j}) = [];
        end
        
        timefieldLabel_Def1.(sessDirs{j}) = [timefieldLabel_Def1.(sessDirs{j}),(DelayField_TrialbyTrial_2Session.(sessDirs{j}).timeField_1(idx)>0)];
        InfoSpike.(sessDirs{j}) = [InfoSpike.(sessDirs{j}),DelayInfoStat_2Session.(sessDirs{j}).info_spike(idx)];
        InfoSec.(sessDirs{j}) = [InfoSec.(sessDirs{j}),DelayInfoStat_2Session.(sessDirs{j}).info_sec(idx)];
        sparsity.(sessDirs{j}) = [sparsity.(sessDirs{j}),DelayInfoStat_2Session.(sessDirs{j}).sparsity(idx)];
        selectivity.(sessDirs{j}) = [selectivity.(sessDirs{j}),DelayInfoStat_2Session.(sessDirs{j}).selectivity(idx)];
        zCoherence.(sessDirs{j}) = [zCoherence.(sessDirs{j}),DelayInfoStat_2Session.(sessDirs{j}).zCoherence(idx)];
        
        sess1 = strcat(sessDirs{j},'_1');
        sess2 = strcat(sessDirs{j},'_2');
        if i == AnalyzeSes(1)
            rate_Block_All.(sess1) = [];
            rate_Block_All.(sess2) = [];
        end
        rate_Block_All.(sess1) = [rate_Block_All.(sess1),ThetaCellRate.(sess1).rate_All(idx)];
        rate_Block_All.(sess2) = [rate_Block_All.(sess2),ThetaCellRate.(sess2).rate_All(idx)];
        
        for k = idx
            timeMap_Def1.(sessDirs{j}) = [timeMap_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate1_Combined_Smooth{k}];
            pretimeMap_Def1.(sessDirs{j}) = [pretimeMap_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).pre_spikeRate1_Combined_Smooth{k}];
            
            fielMapdLabel_Def1.(sessDirs{j}) = [fielMapdLabel_Def1.(sessDirs{j});DelayField_TrialbyTrial_2Session.(sessDirs{j}).fieldLabel_1{k}'];
            fielLabelTemp = DelayField_TrialbyTrial_2Session.(sessDirs{j}).fieldLabel_1{k};
            if sum(fielLabelTemp(end-6:end))>0 
                endField.(sessDirs{j}) = [endField.(sessDirs{j}),1];
            else
                endField.(sessDirs{j}) = [endField.(sessDirs{j}),0];
            end
        end
        
    end     
end

close all
rate_Delay = [rate_Block_All.on10_1;rate_Block_All.on10_2;rate_Block_All.off10_1;rate_Block_All.off10_2;...
    rate_Block_All.on30_1;rate_Block_All.on30_2;rate_Block_All.off30_1;rate_Block_All.off30_2];
delay_onIdx = (sum(rate_Delay>p.avgRateThres,1))>0;

delay_onIdx_1 = sum(rate_Delay([1,2],:)>p.avgRateThres,1)>0;
delay_onIdx_2 = sum(rate_Delay([3,4],:)>p.avgRateThres,1)>0;
delay_onIdx_3 = sum(rate_Delay([5,6],:)>p.avgRateThres,1)>0;
delay_onIdx_4 = sum(rate_Delay([7,8],:)>p.avgRateThres,1)>0;

%% on10, off10, on30, off30

h1 = figure(1);
h1.Position = [100 100 1200 400];
% on10 on30
on10Idx_Time = (timefieldLabel_Def1.on10>0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
pre_cellMapTemp1 = pretimeMap_Def1.on10(on10Idx_Time,:);
[~,peakIdx] = max(cellMapTemp1,[],2);
[~,peak_Sort] = sort(peakIdx);
mapCombined1 = [pre_cellMapTemp1(peak_Sort,:),cellMapTemp1(peak_Sort,:)];
cellMap_Sort_Norm = mapCombined1./max(mapCombined1,[],2);

subplot(1,3,1)
imagesc(cellMap_Sort_Norm)
colormap(jet)
axis on
set(gca, 'xtick', [0 20 size(cellMap_Sort_Norm,2)]);
set(gca, 'xticklabels', {'-3','0','10'});
title('On10 time stable cells')
% on10_mean = mean(cellMapTemp1,2);

on30Idx_Time = (timefieldLabel_Def1.on30>0 & delay_onIdx);
cellMapTemp2 = timeMap_Def1.on30(on10Idx_Time,:);
pre_cellMapTemp2 = pretimeMap_Def1.on30(on10Idx_Time,:);
mapCombined2 = [pre_cellMapTemp2(peak_Sort,:),cellMapTemp2(peak_Sort,:)];
cellMap_Sort_Norm = mapCombined2./max(mapCombined2,[],2);

subplot(1,3,2)
imagesc(cellMap_Sort_Norm)
colormap(jet)
axis on
set(gca, 'xtick', [0 20 size(cellMap_Sort_Norm,2)]);
set(gca, 'xticklabels', {'-3','0','30'});
title('On30 cells sort by 0n 10 cells')

subplot(1,3,3)
infoSpikeTemp1 = selectivity.on10(on10Idx_Time);
infoSpikeTemp1 = infoSpikeTemp1(peak_Sort);
infoSpikeTemp2 = selectivity.on30(on10Idx_Time);
infoSpikeTemp2 = infoSpikeTemp2(peak_Sort);
plot(infoSpikeTemp1,sum(on10Idx_Time):-1:1,'r>-');
hold on
plot(infoSpikeTemp2,sum(on10Idx_Time):-1:1,'k>-');
title('Selectivity in on10&30')

figure(2)
subplot(2,3,1)
plot(InfoSpike.on10(on10Idx_Time),InfoSpike.on30(on10Idx_Time),'o');
title('Info spike in on10 time cells')
xlabel('on10')
ylabel('on30')
xlim([0 5])
x = InfoSpike.on10(on10Idx_Time);
y1 = InfoSpike.on30(on10Idx_Time);
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'r-');
plot(0:5,0:5,'k--')
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
pText = sprintf('%s%1.2f%s%1.3f','r2 = ',r_Val,' pVal = ',p_val);
text(1,4.5,pText);


subplot(2,3,2)
plot(InfoSec.on10(on10Idx_Time),InfoSec.on30(on10Idx_Time),'o');
title('Info sec in on10 time cells')
xlabel('on10')
ylabel('on30')
ylim([0 5])
x = InfoSec.on10(on10Idx_Time);
y1 = InfoSec.on30(on10Idx_Time);
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'r-');
plot(0:5,0:5,'k--')
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
pText = sprintf('%s%1.2f%s%1.3f','r2 = ',r_Val,' pVal = ',p_val);
text(1,4.5,pText);


subplot(2,3,3)
plot(sparsity.on10(on10Idx_Time),sparsity.on30(on10Idx_Time),'o');
title('Sparsity in on10 time cells')
xlabel('on10')
ylabel('on30')
x = sparsity.on10(on10Idx_Time);
y1 = sparsity.on30(on10Idx_Time);
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'r-');
plot(0:1,0:1,'k--')
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
pText = sprintf('%s%1.2f%s%1.3f','r2 = ',r_Val,' pVal = ',p_val);
text(0.1,0.8,pText);


subplot(2,3,4)
plot(selectivity.on10(on10Idx_Time),selectivity.on30(on10Idx_Time),'o');
title('Selectivity in on10 time cells')
xlabel('on10')
ylabel('on30')
xlim([0 50])
x = selectivity.on10(on10Idx_Time);
y1 = selectivity.on30(on10Idx_Time);
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'r-');
plot(0:50,0:50,'k--')
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
pText = sprintf('%s%1.2f%s%1.3f','r2 = ',r_Val,' pVal = ',p_val);
text(30,45,pText);


subplot(2,3,5)
plot(zCoherence.on10(on10Idx_Time),zCoherence.on30(on10Idx_Time),'o');
title('zCoherence in on10 time cells')
xlabel('on10')
ylabel('on30')
xlim([0 5])
ylim([0 5])
x = zCoherence.on10(on10Idx_Time);
y1 = zCoherence.on30(on10Idx_Time);
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'r-');
plot(0:5,0:5,'k--')
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
pText = sprintf('%s%1.2f%s%1.3f','r2 = ',r_Val,' pVal = ',p_val);
text(1,4.5,pText);



h3 = figure(3);
h3.Position = [100 100 1200 400];
% off10 off30
off10Idx_Time = (timefieldLabel_Def1.off10>0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.off10(off10Idx_Time,:);
off10_mean = mean(cellMapTemp1,2);

off30Idx_Time = (timefieldLabel_Def1.off30>0 & delay_onIdx);
cellMapTemp2 = timeMap_Def1.off30(off10Idx_Time,:);
off30_mean = mean(cellMapTemp2,2);

plotTimeMap2(cellMapTemp1,cellMapTemp2,2,[1,4,3],[1,4,4],...
    {'Delay start at barrier','Time cells on10'},{'0','10'},{'Cells off30 select by off10'},{'0','10'});

figure
Violin(on10_mean,1,'ShowData',false,'ViolinColor',[1,0,0]);
Violin(off10_mean,2,'ShowData',false,'ViolinColor',[0,0,0]);
Violin(on30_mean,3,'ShowData',false,'ViolinColor',[1,0,0]);
Violin(off30_mean,4,'ShowData',false,'ViolinColor',[0,0,0]);
set(gca, 'XTick', [1,2,3,4], 'XTickLabel', {'on10','off10','on30','off30'});

[h,pVal] = kstest2(on10_mean,off10_mean)
[h,pVal] = kstest2(on30_mean,off30_mean)
[h,pVal] = kstest2(on10_mean,on30_mean)
[h,pVal] = kstest2(off10_mean,off30_mean)


if p.savePlot == 1
    figName = sprintf('%s%s%s',savedir,'\Onset-Barrier-Time cell in its own session-2Session combined');
    print(figName,'-dpng','-r600');
end


figure
Violin(on10_mean,1,'ShowData',false,'ViolinColor',[1,0,0]);
Violin(off10_mean,2,'ShowData',false,'ViolinColor',[0,0,0]);
Violin(on30_mean,3,'ShowData',false,'ViolinColor',[1,0,0]);
Violin(off30_mean,4,'ShowData',false,'ViolinColor',[0,0,0]);
set(gca, 'XTick', [1,2,3,4], 'XTickLabel', {'on10','off10','on30','off30'});


%%
if p.writeToFile == 1
    save(fullfile(savedir,'timecellQuant_Timebin.mat'), 'timecellQuant_Timebin');
end
clear timecellQuant

end