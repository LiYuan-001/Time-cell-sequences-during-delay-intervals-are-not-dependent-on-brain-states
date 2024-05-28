% it needs to be stable in each session
% and stable in 2 session combined
%
function Fig8TreadmillTimeCell_RateCompareInSection(inFile,AnalyzeSes)

close all

p.savePlot = 0;
p.writeToFile = 0;
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
    sessDirs = {'on10','off10','on30','off30'};
    
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
            endField.(sessDirs{j}) = [];
        end
        
        fieldLabel_Def1.(sessDirs{j}) = [fieldLabel_Def1.(sessDirs{j}),DelayField_TrialbyTrial_2Session.(sessDirs{j}).timeField_1(idx)];
        fieldLabel_Def2.(sessDirs{j}) = [fieldLabel_Def2.(sessDirs{j}),DelayField_TrialbyTrial_2Session.(sessDirs{j}).timeField_2(idx)];
        
        % average correlation value of the 2 blocks
        corr_Def1.(sessDirs{j}) = [corr_Def1.(sessDirs{j}),DelayFireStability_TrialbyTrial_2Session.(sessDirs{j}).delayCorr_Def1_Trial(idx)]; 
        corr_Def2.(sessDirs{j}) = [corr_Def2.(sessDirs{j}),DelayFireStability_TrialbyTrial_2Session.(sessDirs{j}).delayCorr_Def2_Trial(idx)]; 
        
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
            post_Map_Def1.(sessDirs{j}) = [post_Map_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).post_mapMean_Def1_2Session{k}];
            
            if sum(DelayField_TrialbyTrial_2Session.(sessDirs{j}).fieldLabel_1{k}(end-6:end))>0
                endField.(sessDirs{j}) = [endField.(sessDirs{j}),1];
            else
                endField.(sessDirs{j}) = [endField.(sessDirs{j}),0];
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

delay_onIdx_1 = sum(rate_Delay(1,:)>p.avgRateThres,1)>0;
delay_onIdx_2 = sum(rate_Delay(2,:)>p.avgRateThres,1)>0;
delay_onIdx_3 = sum(rate_Delay(3,:)>p.avgRateThres,1)>0;
delay_onIdx_4 = sum(rate_Delay(4,:)>p.avgRateThres,1)>0;
rate_Delay_mean = (rate_Delay(1,:)*10 + rate_Delay(2,:)*10 + rate_Delay(3,:)*30 + rate_Delay(4,:)*30)/80;

% plot wenn of the delay on cell distribution
setListData = {find(delay_onIdx_1==1); find(delay_onIdx_2==1); find(delay_onIdx_3==1); find(delay_onIdx_4==1)};
setLabels = ["on 10"; "off 10"; "on 30"; "off 30"];
% h = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);


%% on10, off10, on30, off30
% all time cells
h2 = figure;
h2.Position = [100 100 1200 400];
% on10 on30
on10Idx_Time = (fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
[~,peakIdx] = max(cellMapTemp1,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = cellMapTemp1(peak_Sort,:);
% on10_mean = mean(cellMapTemp1,2);

% on30 and off30
on30Idx_Time = (fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx);
cellMapTemp2 = timeMap_Def1.on30(on10Idx_Time,:);
cellMapTemp2_Sort = cellMapTemp2(peak_Sort,:);
% on30_mean = mean(cellMapTemp1,2);

plotTimeMap2([pre_Map_Def1.on10(on10Idx_Time,:),cellMapTemp1],[pre_Map_Def1.on30(on10Idx_Time,:),cellMapTemp2],...
    2,[2,2,1],[2,2,2],{'Delay start at barrier','Stable cells on10'},{'-3','0','10','20'},{'on30 by on10'},{'-3','0','10','20'});

on10_sec_1 = mean(cellMapTemp1_Sort(:,1:33),2);
on30_sec_1 = mean(cellMapTemp2_Sort(:,1:33),2);
on10_sec_2 = mean(cellMapTemp1_Sort(:,34:66),2);
on30_sec_2 = mean(cellMapTemp2_Sort(:,34:66),2);
on30_sec_3 = mean(cellMapTemp2_Sort(:,67:99),2);
on30_sec_4 = mean(cellMapTemp2_Sort(:,100:132),2);
on30_sec_5 = mean(cellMapTemp2_Sort(:,133:165),2);
on30_sec_6 = mean(cellMapTemp2_Sort(:,166:198),2);
meanRate_on = (sum(cellMapTemp1_Sort,2) + sum(cellMapTemp2_Sort,2))./(size(cellMapTemp1_Sort,2)+size(cellMapTemp2_Sort,2)); 
% meanRate_on = rate_Delay_mean(on10Idx_Time);
% meanRate_on = meanRate_on(peak_Sort)';
% meanRate_on = mean(cellMapTemp2_Sort,2);

subplot(2,2,3)
yMat = length(on10_sec_1):-1:1;
plot(on10_sec_1-30,yMat,'ko-','MarkerSize',3);
hold on
plot(on30_sec_1-15,yMat,'ko-','MarkerSize',3);

plot(on10_sec_2,yMat,'o-','MarkerSize',3);
plot(on30_sec_2+15,yMat,'o-','MarkerSize',3);
plot(on30_sec_3+30,yMat,'o-','MarkerSize',3);
plot(on30_sec_4+45,yMat,'o-','MarkerSize',3);
plot(on30_sec_5+60,yMat,'o-','MarkerSize',3);
plot(on30_sec_6+75,yMat,'o-','MarkerSize',3);
xlim([-35 95])
set(gca, 'xtick', [-30,-15,0,15,30,45,60,75]);
set(gca, 'xticklabels', {'10s 0-5','30s 0-5','10 s 5-10','30 s 5-10','10-15','15-20','20-25','25-30'},'fontSize',6);
ylim([1 length(on10_sec_1)])

% calculate correlation
[corr_1,pva_1] = corr(on10_sec_2,on30_sec_2,'Type','Pearson');
[corr_2,pva_2] = corr(on30_sec_2,on30_sec_3,'Type','Pearson');
[corr_3,pva_3] = corr(on30_sec_2,on30_sec_4,'Type','Pearson');
[corr_4,pva_4]= corr(on30_sec_2,on30_sec_5,'Type','Pearson');
[corr_5,pva_5] = corr(on30_sec_2,on30_sec_6,'Type','Pearson');
[corr_6,pva_6] = corr(on10_sec_1,on30_sec_1,'Type','Pearson');
[corr_7,pva_7] = corr(on10_sec_1,on10_sec_2,'Type','Pearson');
[corr_8,pva_8] = corr(on30_sec_1,on30_sec_2,'Type','Pearson');

text(-33,length(on10_sec_1),sprintf('%1.2f',corr_6),'fontSize',8);
text(-12,length(on10_sec_1),sprintf('%1.2f',corr_7),'fontSize',8);
text(-5,length(on10_sec_1),sprintf('%1.2f',corr_8),'fontSize',8);

text(5,length(on10_sec_1),sprintf('%1.2f',corr_1),'fontSize',8);
text(20,length(on10_sec_1),sprintf('%1.2f',corr_2),'fontSize',8);
text(35,length(on10_sec_1),sprintf('%1.2f',corr_3),'fontSize',8);
text(50,length(on10_sec_1),sprintf('%1.2f',corr_4),'fontSize',8);
text(65,length(on10_sec_1),sprintf('%1.2f',corr_5),'fontSize',8);


subplot(2,2,4)
yMat = length(on10_sec_1):-1:1;
plot(on10_sec_1./meanRate_on-10,yMat,'ko-','MarkerSize',3);
hold on
plot(on30_sec_1./meanRate_on-5,yMat,'ko-','MarkerSize',3);

plot(on10_sec_2./meanRate_on,yMat,'o-','MarkerSize',3);
plot(on30_sec_2./meanRate_on+2,yMat,'o-','MarkerSize',3);
plot(on30_sec_3./meanRate_on+4,yMat,'o-','MarkerSize',3);
plot(on30_sec_4./meanRate_on+6,yMat,'o-','MarkerSize',3);
plot(on30_sec_5./meanRate_on+8,yMat,'o-','MarkerSize',3);
plot(on30_sec_6./meanRate_on+10,yMat,'o-','MarkerSize',3);
xlim([-13 13])
% set(gca, 'xtick', [-6,-15,0,15,30,45,60,75]);
% set(gca, 'xticklabels', {'10s 0-5','30s 0-5','10 s 5-10','30 s 5-10','10-15','15-20','20-25','25-30'},'fontSize',6);
ylim([1 length(on10_sec_1)])

% calculate correlation
[corr_1,pva_1] = corr(on10_sec_2./meanRate_on,on30_sec_2./meanRate_on,'Type','Pearson');
[corr_2,pva_2] = corr(on30_sec_2./meanRate_on,on30_sec_3./meanRate_on,'Type','Pearson');
[corr_3,pva_3] = corr(on30_sec_2./meanRate_on,on30_sec_4./meanRate_on,'Type','Pearson');
[corr_4,pva_4] = corr(on30_sec_2./meanRate_on,on30_sec_5./meanRate_on,'Type','Pearson');
[corr_5,pva_5] = corr(on30_sec_2./meanRate_on,on30_sec_6./meanRate_on,'Type','Pearson');
[corr_6,pva_6] = corr(on10_sec_1./meanRate_on,on30_sec_1./meanRate_on,'Type','Pearson');
[corr_7,pva_7] = corr(on10_sec_1./meanRate_on,on10_sec_2./meanRate_on,'Type','Pearson');
[corr_8,pva_8] = corr(on30_sec_1./meanRate_on,on30_sec_2./meanRate_on,'Type','Pearson');

text(-8,length(on10_sec_1),sprintf('%1.2f',corr_6),'fontSize',8);
text(-4,length(on10_sec_1),sprintf('%1.2f',corr_7),'fontSize',8);
text(-2,length(on10_sec_1),sprintf('%1.2f',corr_8),'fontSize',8);

text(1,length(on10_sec_1),sprintf('%1.2f',corr_1),'fontSize',8);
text(3,length(on10_sec_1),sprintf('%1.2f',corr_2),'fontSize',8);
text(5,length(on10_sec_1),sprintf('%1.2f',corr_3),'fontSize',8);
text(7,length(on10_sec_1),sprintf('%1.2f',corr_4),'fontSize',8);
text(9,length(on10_sec_1),sprintf('%1.2f',corr_5),'fontSize',8);



off10Idx_Time = (fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.off10(off10Idx_Time,:);
[~,peakIdx] = max(cellMapTemp1,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = cellMapTemp1(peak_Sort,:);
% off10_mean = mean(cellMapTemp2,2);

off30Idx_Time = (fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx);
cellMapTemp2 = timeMap_Def1.off30(off10Idx_Time,:);
cellMapTemp2_Sort = cellMapTemp2(peak_Sort,:);
% off30_mean = mean(cellMapTemp2,2);

figure
plotTimeMap2([pre_Map_Def1.off10(off10Idx_Time,:),cellMapTemp1],[pre_Map_Def1.off30(off10Idx_Time,:),cellMapTemp2],...
    2,[2,2,1],[2,2,2],{'Delay start at barrier','Stable cells off10'},{'-3','0','10','20'},{'off30 by off10'},{'-3','0','10','20'});

off10_sec_1 = mean(cellMapTemp1_Sort(:,1:33),2);
off30_sec_1 = mean(cellMapTemp2_Sort(:,1:33),2);
off10_sec_2 = mean(cellMapTemp1_Sort(:,34:66),2);
off30_sec_2= mean(cellMapTemp2_Sort(:,34:66),2);
off30_sec_3= mean(cellMapTemp2_Sort(:,67:99),2);
off30_sec_4= mean(cellMapTemp2_Sort(:,100:132),2);
off30_sec_5= mean(cellMapTemp2_Sort(:,133:165),2);
off30_sec_6= mean(cellMapTemp2_Sort(:,166:198),2);
meanRate_off = (sum(cellMapTemp1_Sort,2) + sum(cellMapTemp2_Sort,2))./(size(cellMapTemp1_Sort,2)+size(cellMapTemp2_Sort,2));
% meanRate_off = rate_Delay_mean(off10Idx_Time);
% meanRate_off = meanRate_off(peak_Sort);
% meanRate_off = mean(cellMapTemp2_Sort,2);

subplot(2,2,3)
yMat = length(off10_sec_1):-1:1;
plot(off10_sec_1-30,yMat,'ko-','MarkerSize',3);
hold on
plot(off30_sec_1-15,yMat,'ko-','MarkerSize',3);

plot(off10_sec_2,yMat,'o-','MarkerSize',3);
plot(off30_sec_2+15,yMat,'o-','MarkerSize',3);
plot(off30_sec_3+30,yMat,'o-','MarkerSize',3);
plot(off30_sec_4+45,yMat,'o-','MarkerSize',3);
plot(off30_sec_5+60,yMat,'o-','MarkerSize',3);
plot(off30_sec_6+75,yMat,'o-','MarkerSize',3);
xlim([-35 95])
set(gca, 'xtick', [-30,-15,0,15,30,45,60,75]);
set(gca, 'xticklabels', {'10s 0-5','30s 0-5','10 s 5-10','30 s 5-10','10-15','15-20','20-25','25-30'},'fontSize',6);
ylim([1 length(off10_sec_1)])

% calculate correlation
corr_1 = corr(off10_sec_2,off30_sec_2,'Type','Pearson');
corr_2 = corr(off30_sec_2,off30_sec_3,'Type','Pearson');
corr_3 = corr(off30_sec_2,off30_sec_4,'Type','Pearson');
corr_4 = corr(off30_sec_2,off30_sec_5,'Type','Pearson');
corr_5 = corr(off30_sec_2,off30_sec_6,'Type','Pearson');
corr_6 = corr(off10_sec_1,off30_sec_1,'Type','Pearson');
corr_7 = corr(off10_sec_1,off10_sec_2,'Type','Pearson');
corr_8 = corr(off30_sec_1,off30_sec_2,'Type','Pearson');

text(-33,length(off10_sec_1),sprintf('%1.2f',corr_6),'fontSize',8);
text(-12,length(off10_sec_1),sprintf('%1.2f',corr_7),'fontSize',8);
text(-5,length(off10_sec_1),sprintf('%1.2f',corr_8),'fontSize',8);

text(5,length(off10_sec_1),sprintf('%1.2f',corr_1),'fontSize',8);
text(20,length(off10_sec_1),sprintf('%1.2f',corr_2),'fontSize',8);
text(35,length(off10_sec_1),sprintf('%1.2f',corr_3),'fontSize',8);
text(50,length(off10_sec_1),sprintf('%1.2f',corr_4),'fontSize',8);
text(65,length(off10_sec_1),sprintf('%1.2f',corr_5),'fontSize',8);


subplot(2,2,4)
yMat = length(off10_sec_1):-1:1;
plot((off10_sec_1./meanRate_off)-10,yMat,'ko-','MarkerSize',3);
hold on
plot(off30_sec_1./meanRate_off-5,yMat,'ko-','MarkerSize',3);

plot(off10_sec_2./meanRate_off,yMat,'o-','MarkerSize',3);
plot(off30_sec_2./meanRate_off+2,yMat,'o-','MarkerSize',3);
plot(off30_sec_3./meanRate_off+4,yMat,'o-','MarkerSize',3);
plot(off30_sec_4./meanRate_off+6,yMat,'o-','MarkerSize',3);
plot(off30_sec_5./meanRate_off+8,yMat,'o-','MarkerSize',3);
plot(off30_sec_6./meanRate_off+10,yMat,'o-','MarkerSize',3);
xlim([-13 13])
% set(gca, 'xtick', [-6,-15,0,15,30,45,60,75]);
% set(gca, 'xticklabels', {'10s 0-5','30s 0-5','10 s 5-10','30 s 5-10','10-15','15-20','20-25','25-30'},'fontSize',6);
ylim([1 length(off10_sec_1)])

% calculate correlation
[corr_1,pva_1] = corr(off10_sec_2./meanRate_off,off30_sec_2./meanRate_off,'Type','Pearson');
[corr_2,pva_2] = corr(off30_sec_2./meanRate_off,off30_sec_3./meanRate_off,'Type','Pearson');
[corr_3,pva_3] = corr(off30_sec_3./meanRate_off,off30_sec_4./meanRate_off,'Type','Pearson');
[corr_4,pva_4] = corr(off30_sec_4./meanRate_off,off30_sec_5./meanRate_off,'Type','Pearson');
[corr_5,pva_5] = corr(off30_sec_5./meanRate_off,off30_sec_6./meanRate_off,'Type','Pearson');
[corr_6,pva_6] = corr(off10_sec_1./meanRate_off,off30_sec_1./meanRate_off,'Type','Pearson');
[corr_7,pva_7] = corr(off10_sec_1./meanRate_off,off10_sec_2./meanRate_off,'Type','Pearson');
[corr_8,pva_8] = corr(off30_sec_1./meanRate_off,off30_sec_2./meanRate_off,'Type','Pearson');

text(-8,length(off10_sec_1),sprintf('%1.2f',corr_6),'fontSize',8);
text(-4,length(off10_sec_1),sprintf('%1.2f',corr_7),'fontSize',8);
text(-2,length(off10_sec_1),sprintf('%1.2f',corr_8),'fontSize',8);

text(1,length(off10_sec_1),sprintf('%1.2f',corr_1),'fontSize',8);
text(3,length(off10_sec_1),sprintf('%1.2f',corr_2),'fontSize',8);
text(5,length(off10_sec_1),sprintf('%1.2f',corr_3),'fontSize',8);
text(7,length(off10_sec_1),sprintf('%1.2f',corr_4),'fontSize',8);
text(9,length(off10_sec_1),sprintf('%1.2f',corr_5),'fontSize',8);

%%
if p.writeToFile == 1
    save(fullfile(savedir,'timecellQuant_Timebin.mat'), 'timecellQuant_Timebin');
end
clear timecellQuant

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