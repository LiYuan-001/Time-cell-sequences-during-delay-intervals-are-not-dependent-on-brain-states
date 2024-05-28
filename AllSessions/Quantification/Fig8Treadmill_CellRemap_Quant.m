% it needs to be stable in each session
% and stable in 2 session combined
%
function Fig8Treadmill_CellRemap_Quant(inFile,AnalyzeSes)

close all

p.savePlot = 0;
p.writeToFile = 0;
p.avgRateThres = 0.5;

avgRate = [];
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
% plot wenn of the delay on cell distribution
setListData = {find(delay_onIdx_1==1); find(delay_onIdx_2==1); find(delay_onIdx_3==1); find(delay_onIdx_4==1)};
setLabels = ["on 10"; "off 10"; "on 30"; "off 30"];
% h = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);

%% all delay on cells sorted by on10 average rate
h1 = figure;
h1.Position = [100 100 900 800];
% on10 
onIdx = (delay_onIdx==1);
on10_mean = rate_Delay(1,onIdx);
[~,on10_order] = sort(on10_mean);
cellMapTemp1 = timeMap_Def1.on10(onIdx,:);
cellMapTemp1 = cellMapTemp1(on10_order,:);
pre_Map_Temp1 = pre_Map_Def1.on10(onIdx,:);
pre_Map_Temp1 = pre_Map_Temp1(on10_order,:);
% keep 3 sec
post_map_Temp1 = post_Map_Def1.on10(onIdx,1:20);
post_map_Temp1 = post_map_Temp1(on10_order,:);

% on30 
on30_mean = rate_Delay(3,onIdx);
cellMapTemp2 = timeMap_Def1.on30(onIdx,:);
cellMapTemp2 = cellMapTemp2(on10_order,:);
pre_Map_Temp2 = pre_Map_Def1.on30(onIdx,:);
pre_Map_Temp2 = pre_Map_Temp2(on10_order,:);
post_map_Temp2 = post_Map_Def1.on30(onIdx,1:20);
post_map_Temp2 = post_map_Temp2(on10_order,:);

cellMapTemp1_Sort_Norm = [pre_Map_Temp1,cellMapTemp1,post_map_Temp1];
cellMapTemp2_Sort_Norm = [pre_Map_Temp2,cellMapTemp2,post_map_Temp2];

subplot(2,4,1)
imagesc(cellMapTemp1_Sort_Norm)
colormap(jet)
title({'Delay start at barrier','Delay active cells sort by on10'},'Interpreter','None')
axis on
set(gca, 'xtick', [0,size(pre_Map_Temp1,2),size([pre_Map_Temp1,cellMapTemp1],2),size(cellMapTemp1_Sort_Norm,2)]);
set(gca, 'xticklabels', {'-3','0','10','13'});
caxis([0 5])

subplot(2,4,2)
imagesc(cellMapTemp2_Sort_Norm)
colormap(jet)
title({'on30 by on10'},'Interpreter','None')
axis on
set(gca, 'xtick', [0,size(pre_Map_Temp2,2),size([pre_Map_Temp2,cellMapTemp2],2),size(cellMapTemp2_Sort_Norm,2)]);
set(gca, 'xticklabels', {'-3','0','30','33'});
caxis([0 5])


% off10 
off10_mean = rate_Delay(2,onIdx);
cellMapTemp2 = timeMap_Def1.off10(onIdx,:);
cellMapTemp2 = cellMapTemp2(on10_order,:);
pre_Map_Temp2 = pre_Map_Def1.off10(onIdx,:);
pre_Map_Temp2 = pre_Map_Temp2(on10_order,:);
post_map_Temp2 = post_Map_Def1.off10(onIdx,1:20);
post_map_Temp2 = post_map_Temp2(on10_order,:);

cellMapTemp2_Sort_Norm = [pre_Map_Temp2,cellMapTemp2,post_map_Temp2];
subplot(2,4,3)
imagesc(cellMapTemp2_Sort_Norm)
colormap(jet)
title({'off10 by on10'},'Interpreter','None')
axis on
set(gca, 'xtick', [0,size(pre_Map_Temp2,2),size([pre_Map_Temp2,cellMapTemp2],2),size(cellMapTemp2_Sort_Norm,2)]);
set(gca, 'xticklabels', {'-3','0','10','13'});
caxis([0 5])

% off30 
off30_mean = rate_Delay(4,onIdx);
cellMapTemp2 = timeMap_Def1.off30(onIdx,:);
cellMapTemp2 = cellMapTemp2(on10_order,:);
pre_Map_Temp2 = pre_Map_Def1.off30(onIdx,:);
pre_Map_Temp2 = pre_Map_Temp2(on10_order,:);
post_map_Temp2 = post_Map_Def1.off30(onIdx,1:20);
post_map_Temp2 = post_map_Temp2(on10_order,:);

cellMapTemp2_Sort_Norm = [pre_Map_Temp2,cellMapTemp2,post_map_Temp2];
subplot(2,4,4)
imagesc(cellMapTemp2_Sort_Norm)
colormap(jet)
title({'off30 by on10'},'Interpreter','None')
axis on
set(gca, 'xtick', [0,size(pre_Map_Temp2,2),size([pre_Map_Temp2,cellMapTemp2],2),size(cellMapTemp2_Sort_Norm,2)]);
set(gca, 'xticklabels', {'-3','0','13','30'});
caxis([0 5])


% sort by off10
% off 10 
onIdx = (delay_onIdx==1);
[~,off10_order] = sort(off10_mean);
cellMapTemp1 = timeMap_Def1.off10(onIdx,:);
cellMapTemp1 = cellMapTemp1(off10_order,:);
pre_Map_Temp1 = pre_Map_Def1.off10(onIdx,:);
pre_Map_Temp1 = pre_Map_Temp1(off10_order,:);
post_map_Temp1 = post_Map_Def1.off10(onIdx,1:20);
post_map_Temp1 = post_map_Temp1(off10_order,:);

% off 30 
on30_mean = rate_Delay(3,onIdx);
cellMapTemp2 = timeMap_Def1.off30(onIdx,:);
cellMapTemp2 = cellMapTemp2(off10_order,:);
pre_Map_Temp2 = pre_Map_Def1.off30(onIdx,:);
pre_Map_Temp2 = pre_Map_Temp2(off10_order,:);
post_map_Temp2 = post_Map_Def1.off30(onIdx,1:20);
post_map_Temp2 = post_map_Temp2(off10_order,:);

cellMapTemp1_Sort_Norm = [pre_Map_Temp1,cellMapTemp1,post_map_Temp1];
cellMapTemp2_Sort_Norm = [pre_Map_Temp2,cellMapTemp2,post_map_Temp2];

subplot(2,4,5)
imagesc(cellMapTemp1_Sort_Norm)
colormap(jet)
title({'Delay start at barrier','Delay active cells sort by off10'},'Interpreter','None')
axis on
set(gca, 'xtick', [0,size(pre_Map_Temp1,2),size([pre_Map_Temp1,cellMapTemp1],2),size(cellMapTemp1_Sort_Norm,2)]);
set(gca, 'xticklabels', {'-3','0','10','13'});
caxis([0 5])

subplot(2,4,6)
imagesc(cellMapTemp2_Sort_Norm)
colormap(jet)
title({'off30 by off10'},'Interpreter','None')
axis on
set(gca, 'xtick', [0,size(pre_Map_Temp2,2),size([pre_Map_Temp2,cellMapTemp2],2),size(cellMapTemp2_Sort_Norm,2)]);
set(gca, 'xticklabels', {'-3','0','30','33'});
caxis([0 5])

% on10 
cellMapTemp2 = timeMap_Def1.on10(onIdx,:);
cellMapTemp2 = cellMapTemp2(off10_order,:);
pre_Map_Temp2 = pre_Map_Def1.on10(onIdx,:);
pre_Map_Temp2 = pre_Map_Temp2(off10_order,:);
post_map_Temp2 = post_Map_Def1.on10(onIdx,1:20);
post_map_Temp2 = post_map_Temp2(off10_order,:);

cellMapTemp2_Sort_Norm = [pre_Map_Temp2,cellMapTemp2,post_map_Temp2];
subplot(2,4,7)
imagesc(cellMapTemp2_Sort_Norm)
colormap(jet)
title({'on10 by off10'},'Interpreter','None')
axis on
set(gca, 'xtick', [0,size(pre_Map_Temp2,2),size([pre_Map_Temp2,cellMapTemp2],2),size(cellMapTemp2_Sort_Norm,2)]);
set(gca, 'xticklabels', {'-3','0','10','13'});
caxis([0 5])

% on30 
cellMapTemp2 = timeMap_Def1.on30(onIdx,:);
cellMapTemp2 = cellMapTemp2(off10_order,:);
pre_Map_Temp2 = pre_Map_Def1.on30(onIdx,:);
pre_Map_Temp2 = pre_Map_Temp2(off10_order,:);
post_map_Temp2 = post_Map_Def1.on30(onIdx,1:20);
post_map_Temp2 = post_map_Temp2(off10_order,:);

% cellMapTemp2_Sort_Norm = [pre_Map_Temp2,cellMapTemp2,post_map_Temp2]./max([pre_Map_Temp2,cellMapTemp2,post_map_Temp2],[],2);
cellMapTemp2_Sort_Norm = [pre_Map_Temp2,cellMapTemp2,post_map_Temp2];

subplot(2,4,8)
imagesc(cellMapTemp2_Sort_Norm)
colormap(jet)
title({'on30 by off10'},'Interpreter','None')
axis on
set(gca, 'xtick', [0,size(pre_Map_Temp2,2),size([pre_Map_Temp2,cellMapTemp2],2),size(cellMapTemp2_Sort_Norm,2)]);
set(gca, 'xticklabels', {'-3','0','30','33'});
caxis([0 5])

% figure
figure
subplot(3,1,1)
plot(on10_mean,on30_mean,'ro','MarkerSize',3)
hold on
plot(off10_mean,off30_mean,'k^','MarkerSize',3)
xlim([0 15])
ylim([0 15])
title('rate 10 vs 30')
xlabel('10 sec')
ylabel('30 sec')

% fit on session
x = on10_mean;
y1 = on30_mean;
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'r-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
[Rho,P] = corr(x',y1','Type','Spearman');
pText = sprintf('%s%1.2f%s%1.3f','r = ',Rho,' pVal = ',P);
text(5,10,pText,'Color','red');
% fit off session
x = off10_mean;
y1 = off30_mean;
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'k-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
[Rho,P] = corr(x',y1','Type','Spearman');
pText = sprintf('%s%1.2f%s%1.3f','r = ',Rho,' pVal = ',P);
text(5,8,pText,'Color','k');


subplot(3,1,2)
plot(on10_mean,off10_mean,'mo','MarkerSize',3)
hold on
plot(on30_mean,off30_mean,'b^','MarkerSize',3)
xlim([0 15])
ylim([0 15])
title('rate on vs off')
xlabel('on')
ylabel('off')

% fit 10 sec session
x = on10_mean;
y1 = off10_mean;
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'m-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
[Rho,P] = corr(x',y1','Type','Spearman');
pText = sprintf('%s%1.2f%s%1.3f','r = ',Rho,' pVal = ',P);
text(5,10,pText,'Color','red');
% fit 30 sec session
x = on30_mean;
y1 = off30_mean;
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'b-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
[Rho,P] = corr(x',y1','Type','Spearman');
pText = sprintf('%s%1.2f%s%1.3f','r = ',Rho,' pVal = ',P);
text(5,8,pText,'Color','b');

subplot(3,1,3)
temp1 = (on10_mean-on30_mean)./(on10_mean+on30_mean);
temp2 = (off10_mean-off30_mean)./(off10_mean+off30_mean);
temp3 = (on10_mean-off10_mean)./(on10_mean+off10_mean);
temp4 = (on30_mean-off30_mean)./(on30_mean+off30_mean);
Violin(abs([temp1,temp2]),1,'ShowData',false);
Violin(abs([temp3,temp4]),2,'ShowData',false);
title('Diff / Sum')
set(gca,'XTick',[1,2],'XTickLabel',{'10 vs 30','on vs off'});

% fisher's Z
x1 = [on10_mean,off10_mean];
y1 = [on30_mean,off30_mean];
[Rho_1,P_1] = corr(x1',y1','Type','Spearman');

x2 = [on10_mean,on30_mean];
y2 = [off10_mean,off30_mean];
[Rho_2,P_2] = corr(x2',y2','Type','Spearman');

%[p, z, za, zb] = corr_rtest(ra, rb, na, nb)

%
% Ryosuke F Takeuchi 2017/02/02 - 
	za  =  0.5*log((1+Rho_1)/(1-Rho_1));
	zb  = 0.5*log((1+Rho_2)/(1-Rho_2));
	szab = sqrt(1/(length(x1)-3) + 1/(length(x2)-3));
	z = abs(za-zb)/szab;
	pVal(1) = 1-normcdf(z, 0, 1);
	pVal(2) = 2*pVal(1);	

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