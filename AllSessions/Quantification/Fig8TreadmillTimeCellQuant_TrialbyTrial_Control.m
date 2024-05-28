% it needs to be stable in each session
% and stable in 2 session combined
%
function Fig8TreadmillTimeCellQuant_TrialbyTrial_Control(inFile,AnalyzeSes)

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
    armRate = load(armRateFile);
    
    
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
            post_Map_Def1.(sessDirs{j}) = [post_Map_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).post_mapMean_Def1_2Session{k}(1:20)];
            
            if sum(DelayField_TrialbyTrial_2Session.(sessDirs{j}).fieldLabel_1{k}(end-6:end))>0
                endField.(sessDirs{j}) = [endField.(sessDirs{j}),1];
            else
                endField.(sessDirs{j}) = [endField.(sessDirs{j}),0];
            end
        end        
    end 
    rate_DelayTemp = [armRate.Fig8TreadmillArmRate.(sessDirs{1}).rateDelay(idx,1)';armRate.Fig8TreadmillArmRate.(sessDirs{2}).rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];
    delay_onIdxTemp = (sum(rate_DelayTemp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
    
end

close all

delay_onIdx_1 = sum(rate_Delay(1,:)>p.avgRateThres,1)>0;
delay_onIdx_2 = sum(rate_Delay(2,:)>p.avgRateThres,1)>0;
% plot wenn of the delay on cell distribution
% setListData = {find(delay_onIdx_1==1); find(delay_onIdx_2==1); find(delay_onIdx_3==1); find(delay_onIdx_4==1)};
% setLabels = ["on 10"; "off 10"; "on 30"; "off 30"];
% h = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);

%% on10, off10, on30, off30

% on10 on30
if contains(sessDirs{1},'on')
    % all delay active cells
    h1 = figure;
    h1.Position = [100 100 600 400];
    on10Idx_Time = (delay_onIdx==1);
    cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
    
    on30Idx_Time = (delay_onIdx==1);
    cellMapTemp2 = timeMap_Def1.on30(on30Idx_Time,:);
    on30_mean = mean(cellMapTemp2,2);

    plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.on30(on30Idx_Time,:),post_Map_Def1.on30(on30Idx_Time,:),...
        1,[1,2,1],[1,2,2],{'Delay start at barrier','Delay on10'},{'-3','0','10','13'},{'Delay on30'},{'-3','0','30','33'});
    % figure
    figure       
    on10_mean = mean(timeMap_Def1.on10(delay_onIdx==1,:),2);
    on30_mean = mean(timeMap_Def1.on30(delay_onIdx==1,:),2);
    plot(on10_mean,on30_mean,'ro','MarkerSize',3)
    xlabel('on10 rate')
    ylabel('on30 rate')
    title('Delay active cell control on group')
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
    [Rho,P] = corr(x,y1,'Type','Spearman');
    pText = sprintf('%s%1.2f%s%1.3f','r = ',Rho,' pVal = ',P);
    text(5,5,pText,'Color','red');
     
    % all time cells
    h1 = figure;
    h1.Position = [100 100 600 400];
%     on10Idx_Time = (fieldLabel_Def1.on10>0 & delay_onIdx & ~endField.on10);
    on10Idx_Time = (fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & endField.on10==0 & delay_onIdx);
    cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
    on10_mean = mean(cellMapTemp1,2);
    
    on30Idx_Time = (fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & endField.on30==0 & delay_onIdx);
    cellMapTemp2 = timeMap_Def1.on30(on30Idx_Time,:);
    on30_mean = mean(cellMapTemp2,2);
    
    plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.on30(on30Idx_Time,:),post_Map_Def1.on30(on30Idx_Time,:),...
        1,[1,2,1],[1,2,2],{'Delay start at barrier','Time cells on10'},{'-3','0','10','13'},{'Time cells on30'},{'-3','0','30','33'});
    
    % on30 sorted by on10
    h1 = figure;
    h1.Position = [100 100 600 400];
    cellMapTemp2 = timeMap_Def1.on30(on10Idx_Time,:);
    plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.on30(on10Idx_Time,:),post_Map_Def1.on30(on10Idx_Time,:),...
        2,[1,2,1],[1,2,2],{'Delay start at barrier','Time cells on10'},{'-3','0','10','13'},{'on30 sort by on10'},{'-3','0','30','33'});
    
    % peak time corr
    figure
    [~,peakIdx] = max(cellMapTemp1,[],2);
    [peakIdxSort,peak_Sort] = sort(peakIdx);
    map1 = cellMapTemp1(peak_Sort,:);
    [~,peakIdx1] = max(map1,[],2);
    plot(log2(peakIdx1*0.15),length(peak_Sort):-1:1,'r^-');
    hold on
    map2 = cellMapTemp2(peak_Sort,:);
    [~,peakIdx2] = max(map2,[],2);
    plot(log2(peakIdx2*0.15),length(peak_Sort):-1:1,'ro--');
    [Rho,P] = corr(peakIdx1,peakIdx2,'Type','Spearman');
    pText = sprintf('%s%1.2f%s%1.3f','r = ',Rho,' pVal = ',P);
    text(1,1,pText,'Color','red');
    xlabel('Time (log2)')
else
    % off10 off30   
    % all delay active cells
    h1 = figure;
    h1.Position = [100 100 600 400];
    off10Idx_Time = (delay_onIdx==1);
    cellMapTemp1 = timeMap_Def1.off10(off10Idx_Time,:);
    off10_mean = mean(cellMapTemp1,2);
    
    off30Idx_Time = (delay_onIdx==1);
    cellMapTemp2 = timeMap_Def1.off30(off30Idx_Time,:);
    off30_mean = mean(cellMapTemp2,2);
    
    plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),pre_Map_Def1.off30(off30Idx_Time,:),post_Map_Def1.off30(off30Idx_Time,:),...
        1,[1,2,1],[1,2,2],{'Delay start at barrier','Delay off10'},{'-3','0','10','13'},{'Delay off30'},{'-3','0','30','33'});
    
    figure
    off10_mean = mean(timeMap_Def1.off10(delay_onIdx==1,:),2);
    off30_mean = mean(timeMap_Def1.off30(delay_onIdx==1,:),2);
    plot(off10_mean,off30_mean,'k^','MarkerSize',3)
    xlabel('off10 rate')
    ylabel('off30 rate')
    title('Delay active cell control off group')
    % fit on session
    x = off10_mean;
    y1 = off30_mean;
    [P,S] = polyfit(x,y1,1);
    yfit = P(1)*x+P(2);
    hold on;
    plot(x,yfit,'k-');
    lm = fitlm(x,y1);
    p_val = lm.coefTest;
    r_Val = lm.Rsquared.Ordinary;
    [Rho,P] = corr(x,y1,'Type','Spearman');
    pText = sprintf('%s%1.2f%s%1.3f','r = ',Rho,' pVal = ',P);
    text(5,5,pText,'Color','red');
    
    % all time cells
    h1 = figure;
    h1.Position = [100 100 600 400];
    off10Idx_Time = (fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & endField.off10==0 & delay_onIdx);
    cellMapTemp1 = timeMap_Def1.off10(off10Idx_Time,:);
    
    off30Idx_Time = (fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & endField.off30==0 & delay_onIdx);
    cellMapTemp2 = timeMap_Def1.off30(off30Idx_Time,:);
    
    plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),pre_Map_Def1.off30(off30Idx_Time,:),post_Map_Def1.off30(off30Idx_Time,:),...
        1,[1,2,1],[1,2,2],{'Delay start at barrier','Time cells off10'},{'-3','0','10','13'},{'Time cells off30'},{'-3','0','30','33'});
    
    % on30 sorted by on10
    h1 = figure;
    h1.Position = [100 100 600 400];
    cellMapTemp2 = timeMap_Def1.off30(off10Idx_Time,:);
    plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),pre_Map_Def1.off30(off10Idx_Time,:),post_Map_Def1.off30(off10Idx_Time,:),...
        2,[1,2,1],[1,2,2],{'Delay start at barrier','Time cells off10'},{'-3','0','10','13'},{'off30 sort by off10'},{'-3','0','30','33'});
    
    % peak time corr
    figure
    [~,peakIdx] = max(cellMapTemp1,[],2);
    [peakIdxSort,peak_Sort] = sort(peakIdx);
    map1 = cellMapTemp1(peak_Sort,:);
    [~,peakIdx1] = max(map1,[],2);
    plot(log2(peakIdx1*0.15),length(peak_Sort):-1:1,'k^-');
    hold on
    map2 = cellMapTemp2(peak_Sort,:);
    [~,peakIdx2] = max(map2,[],2);
    plot(log2(peakIdx2*0.15),length(peak_Sort):-1:1,'ko--');
    [Rho,P] = corr(peakIdx1,peakIdx2,'Type','Spearman');
    pText = sprintf('%s%1.2f%s%1.3f','r = ',Rho,' pVal = ',P);
    text(1,1,pText,'Color','red');
    xlabel('Time (log2)')
    
end
print(figName,'-depsc', '-painters');




figure
x = 0:0.05:10;
mu = mean(on10_mean);
stdTemp = std(on10_mean);
pdf_normal = pdf('Beta',x,mu,stdTemp);
plot(x,pdf_normal,'LineWidth',2)
hold on
mu = mean(off10_mean);
stdTemp = std(off10_mean);
pdf_normal = pdf('Beta',x,mu,stdTemp);
plot(x,pdf_normal,'LineWidth',2)
title('10 sec rate pdf')

figure
mu = mean(on30_mean);
stdTemp = std(on30_mean);
pdf_normal = pdf('Beta',x,mu,stdTemp);
plot(x,pdf_normal,'LineWidth',2)
hold on

mu = mean(off30_mean);
stdTemp = std(off30_mean);
pdf_normal = pdf('Beta',x,mu,stdTemp);
plot(x,pdf_normal,'LineWidth',2)
title('30 sec rate pdf')

figure
boxplot(on10_mean,'Position',1,'Width',0.4)
hold on
boxplot(off10_mean,'Position',2,'Width',0.4)
boxplot(on30_mean,'Position',3,'Width',0.4)
boxplot(off30_mean,'Position',4,'Width',0.4)
set(gca, 'XTick', [1,2,3,4], 'XTickLabel', {'on10','off10','on30','off30'});




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