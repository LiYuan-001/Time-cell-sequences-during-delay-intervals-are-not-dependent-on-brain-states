% it needs to be stable in each session
% and stable in 2 session combined
%
function Fig8TreadmillTimeCell_Cell_Theta_Correlation(inFile,AnalyzeSes)

close all

p.savePlot = 0;
p.writeToFile = 0;
p.avgRateThres = 0.5;

% % Read in input information
sessInfo = SessInfoImport(inFile);

% if p.savePlot
%     % directory for plot figures
%     % generate a folder for each rat eah day under the current folder
%     savedir = sprintf('%s%s',cd,'\Figures\Quantification Figures\Time cell num-dist quant');
%     if ~exist(savedir, 'dir')
%         mkdir(savedir);
%     end
% end
    

for i = AnalyzeSes(1:end)
    mainDir = sessInfo(i).mainDir;
    
    % get theta delta file and channels
    eegCh = str2num(sessInfo(i).avgEEGch{1});    
    tdrFile = fullfile(mainDir,'Cell Property', 'ThetaDeltaRatioAll.mat');
    load(tdrFile);
    eegCh_matrix = nan(length(ThetaDeltaRatioAll.EEGch),1);
    for m = 1:length(eegCh_matrix)
        eegCh_matrix(m) = str2num(ThetaDeltaRatioAll.EEGch{m});
    end
    [~,eegCh_matrix_Ind] = ismember(eegCh,eegCh_matrix);
    
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load time field file
    timeFieldFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayField_TrialbyTrial_2Session.mat');
    load(timeFieldFile);
    % load delayFire map
    delayMapFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap_2Session.mat');
    load(delayMapFile);
    % load delay stability
    delayStabilityFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayFireStability_TrialbyTrial_2Session.mat');
    load(delayStabilityFile);
    % load arm rate file
    armRateFile = fullfile(sessInfo(i).mainDir,'Cell Property', 'Fig8TreadmillArmRate.mat');
    load(armRateFile);
    
    
    clusterNum = length(SpikeProp.AvgRate.Fig8Rate);
    idx = find(SpikeProp.AvgRate.Fig8Rate>0.1 & SpikeProp.AvgRate.Fig8Rate<5);
    
    % get each phase names (no delay etc)
    sessDirs = {'on10','off10','on30','off30'};
    
    rate_Delay = [Fig8TreadmillArmRate.on10.rateDelay(:,1)';Fig8TreadmillArmRate.off10.rateDelay(:,1)';...
        Fig8TreadmillArmRate.on30.rateDelay(:,1)';Fig8TreadmillArmRate.off30.rateDelay(:,1)'];
    delay_onIdx = (sum(rate_Delay>p.avgRateThres,1))>0;
    delay_onIdx = delay_onIdx(idx);
    
    for j = 1:length(sessDirs)
        if i == AnalyzeSes(1)
            timeMap.(sessDirs{j}) = [];
            timeStability.(sessDirs{j}) = [];   
            timeFieldStability.(sessDirs{j}) = []; 
            thetaDeltaMean.(sessDirs{j}) = [];
            peakTime_TimeCell.(sessDirs{j}) = [];
        end
        
        thetaDeltaMean_Session.(sessDirs{j}) = [];
        
        mapStability = DelayFireStability_TrialbyTrial_2Session.(sessDirs{j}).delayCorr_Def1_Trial(idx);
        
        % get rate for every session
        sess1 = strcat(sessDirs{j},'_1');
        sess2 = strcat(sessDirs{j},'_2');
        tdr_meanTemp = [];
        % get delay start for every trial
        % load analyzed positions               
        delayFile = fullfile(sessInfo(i).mainDir,sess1, 'Fig8DelayZonePos.mat');
        load(delayFile);
        
        % get delay start time in each session
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        delayTstart = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
        
        timeStamp = ThetaDeltaRatioAll.(sess1).timeStamp;
        % get from every channel being chosen
        tdrTemp = nanmean(ThetaDeltaRatioAll.(sess1).tdr_Smooth(eegCh_matrix_Ind,:),1);  
        
        
        dt = timeStamp(2)-timeStamp(1);
        
        if contains(sess1,'10')
            maxT = 10;            
        elseif contains(sess1,'30')
            maxT = 28;
        else
            error('Delay time is wrong')
        end
        matLength = round(maxT./dt);
        
        
        [~,startInd] = min(abs(delayTstart - timeStamp));        
        % get for every trial
        for k = 1:length(startInd)
            tdr_meanTemp = [tdr_meanTemp;tdrTemp(startInd(k):startInd(k)+matLength-1)];
        end
        tdrMean.(sessDirs{j}) = tdr_meanTemp;
        
        
        % get the matrix      
        delayFile = fullfile(sessInfo(i).mainDir,sess2, 'Fig8DelayZonePos.mat');
        load(delayFile);
        tdr_meanTemp = [];
        % get delay start time in each session
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        delayTstart = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
        
        timeStamp = ThetaDeltaRatioAll.(sess2).timeStamp;
        % get from every channel being chosen
        tdrTemp = nanmean(ThetaDeltaRatioAll.(sess2).tdr_Smooth(eegCh_matrix_Ind,:),1);              
        
        [~,startInd] = min(abs(delayTstart - timeStamp));        
        % get for every trial
        for k = 1:length(startInd)
            tdr_meanTemp = [tdr_meanTemp;tdrTemp(startInd(k):startInd(k)+matLength-1)];
        end
        tdrMean.(sessDirs{j}) = [tdrMean.(sessDirs{j});tdr_meanTemp];
        
%         tdr_timestamp = timeStamp(startInd(1):startInd(1)+matLength-1) - timeStamp(startInd(1));
%         
        timeMap_Def1.(sessDirs{j}) = [];
        timeMap_Trial_Def1.(sessDirs{j}) = [];
        timeField_Def1.(sessDirs{j}) = [];
        endField.(sessDirs{j}) = [];           
        fieldLabel_Def1.(sessDirs{j}) = DelayField_TrialbyTrial_2Session.(sessDirs{j}).timeField_1(idx);
        
        cellCount = 0;
        for k = idx
            cellCount = cellCount + 1;
            timeMap_Def1.(sessDirs{j}) = [timeMap_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate1_Combined_Smooth{k}];
            timeField_Def1.(sessDirs{j}) = [timeField_Def1.(sessDirs{j});DelayField_TrialbyTrial_2Session.(sessDirs{j}).fieldLabel_1{k}'];
            timeMap_Trial_Def1.(sessDirs{j}){cellCount} = Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate1_Smooth{k};
            
            if sum(DelayField_TrialbyTrial_2Session.(sessDirs{j}).fieldLabel_1{k}(end-6:end))>0
                endField.(sessDirs{j}) = [endField.(sessDirs{j}),1];
            else
                endField.(sessDirs{j}) = [endField.(sessDirs{j}),0];
            end
        end
        
        % find time cell and delay on cell
        % find their matching theta-delta ratio
        Idx_Time = (fieldLabel_Def1.(sessDirs{j})>0 & delay_onIdx & ~endField.(sessDirs{j}));
        timeStability.(sessDirs{j}) = [timeStability.(sessDirs{j}),mapStability(Idx_Time)];
        
        if sum(Idx_Time) > 0
            cellMapTemp1 = timeMap_Def1.(sessDirs{j})(Idx_Time,:);
            fieldTemp1 = timeField_Def1.(sessDirs{j})(Idx_Time,:);
            Ind_Time = find(Idx_Time == 1);
            
            % find each cell peak and their matching theta power
            % dekay firing mat bin width
            dt_Bin = Fig8DelayTimeMap_2Session.timeBin/10^3;
            delay_TimeBin = (1:size(timeMap_Def1.(sessDirs{j}),2))*dt_Bin;
            
            
            
            for m = 1:sum(Idx_Time)
                % take time map
                timeMap.(sessDirs{j}) = [timeMap.(sessDirs{j});cellMapTemp1(m,:)];
                % take field time                
                fieldTime = delay_TimeBin(fieldTemp1(m,:)==1);               
                [~,peakTime] = max(cellMapTemp1(m,:));               
                peakTime = peakTime * dt_Bin;
                peakTime_TimeCell.(sessDirs{j}) = [peakTime_TimeCell.(sessDirs{j}),peakTime];
                
                % if time field >= 0.5 sec, use all field length
                % else use peak time +- 0.25 sec
                if fieldTime(end) - fieldTime(1) >= 0.5
                    thetaIndTemp = round(fieldTime(1)/dt):round(fieldTime(end)/dt);
                else
                    thetaIndTemp = round((peakTime-0.25)/dt):round((peakTime+0.25)/dt);
                end
%                 thetaIndTemp = round((peakTime-0.25)/dt):round((peakTime+0.25)/dt);
                thetaIndTemp = thetaIndTemp(thetaIndTemp>0);
                thetaIndTemp = thetaIndTemp(thetaIndTemp<size(tdr_meanTemp,2));
                
                tdr_Field_Temp = tdr_meanTemp(:,thetaIndTemp);
                % calculate the std/mean of the tdr inside the time field
                tdr_Field_Temp = reshape(tdr_Field_Temp,size(tdr_Field_Temp,1)*size(tdr_Field_Temp,2),1);
                tdr_mean = nanmean(tdr_Field_Temp);
%                 tdr_std = std(tdr_Field_Temp);
%                 tdr_Var = tdr_std/tdr_mean;               
                thetaDeltaMean_Session.(sessDirs{j}) = [thetaDeltaMean_Session.(sessDirs{j}),tdr_mean];
                
                % get field time map corr 
                field_Trial_Map = timeMap_Trial_Def1.(sessDirs{j}){Ind_Time(m)}(:,fieldTemp1(m,:)==1);
                corrTemp = nan(nchoosek(size(field_Trial_Map,1),2),1);
                countTemp = 0;
                for nn = 1:size(field_Trial_Map,1)-1
                    for mm = nn+1:size(field_Trial_Map,1)
                        countTemp = countTemp+1;
                        corrTemp(countTemp) = Xcorrelate(field_Trial_Map(nn,:),field_Trial_Map(mm,:));
                    end
                end
                timeFieldStability.(sessDirs{j}) = [timeFieldStability.(sessDirs{j}),median(corrTemp,'omitnan')];
            end
        end
    end
    tdrMean_AllDelay = nanmean([tdrMean.on10(:);tdrMean.off10(:);tdrMean.on30(:);tdrMean.off30(:)]);
    thetaDeltaMean.on10 = [thetaDeltaMean.on10,thetaDeltaMean_Session.on10/tdrMean_AllDelay];
    thetaDeltaMean.off10 = [thetaDeltaMean.off10,thetaDeltaMean_Session.off10/tdrMean_AllDelay];
    thetaDeltaMean.on30 = [thetaDeltaMean.on30,thetaDeltaMean_Session.on30/tdrMean_AllDelay];
    thetaDeltaMean.off30 = [thetaDeltaMean.off30,thetaDeltaMean_Session.off30/tdrMean_AllDelay];
end

% figure
% plot(timeStability.on10,thetaDeltaMean.on10,'r^');
% hold on
% plot(timeStability.off10,thetaDeltaMean.off10,'k^');
% xlabel('Time cell stability')
% ylabel('Theta delta norm')
% title('10 SEC DELAY')
% 
% figure
% plot(timeStability.on30,thetaDeltaMean.on30,'r^');
% hold on
% plot(timeStability.off30,thetaDeltaMean.off30,'k^');
% xlabel('Time cell stability')
% ylabel('Theta delta norm')
% title('30 SEC DELAY')
%%
figure
subplot(3,2,1)
for k = 1:length(timeFieldStability.on10)
    if peakTime_TimeCell.on10(k) <= 0.5
        plotColor = [0.5,0,0];
    else
        plotColor = [1,0,0];
    end
%     plot(timeFieldStability.on10(k),thetaDeltaMean.on10(k),'o','MarkerEdgeColor',plotColor,'MarkerFaceColor',plotColor);
    plot(timeFieldStability.on10(k),thetaDeltaMean.on10(k),'ro')
hold on
end
for k = 1:length(timeFieldStability.on30)
    if peakTime_TimeCell.on30(k) <= 0.5
        plotColor = [0.5,0,0];
    else
        plotColor = [1,0,0];
    end
%     plot(timeFieldStability.on30(k),thetaDeltaMean.on30(k),'o','MarkerEdgeColor',plotColor,'MarkerFaceColor',plotColor);
    plot(timeFieldStability.on30(k),thetaDeltaMean.on30(k),'ro')
hold on
end
for k = 1:length(timeFieldStability.off10)
    if peakTime_TimeCell.off10(k) <= 0.5
        plotColor = [0,0,0];
    else
        plotColor = [0.5,0.5,0.5];
    end
%     plot(timeFieldStability.off10(k),thetaDeltaMean.off10(k),'o','MarkerEdgeColor',plotColor,'MarkerFaceColor',plotColor);
    plot(timeFieldStability.off10(k),thetaDeltaMean.off10(k),'ko')
hold on
end
for k = 1:length(timeFieldStability.off30)
    if peakTime_TimeCell.off30(k) <= 0.5
        plotColor = [0,0,0];
    else
        plotColor = [0.5,0.5,0.5];
    end
%     plot(timeFieldStability.off30(k),thetaDeltaMean.off30(k),'o','MarkerEdgeColor',plotColor,'MarkerFaceColor',plotColor);
    plot(timeFieldStability.off30(k),thetaDeltaMean.off30(k),'ko')
    hold on
end
% fit on session
x = [timeFieldStability.on10,timeFieldStability.on30];
y1 = [thetaDeltaMean.on10,thetaDeltaMean.on30];
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'r-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
[Rho,P] = corr(x',y1','Type','Spearman');
pText = sprintf('%s%1.2f%s%1.3f','r = ',Rho,' pVal = ',P);
text(0.4,1.5,pText,'Color','red');
% fit off session
x = [timeFieldStability.off10,timeFieldStability.off30];
y1 = [thetaDeltaMean.off10,thetaDeltaMean.off30];
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'k-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
[Rho,P] = corr(x',y1','Type','Spearman');
pText = sprintf('%s%1.2f%s%1.3f','r = ',Rho,' pVal = ',P);
text(0.4,1.6,pText,'Color',[0,0,0]);
xlabel('Time field stability')
ylabel('Theta delta norm')
ylim([0.8 1.6])
    
subplot(3,2,3)
plot([timeFieldStability.on10,timeFieldStability.on30],log2([peakTime_TimeCell.on10,peakTime_TimeCell.on30]),'ro');
hold on
plot([timeFieldStability.off10,timeFieldStability.off30],log2([peakTime_TimeCell.off10,peakTime_TimeCell.off30]),'ko');
xlabel('Time field stability')
ylabel('Peak time (log2)')
% fit on session
x = [timeFieldStability.on10,timeFieldStability.on30];
y1 = log2([peakTime_TimeCell.on10,peakTime_TimeCell.on30]);
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'r-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
[Rho,P] = corr(x',y1','Type','Spearman');
pText = sprintf('%s%1.2f%s%1.3f','r = ',Rho,' pVal = ',P);
text(0.4,4,pText,'Color','red');
% fit off session
x = [timeFieldStability.off10,timeFieldStability.off30];
y1 = log2([peakTime_TimeCell.off10,peakTime_TimeCell.off30]);
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'k-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
[Rho,P] = corr(x',y1','Type','Spearman');
pText = sprintf('%s%1.2f%s%1.3f','r = ',Rho,' pVal = ',P);
text(0.4,5,pText,'Color',[0,0,0]);



subplot(3,2,2)
edges=0.4:0.05:2.0;
histplotCount = histcounts([thetaDeltaMean.on10,thetaDeltaMean.on30],edges,'Normalization','cumcount');
normCount = histplotCount/length([thetaDeltaMean.on10,thetaDeltaMean.on30]);
plot(edges(1:end-1),normCount,'r','LineWidth',2)
hold on
histplotCount = histcounts([thetaDeltaMean.off10,thetaDeltaMean.off30],edges,'Normalization','cumcount');
normCount = histplotCount/length([thetaDeltaMean.off10,thetaDeltaMean.off30]);
plot(edges(1:end-1),normCount,'k','LineWidth',2)
xlabel('Theta delta norm')
ylabel('Normalized count')


subplot(3,2,4)
edges=-0.4:0.05:1;
histplotCount = histcounts([timeFieldStability.on10,timeFieldStability.on30],edges,'Normalization','cumcount');
normCount = histplotCount/length([timeFieldStability.on10,timeFieldStability.on30]);
plot(edges(1:end-1),normCount,'r','LineWidth',2)
hold on
histplotCount = histcounts([timeFieldStability.off10,timeFieldStability.off30],edges,'Normalization','cumcount');
normCount = histplotCount/length([timeFieldStability.off10,timeFieldStability.off30]);
plot(edges(1:end-1),normCount,'k','LineWidth',2)
xlabel('Time field stability')
ylabel('Normalized count')

subplot(3,2,5)
plot(log2([peakTime_TimeCell.on10,peakTime_TimeCell.on30]),[thetaDeltaMean.on10,thetaDeltaMean.on30],'ro');
hold on
plot(log2([peakTime_TimeCell.off10,peakTime_TimeCell.off30]),[thetaDeltaMean.off10,thetaDeltaMean.off30],'ko');
xlabel('Peak time (log2)')
ylabel('Theta delta norm')

% fit on session
x = log2([peakTime_TimeCell.on10,peakTime_TimeCell.on30]);
y1 = [thetaDeltaMean.on10,thetaDeltaMean.on30];
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'r-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
[Rho,P] = corr(x',y1','Type','Spearman');
pText = sprintf('%s%1.2f%s%1.3f','r = ',Rho,' pVal = ',P);
text(0.4,1.5,pText,'Color','red');
% fit off session
x = log2([peakTime_TimeCell.off10,peakTime_TimeCell.off30]);
y1 = [thetaDeltaMean.off10,thetaDeltaMean.off30];
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'k-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
[Rho,P] = corr(x',y1','Type','Spearman');
pText = sprintf('%s%1.2f%s%1.3f','r = ',Rho,' pVal = ',P);
text(0.4,2,pText,'Color',[0,0,0]);

subplot(3,2,6)
edges=0:0.1:10;
histplotCount = histcounts([peakTime_TimeCell.on10,peakTime_TimeCell.on30],edges,'Normalization','cumcount');
normCount = histplotCount/length([peakTime_TimeCell.on10,peakTime_TimeCell.on30]);
plot(edges(1:end-1),normCount,'r','LineWidth',2)
hold on
histplotCount = histcounts([peakTime_TimeCell.off10,peakTime_TimeCell.off30],edges,'Normalization','cumcount');
normCount = histplotCount/length([peakTime_TimeCell.off10,peakTime_TimeCell.off30]);
plot(edges(1:end-1),normCount,'k','LineWidth',2)
xlabel('Peak time')
ylabel('Normalized count')


%%    
figure
Violin([timeFieldStability.on10,timeFieldStability.on30],1,'ViolinColor',[1,0,0],'ShowData',false)
Violin([timeFieldStability.off10,timeFieldStability.off30],2,'ViolinColor',[0,0,0],'ShowData',false)


figure
plot(timeFieldStability.off10,thetaDeltaMean.off10,'k^');
xlabel('Time field stability')
ylabel('Theta delta norm')
title('10 SEC DELAY')
% fit on session
x = timeFieldStability.on10;
y1 = thetaDeltaMean.on10;
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'r-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
pText = sprintf('%s%1.2f%s%1.3f','r = ',r_Val,' pVal = ',p_val);
text(0.5,1.3,pText,'Color','red');
% fit off session
x = timeFieldStability.off10;
y1 = thetaDeltaMean.off10;
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'k-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
pText = sprintf('%s%1.2f%s%1.3f','r = ',r_Val,' pVal = ',p_val);
text(0.5,1,pText,'Color',[0,0,0]);




figure
plot(timeFieldStability.on30,thetaDeltaMean.on30,'r^');
hold on
plot(timeFieldStability.off30,thetaDeltaMean.off30,'k^');
xlabel('Time field stability')
ylabel('Theta delta norm')
title('30 SEC DELAY')
% fit on session
x = timeFieldStability.on30;
y1 = thetaDeltaMean.on30;
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'r-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
pText = sprintf('%s%1.2f%s%1.3f','r = ',r_Val,' pVal = ',p_val);
text(0.5,1.3,pText,'Color','red');
% fit off session
x = timeFieldStability.off30;
y1 = thetaDeltaMean.off30;
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'k-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
pText = sprintf('%s%1.2f%s%1.3f','r = ',r_Val,' pVal = ',p_val);
text(0.5,1,pText,'Color',[0,0,0]);



figure
plot([timeFieldStability.on10,timeFieldStability.on30],[thetaDeltaMean.on10,thetaDeltaMean.on30],'r^');
hold on
plot([timeFieldStability.off10,timeFieldStability.off30],[thetaDeltaMean.off10,thetaDeltaMean.off30],'k^');
xlabel('Time field stability')
ylabel('Theta delta norm')
title('10 SEC DELAY')
% fit on session
x = [timeFieldStability.on10,timeFieldStability.on30];
y1 = [thetaDeltaMean.on10,thetaDeltaMean.on30];
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'r-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
pText = sprintf('%s%1.2f%s%1.3f','r = ',r_Val,' pVal = ',p_val);
text(0.4,1.3,pText,'Color','red');
% fit off session
x = [timeFieldStability.off10,timeFieldStability.off30];
y1 = [thetaDeltaMean.off10,thetaDeltaMean.off30];
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'k-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
pText = sprintf('%s%1.2f%s%1.3f','r = ',r_Val,' pVal = ',p_val);
text(0.4,0.9,pText,'Color',[0,0,0]);




figure
[~,peakIdx] = max(timeMap.on10,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = timeMap.on10(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);
subplot(2,2,1)
imagesc(cellMapTemp1_Sort_Norm)
colormap(jet)
title('on10 time cell','Interpreter','None')
axis on
set(gca, 'xtick', [0 size(cellMapTemp1_Sort_Norm,2)]);
set(gca, 'xticklabels', [0 10]);


subplot(2,2,2)
errorbar(1:length(peakIdx),thetaDeltaMean_Session.on10(peak_Sort),thetaDeltaPowVar.on10(peak_Sort),'-s','Color',[1 0 0])
hold on
errorbar(1:length(peakIdx),thetaDeltaMean_WholeDelay.on10(peak_Sort),thetaDeltaVar_WholeDelay.on10(peak_Sort),'o','Color',[0.5 0.5 0.5])     
view([90 -90])
set(gca,'XDir','Reverse')
title('Theta-Delta ratio')
xlabel('Cell')
ylabel('Theta-Delta ratio')
ylim([0 6])


[~,peakIdx] = max(timeMap.off10,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = timeMap.off10(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);
subplot(2,2,3)
imagesc(cellMapTemp1_Sort_Norm)
colormap(jet)
title('off10 time cell','Interpreter','None')
axis on
set(gca, 'xtick', [0 size(cellMapTemp1_Sort_Norm,2)]);
set(gca, 'xticklabels', [0 10]);


subplot(2,2,4)
errorbar(1:length(peakIdx),thetaDeltaMean_Session.off10(peak_Sort),thetaDeltaPowVar.off10(peak_Sort),'-s','Color',[1 0 0])
hold on
errorbar(1:length(peakIdx),thetaDeltaMean_WholeDelay.off10(peak_Sort),thetaDeltaVar_WholeDelay.off10(peak_Sort),'o','Color',[0.5 0.5 0.5])    
ylim([0 6])
view([90 -90])
set(gca,'XDir','Reverse')
title('Theta-Delta ratio')
xlabel('Cell')
ylabel('Theta-Delta ratio')


figure
[~,peakIdx] = max(timeMap.on30,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = timeMap.on30(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);
subplot(2,2,1)
imagesc(cellMapTemp1_Sort_Norm)
colormap(jet)
title('on30 time cell','Interpreter','None')
axis on
set(gca, 'xtick', [0 size(cellMapTemp1_Sort_Norm,2)]);
set(gca, 'xticklabels', [0 30]);

subplot(2,2,2)
errorbar(1:length(peakIdx),thetaDeltaMean_Session.on30(peak_Sort),thetaDeltaPowVar.on30(peak_Sort),'-s','Color',[1 0 0])
hold on
errorbar(1:length(peakIdx),thetaDeltaMean_WholeDelay.on30(peak_Sort),thetaDeltaVar_WholeDelay.on30(peak_Sort),'o','Color',[0.5 0.5 0.5])    
ylim([0 6])

view([90 -90])
set(gca,'XDir','Reverse')
title('Theta-Delta ratio')
xlabel('Cell')
ylabel('Theta-Delta ratio')


[~,peakIdx] = max(timeMap.off30,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = timeMap.off30(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);
subplot(2,2,3)
imagesc(cellMapTemp1_Sort_Norm)
colormap(jet)
title('off30 time cell','Interpreter','None')
axis on
set(gca, 'xtick', [0 size(cellMapTemp1_Sort_Norm,2)]);
set(gca, 'xticklabels', [0 10]);

subplot(2,2,4)
errorbar(1:length(peakIdx),thetaDeltaMean_Session.off30(peak_Sort),thetaDeltaPowVar.off30(peak_Sort),'-s','Color',[1 0 0])
hold on
errorbar(1:length(peakIdx),thetaDeltaMean_WholeDelay.off30(peak_Sort),thetaDeltaVar_WholeDelay.off30(peak_Sort),'o','Color',[0.5 0.5 0.5])    
ylim([0 6])

view([90 -90])
set(gca,'XDir','Reverse')
title('Theta-Delta ratio')
xlabel('Cell')
ylabel('Theta-Delta ratio')

figure
Violin(thetaDeltaMean_Session.on10,1,'ShowData',false);
Violin(thetaDeltaMean_Session.off10,2,'ShowData',false);
Violin(thetaDeltaMean_Session.on30,3,'ShowData',false);
Violin(thetaDeltaMean_Session.off30,4,'ShowData',false);
%% on10, off10, on30, off30




if p.savePlot == 1
    figName = sprintf('%s%s%s',savedir,'\Onset-Barrier-Time cell in its own session-2Session combined');
    print(figName,'-dpng','-r600');
end



%
if p.writeToFile == 1
    save(fullfile(savedir,'timecellQuant_Timebin.mat'), 'timecellQuant_Timebin');
end
clear timecellQuant

end