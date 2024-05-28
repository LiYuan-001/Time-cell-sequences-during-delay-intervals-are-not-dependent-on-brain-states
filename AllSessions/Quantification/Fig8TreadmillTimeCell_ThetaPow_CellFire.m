% it needs to be stable in each session
% and stable in 2 session combined
%
function Fig8TreadmillTimeCell_ThetaPow_CellFire(inFile,AnalyzeSes)

close all

p.savePlot = 0;
p.writeToFile = 0;
p.avgRateThres = 0.5;

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
    delayMapFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap_2Session.mat');
    load(delayMapFile);
    % load arm rate file
    armRateFile = fullfile(sessInfo(i).mainDir,'Cell Property', 'Fig8TreadmillArmRate.mat');
    load(armRateFile);
    % load theta delta ratio
    thetaDeltaRatioFile = fullfile(sessInfo(i).mainDir,'Cell Property', 'ThetaDeltaRatio.mat');
    load(thetaDeltaRatioFile);
    
    
    clusterNum = length(SpikeProp.AvgRate.Fig8Rate);
    idx = find(SpikeProp.AvgRate.Fig8Rate>0.1 & SpikeProp.AvgRate.Fig8Rate<5);
    
    % get each phase names (no delay etc)
    sessDirs = {'on10','off10','on30','off30'};
    fieldLabelTemp = zeros(1,length(idx));
    
  
    rate_DelayTemp = [Fig8TreadmillArmRate.on10.rateDelay(idx,1)';Fig8TreadmillArmRate.off10.rateDelay(idx,1)';...
        Fig8TreadmillArmRate.on30.rateDelay(idx,1)';Fig8TreadmillArmRate.off30.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];
    delay_onIdxTemp = (sum(rate_DelayTemp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
    
    for j = 1:length(sessDirs)
        if i == AnalyzeSes(1)
            timeMap.(sessDirs{j}) = [];
            thetaPowMean.(sessDirs{j}) = [];
            thetaPowVar.(sessDirs{j}) = [];
            thetaDeltaMean.(sessDirs{j}) = [];
            thetaDeltaPowVar.(sessDirs{j}) = [];
            
            thetaPowMean_WholeDelay.(sessDirs{j}) = [];
            thetaPowVar_WholeDelay.(sessDirs{j}) = [];            
            thetaDeltaMean_WholeDelay.(sessDirs{j}) = [];
            thetaDeltaVar_WholeDelay.(sessDirs{j}) = [];
        end
        
        
        % get rate for every session
        sess1 = strcat(sessDirs{j},'_1');
        sess2 = strcat(sessDirs{j},'_2');
        
        % get delay start for every trial
        % load analyzed positions               
        delayFile = fullfile(sessInfo(i).mainDir,sess1, 'Fig8DelayZonePos.mat');
        load(delayFile);
        
        % get delay start time in each session
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        delayTstart = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
        dt = ThetaDeltaRatio.(sess1).timeStamp(2)-ThetaDeltaRatio.(sess1).timeStamp(1);
        if contains(sess1,'10')
            maxT = 10;
            matLength = round(10./dt);
        elseif contains(sess1,'30')
            maxT = 25;
            matLength = round(25./dt);
        else
            error('Delay time is wrong')
        end
        
        thetaPowMat_1 = zeros(length(delayTstart),matLength);
        theta_deltaMat_1 = zeros(length(delayTstart),matLength);
        
        tsTemp = ThetaDeltaRatio.(sess1).timeStamp;
        for m = 1:length(delayTstart)
            [min_TsDiff,tsStartInd] = min(abs(delayTstart(m) - tsTemp));
            if min_TsDiff > dt
                error('LFP assignment is wrong')
            end
            thetaPowMat_1(m,:) = ThetaDeltaRatio.(sess1).thetaPow(tsStartInd:tsStartInd+matLength-1);
            theta_deltaMat_1(m,:) = ThetaDeltaRatio.(sess1).thetaDeltaRatio(tsStartInd:tsStartInd+matLength-1);
        end

        
        
        % get the matrix
        
        delayFile = fullfile(sessInfo(i).mainDir,sess2, 'Fig8DelayZonePos.mat');
        load(delayFile);
        
        % get delay start time in each session
        % def1: delay starts at entrance
        % def2: delay starts at barrier
        delayTstart = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
        thetaPowMat_2 = zeros(length(delayTstart),matLength);
        theta_deltaMat_2 = zeros(length(delayTstart),matLength);
        
        tsTemp = ThetaDeltaRatio.(sess2).timeStamp;
        for m = 1:length(delayTstart)
            [min_TsDiff,tsStartInd] = min(abs(delayTstart(m) - tsTemp));
            if min_TsDiff > dt
                error('LFP assignment is wrong')
            end
            thetaPowMat_2(m,:) = ThetaDeltaRatio.(sess2).thetaPow(tsStartInd:tsStartInd+matLength-1);
            theta_deltaMat_2(m,:) = ThetaDeltaRatio.(sess2).thetaDeltaRatio(tsStartInd:tsStartInd+matLength-1);
        end        
        
        thetaPowMat = [thetaPowMat_1;thetaPowMat_2];
        theta_deltaMat = [theta_deltaMat_1;theta_deltaMat_2];
              
        timeMap_Def1.(sessDirs{j}) = [];
        timeField_Def1.(sessDirs{j}) = [];
        endField.(sessDirs{j}) = [];           
        fieldLabel_Def1.(sessDirs{j}) = DelayField_TrialbyTrial_2Session.(sessDirs{j}).timeField_1(idx);
        
        for k = idx
            timeMap_Def1.(sessDirs{j}) = [timeMap_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate1_Combined_Smooth{k}];
            timeField_Def1.(sessDirs{j}) = [timeField_Def1.(sessDirs{j});DelayField_TrialbyTrial_2Session.(sessDirs{j}).fieldLabel_1{k}'];
            if sum(DelayField_TrialbyTrial_2Session.(sessDirs{j}).fieldLabel_1{k}(end-6:end))>0
                endField.(sessDirs{j}) = [endField.(sessDirs{j}),1];
            else
                endField.(sessDirs{j}) = [endField.(sessDirs{j}),0];
            end
        end
        
        % find time cell and delay on cell
        % find their matching theta-delta ratio
        Idx_Time = (fieldLabel_Def1.(sessDirs{j})>0 & delay_onIdxTemp & ~endField.(sessDirs{j}));
        
        if sum(Idx_Time) > 0
            cellMapTemp1 = timeMap_Def1.(sessDirs{j})(Idx_Time,:);
            fieldTemp1 = timeField_Def1.(sessDirs{j})(Idx_Time,:);

            
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
                
                % if time field >= 0.5 sec, use all field length
                % else use peak time +- 0.25 sec
                if fieldTime(end) - fieldTime(1) >= 0.5
                    thetaIndTemp = round(fieldTime(1)/dt):round(fieldTime(end)/dt);
                else
                    thetaIndTemp = round((peakTime-0.25)/dt):round((peakTime+0.25)/dt);
                end
%                 thetaIndTemp = round((peakTime-0.25)/dt):round((peakTime+0.25)/dt);
                thetaIndTemp = thetaIndTemp(thetaIndTemp>0);
                thetaIndTemp = thetaIndTemp(thetaIndTemp<size(thetaPowMat,2));
                
                thetaPowTemp = thetaPowMat(:,thetaIndTemp);
                thetaPowVarTemp = std(nanmean(thetaPowTemp,2))/sqrt(length(nanmean(thetaPowTemp,2))-1);
                thetaPowMean.(sessDirs{j}) = [thetaPowMean.(sessDirs{j}),nanmean(nanmean(thetaPowTemp))];
                thetaPowVar.(sessDirs{j}) = [thetaPowVar.(sessDirs{j}),thetaPowVarTemp];
                
                thetaDeltaTemp = theta_deltaMat(:,thetaIndTemp);
                thetaDeltaVarTemp = std(nanmean(thetaDeltaTemp,2))/sqrt(length(nanmean(thetaPowTemp,2))-1);
                thetaDeltaMean.(sessDirs{j}) = [thetaDeltaMean.(sessDirs{j}),nanmean(nanmean(thetaDeltaTemp))];
                thetaDeltaPowVar.(sessDirs{j}) = [thetaDeltaPowVar.(sessDirs{j}),thetaDeltaVarTemp];
                
                thetaPowMean_WholeDelay.(sessDirs{j}) = [thetaPowMean_WholeDelay.(sessDirs{j}),nanmean(nanmean(thetaPowMat))];
                thetaPowVar_WholeDelay_Temp = std(nanmean(thetaPowMat,2))/sqrt(length(nanmean(thetaPowMat,2))-1);
                thetaPowVar_WholeDelay.(sessDirs{j}) = [thetaPowVar_WholeDelay.(sessDirs{j}),thetaPowVar_WholeDelay_Temp];
                
                thetaDeltaMean_WholeDelay.(sessDirs{j}) = [thetaDeltaMean_WholeDelay.(sessDirs{j}),nanmean(nanmean(theta_deltaMat))];
                thetaDeltaVar_WholeDelay_Temp = std(nanmean(theta_deltaMat,2))/sqrt(length(nanmean(theta_deltaMat,2))-1);
                thetaDeltaVar_WholeDelay.(sessDirs{j}) = [thetaDeltaVar_WholeDelay.(sessDirs{j}),thetaDeltaVar_WholeDelay_Temp];
            end
        end
    end
end

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
errorbar(1:length(peakIdx),thetaDeltaMean.on10(peak_Sort),thetaDeltaPowVar.on10(peak_Sort),'-s','Color',[1 0 0])
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
errorbar(1:length(peakIdx),thetaDeltaMean.off10(peak_Sort),thetaDeltaPowVar.off10(peak_Sort),'-s','Color',[1 0 0])
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
errorbar(1:length(peakIdx),thetaDeltaMean.on30(peak_Sort),thetaDeltaPowVar.on30(peak_Sort),'-s','Color',[1 0 0])
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
errorbar(1:length(peakIdx),thetaDeltaMean.off30(peak_Sort),thetaDeltaPowVar.off30(peak_Sort),'-s','Color',[1 0 0])
hold on
errorbar(1:length(peakIdx),thetaDeltaMean_WholeDelay.off30(peak_Sort),thetaDeltaVar_WholeDelay.off30(peak_Sort),'o','Color',[0.5 0.5 0.5])    
ylim([0 6])

view([90 -90])
set(gca,'XDir','Reverse')
title('Theta-Delta ratio')
xlabel('Cell')
ylabel('Theta-Delta ratio')


figure
subplot(2,2,1)
[~,peakIdx] = max(timeMap.on10,[],2);
plot(peakIdx,thetaDeltaMean.on10,'ro')
hold on
[~,peakIdx] = max(timeMap.off10,[],2);
plot(peakIdx,thetaDeltaMean.off10,'ko')



figure
Violin(thetaDeltaMean.on10,1,'ShowData',false);
Violin(thetaDeltaMean.off10,2,'ShowData',false);
Violin(thetaDeltaMean.on30,3,'ShowData',false);
Violin(thetaDeltaMean.off30,4,'ShowData',false);
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