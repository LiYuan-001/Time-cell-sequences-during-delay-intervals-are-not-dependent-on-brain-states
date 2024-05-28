% it needs to be stable in each session
% and stable in 2 session combined
%
clear all
inFile = 'Fig8Treadmill_OnOff.xlsx';

close all
p.avgRateThres = 0.5;
% specifiy each rat session
rat1043 = 1:4;
rat1044 = 5:8;
rat1046 = 9:12;
rat1058 = 13:15;
rat1079 = 16:18;
ratAll = 1:18;

% % Read in input information
sessInfo = SessInfoImport(inFile);

% rat 1043
AnalyzeSes = rat1043;    
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
%     % avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
        
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
%% on10, off10, on30, off30
% all time cells
on10Idx_Time = (fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & endField_1.on10==0);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.on10 = sum(peakInd<34)/length(peakInd)*100;

off10Idx_Time = (fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & endField_1.off10==0);
cellMapTemp1 = timeMap_Def1.off10(off10Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.off10 = sum(peakInd<34)/length(peakInd)*100;

on30Idx_Time = (fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & endField_1.on30==0);
cellMapTemp1 = timeMap_Def1.on30(on30Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.on30 = sum(peakInd<34)/length(peakInd)*100;

off30Idx_Time = (fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & endField_1.off30==0);
cellMapTemp1 = timeMap_Def1.off30(off30Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.off30 = sum(peakInd<34)/length(peakInd)*100;
    
peakTime.rat1043.on10_mean = nanmean(within5s.on10);
peakTime.rat1043.off10_mean = nanmean(within5s.off10);
peakTime.rat1043.on30_mean = nanmean(within5s.on30);
peakTime.rat1043.off30_mean = nanmean(within5s.off30);

peakTime.rat1043.on10_std = std(within5s.on10,'omitnan');
peakTime.rat1043.off10_std = std(within5s.off10,'omitnan');
peakTime.rat1043.on30_std = std(within5s.on30,'omitnan');
peakTime.rat1043.off30_std = std(within5s.off30,'omitnan');

% rat 1044
AnalyzeSes = rat1044;   
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
%     % avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
        
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
%% on10, off10, on30, off30
% all time cells
on10Idx_Time = (fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & endField_1.on10==0);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.on10 = sum(peakInd<34)/length(peakInd)*100;

off10Idx_Time = (fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & endField_1.off10==0);
cellMapTemp1 = timeMap_Def1.off10(off10Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.off10 = sum(peakInd<34)/length(peakInd)*100;

on30Idx_Time = (fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & endField_1.on30==0);
cellMapTemp1 = timeMap_Def1.on30(on30Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.on30 = sum(peakInd<34)/length(peakInd)*100;

off30Idx_Time = (fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & endField_1.off30==0);
cellMapTemp1 = timeMap_Def1.off30(off30Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.off30 = sum(peakInd<34)/length(peakInd)*100;


peakTime.rat1044.on10_mean = nanmean(within5s.on10);
peakTime.rat1044.off10_mean = nanmean(within5s.off10);
peakTime.rat1044.on30_mean = nanmean(within5s.on30);
peakTime.rat1044.off30_mean = nanmean(within5s.off30);

peakTime.rat1044.on10_std = std(within5s.on10,'omitnan');
peakTime.rat1044.off10_std = std(within5s.off10,'omitnan');
peakTime.rat1044.on30_std = std(within5s.on30,'omitnan');
peakTime.rat1044.off30_std = std(within5s.off30,'omitnan');


% rat 1046
AnalyzeSes = rat1046;    
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
    % avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
        
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
%% on10, off10, on30, off30
% all time cells
on10Idx_Time = (fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & endField_1.on10==0);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.on10 = sum(peakInd<34)/length(peakInd)*100;

off10Idx_Time = (fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & endField_1.off10==0);
cellMapTemp1 = timeMap_Def1.off10(off10Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.off10 = sum(peakInd<34)/length(peakInd)*100;

on30Idx_Time = (fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & endField_1.on30==0);
cellMapTemp1 = timeMap_Def1.on30(on30Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.on30 = sum(peakInd<34)/length(peakInd)*100;

off30Idx_Time = (fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & endField_1.off30==0);
cellMapTemp1 = timeMap_Def1.off30(off30Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.off30 = sum(peakInd<34)/length(peakInd)*100;


peakTime.rat1046.on10_mean = nanmean(within5s.on10);
peakTime.rat1046.off10_mean = nanmean(within5s.off10);
peakTime.rat1046.on30_mean = nanmean(within5s.on30);
peakTime.rat1046.off30_mean = nanmean(within5s.off30);

peakTime.rat1046.on10_std = std(within5s.on10,'omitnan');
peakTime.rat1046.off10_std = std(within5s.off10,'omitnan');
peakTime.rat1046.on30_std = std(within5s.on30,'omitnan');
peakTime.rat1046.off30_std = std(within5s.off30,'omitnan');


% rat 1058
AnalyzeSes = rat1058;    
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
    % avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
        
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
%% on10, off10, on30, off30
% all time cells
on10Idx_Time = (fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & endField_1.on10==0);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.on10 = sum(peakInd<34)/length(peakInd)*100;

off10Idx_Time = (fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & endField_1.off10==0);
cellMapTemp1 = timeMap_Def1.off10(off10Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.off10 = sum(peakInd<34)/length(peakInd)*100;

on30Idx_Time = (fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & endField_1.on30==0);
cellMapTemp1 = timeMap_Def1.on30(on30Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.on30 = sum(peakInd<34)/length(peakInd)*100;

off30Idx_Time = (fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & endField_1.off30==0);
cellMapTemp1 = timeMap_Def1.off30(off30Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.off30 = sum(peakInd<34)/length(peakInd)*100;


peakTime.rat1058.on10_mean = nanmean(within5s.on10);
peakTime.rat1058.off10_mean = nanmean(within5s.off10);
peakTime.rat1058.on30_mean = nanmean(within5s.on30);
peakTime.rat1058.off30_mean = nanmean(within5s.off30);

peakTime.rat1058.on10_std = std(within5s.on10,'omitnan');
peakTime.rat1058.off10_std = std(within5s.off10,'omitnan');
peakTime.rat1058.on30_std = std(within5s.on30,'omitnan');
peakTime.rat1058.off30_std = std(within5s.off30,'omitnan');


% rat 1079
AnalyzeSes = rat1079;    
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
    % avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
        
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
%% on10, off10, on30, off30
% all time cells
on10Idx_Time = (fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & endField_1.on10==0);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.on10 = sum(peakInd<34)/length(peakInd)*100;

off10Idx_Time = (fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & endField_1.off10==0);
cellMapTemp1 = timeMap_Def1.off10(off10Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.off10 = sum(peakInd<34)/length(peakInd)*100;

on30Idx_Time = (fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & endField_1.on30==0);
cellMapTemp1 = timeMap_Def1.on30(on30Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.on30 = sum(peakInd<34)/length(peakInd)*100;

off30Idx_Time = (fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & endField_1.off30==0);
cellMapTemp1 = timeMap_Def1.off30(off30Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.off30 = sum(peakInd<34)/length(peakInd)*100;


peakTime.rat1079.on10_mean = nanmean(within5s.on10);
peakTime.rat1079.off10_mean = nanmean(within5s.off10);
peakTime.rat1079.on30_mean = nanmean(within5s.on30);
peakTime.rat1079.off30_mean = nanmean(within5s.off30);

peakTime.rat1079.on10_std = std(within5s.on10,'omitnan');
peakTime.rat1079.off10_std = std(within5s.off10,'omitnan');
peakTime.rat1079.on30_std = std(within5s.on30,'omitnan');
peakTime.rat1079.off30_std = std(within5s.off30,'omitnan');


% rat all
AnalyzeSes = ratAll;    
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
    % avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
        
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
%% on10, off10, on30, off30
% all time cells
on10Idx_Time = (fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & endField_1.on10==0);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.on10 = sum(peakInd<34)/length(peakInd)*100;

off10Idx_Time = (fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & endField_1.off10==0);
cellMapTemp1 = timeMap_Def1.off10(off10Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.off10 = sum(peakInd<34)/length(peakInd)*100;

on30Idx_Time = (fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & endField_1.on30==0);
cellMapTemp1 = timeMap_Def1.on30(on30Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.on30 = sum(peakInd<34)/length(peakInd)*100;

off30Idx_Time = (fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & endField_1.off30==0);
cellMapTemp1 = timeMap_Def1.off30(off30Idx_Time,:);
[~,peakInd] = max(cellMapTemp1,[],2);
within5s.off30 = sum(peakInd<34)/length(peakInd)*100;


peakTime.ratAll.on10_mean = nanmean(within5s.on10);
peakTime.ratAll.off10_mean = nanmean(within5s.off10);
peakTime.ratAll.on30_mean = nanmean(within5s.on30);
peakTime.ratAll.off30_mean = nanmean(within5s.off30);

peakTime.ratAll.on10_std = std(within5s.on10,'omitnan');
peakTime.ratAll.off10_std = std(within5s.off10,'omitnan');
peakTime.ratAll.on30_std = std(within5s.on30,'omitnan');
peakTime.ratAll.off30_std = std(within5s.off30,'omitnan');