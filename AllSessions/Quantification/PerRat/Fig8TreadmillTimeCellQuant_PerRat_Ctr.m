% it needs to be stable in each session
% and stable in 2 session combined
%
clear all
inFile = 'Fig8Treadmill_On.xlsx';

close all

p.savePlot = 0;
p.writeToFile = 0;

% specifiy each rat session
rat1056 = 1:4;
rat1075 = 5:8;
rat1077 = 9:14;
rat1081 = 15:18;
ratAll = 1:18;

p.savePlot = 1;
p.writeToFile = 1;
p.avgRateThres = 0.5;

avgRate = [];
fieldLabel.allSession = [];
allClusterNum = 0;

% % Read in input information
sessInfo = SessInfoImport(inFile);
  
% rat 1056
AnalyzeSes = rat1056;
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
    
    rate_DelayTemp = [Fig8TreadmillArmRate.on10.rateDelay(idx,1)';Fig8TreadmillArmRate.on30.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];
    delay_onIdxTemp = (sum(rate_DelayTemp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
    
end

delay_onIdx_1 = sum(rate_Delay(1,:)>p.avgRateThres,1)>0;
delay_onIdx_2 = sum(rate_Delay(2,:)>p.avgRateThres,1)>0;

delayActiveRatio.rat1056.on10 = 100*sum(delay_onIdx_1)/length(delay_onIdx_1);
delayActiveRatio.rat1056.on30 = 100*sum(delay_onIdx_2)/length(delay_onIdx_2);
delayActiveRatio.rat1056.All = 100*sum(delay_onIdx)/length(delay_onIdx);
delayActiveRatio.rat1056.overlap.on = 2*100*sum((delay_onIdx_1+delay_onIdx_2)==2)/sum(delay_onIdx_1+delay_onIdx_2);

timeCell.rat1056.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx)/length(delay_onIdx);
timeCell.rat1056.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx)/length(delay_onIdx);

time_limited_Cell.rat1056.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & ~endField_1.on10)/length(delay_onIdx);
time_limited_Cell.rat1056.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & ~endField_1.on30)/length(delay_onIdx);

consistentCell.rat1056.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & endField_1.on10)/length(delay_onIdx);
consistentCell.rat1056.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & endField_1.on30)/length(delay_onIdx);


% rat 1075
AnalyzeSes = rat1075;
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
    
    rate_DelayTemp = [Fig8TreadmillArmRate.on10.rateDelay(idx,1)';Fig8TreadmillArmRate.on30.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];
    delay_onIdxTemp = (sum(rate_DelayTemp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
    
end

delay_onIdx_1 = sum(rate_Delay(1,:)>p.avgRateThres,1)>0;
delay_onIdx_2 = sum(rate_Delay(2,:)>p.avgRateThres,1)>0;

delayActiveRatio.rat1075.on10 = 100*sum(delay_onIdx_1)/length(delay_onIdx_1);
delayActiveRatio.rat1075.on30 = 100*sum(delay_onIdx_2)/length(delay_onIdx_2);
delayActiveRatio.rat1075.All = 100*sum(delay_onIdx)/length(delay_onIdx);
delayActiveRatio.rat1075.overlap.on = 2*100*sum((delay_onIdx_1+delay_onIdx_2)==2)/sum(delay_onIdx_1+delay_onIdx_2);

timeCell.rat1075.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx)/length(delay_onIdx);
timeCell.rat1075.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx)/length(delay_onIdx);

time_limited_Cell.rat1075.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & ~endField_1.on10)/length(delay_onIdx);
time_limited_Cell.rat1075.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & ~endField_1.on30)/length(delay_onIdx);

consistentCell.rat1075.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & endField_1.on10)/length(delay_onIdx);
consistentCell.rat1075.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & endField_1.on30)/length(delay_onIdx);


% rat 1077
AnalyzeSes = rat1077;
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
    
    rate_DelayTemp = [Fig8TreadmillArmRate.on10.rateDelay(idx,1)';Fig8TreadmillArmRate.on30.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];
    delay_onIdxTemp = (sum(rate_DelayTemp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
    
end

delay_onIdx_1 = sum(rate_Delay(1,:)>p.avgRateThres,1)>0;
delay_onIdx_2 = sum(rate_Delay(2,:)>p.avgRateThres,1)>0;

delayActiveRatio.rat1077.on10 = 100*sum(delay_onIdx_1)/length(delay_onIdx_1);
delayActiveRatio.rat1077.on30 = 100*sum(delay_onIdx_2)/length(delay_onIdx_2);
delayActiveRatio.rat1077.All = 100*sum(delay_onIdx)/length(delay_onIdx);
delayActiveRatio.rat1077.overlap.on = 2*100*sum((delay_onIdx_1+delay_onIdx_2)==2)/sum(delay_onIdx_1+delay_onIdx_2);

timeCell.rat1077.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx)/length(delay_onIdx);
timeCell.rat1077.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx)/length(delay_onIdx);

time_limited_Cell.rat1077.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & ~endField_1.on10)/length(delay_onIdx);
time_limited_Cell.rat1077.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & ~endField_1.on30)/length(delay_onIdx);

consistentCell.rat1077.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & endField_1.on10)/length(delay_onIdx);
consistentCell.rat1077.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & endField_1.on30)/length(delay_onIdx);


% rat 1081
AnalyzeSes = rat1081;
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
    
    rate_DelayTemp = [Fig8TreadmillArmRate.on10.rateDelay(idx,1)';Fig8TreadmillArmRate.on30.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];
    delay_onIdxTemp = (sum(rate_DelayTemp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
    
end

delay_onIdx_1 = sum(rate_Delay(1,:)>p.avgRateThres,1)>0;
delay_onIdx_2 = sum(rate_Delay(2,:)>p.avgRateThres,1)>0;

delayActiveRatio.rat1081.on10 = 100*sum(delay_onIdx_1)/length(delay_onIdx_1);
delayActiveRatio.rat1081.on30 = 100*sum(delay_onIdx_2)/length(delay_onIdx_2);
delayActiveRatio.rat1081.All = 100*sum(delay_onIdx)/length(delay_onIdx);
delayActiveRatio.rat1081.overlap.on = 2*100*sum((delay_onIdx_1+delay_onIdx_2)==2)/sum(delay_onIdx_1+delay_onIdx_2);

timeCell.rat1081.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx)/length(delay_onIdx);
timeCell.rat1081.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx)/length(delay_onIdx);

time_limited_Cell.rat1081.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & ~endField_1.on10)/length(delay_onIdx);
time_limited_Cell.rat1081.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & ~endField_1.on30)/length(delay_onIdx);

consistentCell.rat1081.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & endField_1.on10)/length(delay_onIdx);
consistentCell.rat1081.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & endField_1.on30)/length(delay_onIdx);


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
    
    rate_DelayTemp = [Fig8TreadmillArmRate.on10.rateDelay(idx,1)';Fig8TreadmillArmRate.on30.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];
    delay_onIdxTemp = (sum(rate_DelayTemp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
    
end

delay_onIdx_1 = sum(rate_Delay(1,:)>p.avgRateThres,1)>0;
delay_onIdx_2 = sum(rate_Delay(2,:)>p.avgRateThres,1)>0;

delayActiveRatio.ratAll.on10 = 100*sum(delay_onIdx_1)/length(delay_onIdx_1);
delayActiveRatio.ratAll.on30 = 100*sum(delay_onIdx_2)/length(delay_onIdx_2);
delayActiveRatio.ratAll.All = 100*sum(delay_onIdx)/length(delay_onIdx);
delayActiveRatio.ratAll.overlap.on = 2*100*sum((delay_onIdx_1+delay_onIdx_2)==2)/sum(delay_onIdx_1+delay_onIdx_2);

timeCell.ratAll.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx)/length(delay_onIdx);
timeCell.ratAll.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx)/length(delay_onIdx);

time_limited_Cell.ratAll.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & ~endField_1.on10)/length(delay_onIdx);
time_limited_Cell.ratAll.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & ~endField_1.on30)/length(delay_onIdx);

consistentCell.ratAll.on10 = 100*sum(fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & endField_1.on10)/length(delay_onIdx);
consistentCell.ratAll.on30 = 100*sum(fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & endField_1.on30)/length(delay_onIdx);


%% off control group

clear all

rat1069 = 1:6;
rat1071 = 7:13;
rat1072 = 14:19;
rat1076 = 20:26;
ratAll = 1:26;

inFile = 'Fig8Treadmill_Off.xlsx';

p.avgRateThres = 0.5;

avgRate = [];
fieldLabel.allSession = [];
allClusterNum = 0;

% % Read in input information
sessInfo = SessInfoImport(inFile);
  
% rat 1069
AnalyzeSes = rat1069;
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
    
    rate_DelayTemp = [Fig8TreadmillArmRate.off10.rateDelay(idx,1)';Fig8TreadmillArmRate.off30.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];
    delay_onIdxTemp = (sum(rate_DelayTemp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
    
end

delay_onIdx_1 = sum(rate_Delay(1,:)>p.avgRateThres,1)>0;
delay_onIdx_2 = sum(rate_Delay(2,:)>p.avgRateThres,1)>0;

delayActiveRatio.rat1069.off10 = 100*sum(delay_onIdx_1)/length(delay_onIdx_1);
delayActiveRatio.rat1069.off30 = 100*sum(delay_onIdx_2)/length(delay_onIdx_2);
delayActiveRatio.rat1069.All = 100*sum(delay_onIdx)/length(delay_onIdx);
delayActiveRatio.rat1069.overlap.on = 2*100*sum((delay_onIdx_1+delay_onIdx_2)==2)/sum(delay_onIdx_1+delay_onIdx_2);

timeCell.rat1069.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx)/length(delay_onIdx);
timeCell.rat1069.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx)/length(delay_onIdx);

time_limited_Cell.rat1069.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & ~endField_1.off10)/length(delay_onIdx);
time_limited_Cell.rat1069.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & ~endField_1.off30)/length(delay_onIdx);

consistentCell.rat1069.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & endField_1.off10)/length(delay_onIdx);
consistentCell.rat1069.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & endField_1.off30)/length(delay_onIdx);


% rat 1071
AnalyzeSes = rat1071;
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
    
    rate_DelayTemp = [Fig8TreadmillArmRate.off10.rateDelay(idx,1)';Fig8TreadmillArmRate.off30.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];
    delay_onIdxTemp = (sum(rate_DelayTemp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
    
end

delay_onIdx_1 = sum(rate_Delay(1,:)>p.avgRateThres,1)>0;
delay_onIdx_2 = sum(rate_Delay(2,:)>p.avgRateThres,1)>0;

delayActiveRatio.rat1071.off10 = 100*sum(delay_onIdx_1)/length(delay_onIdx_1);
delayActiveRatio.rat1071.off30 = 100*sum(delay_onIdx_2)/length(delay_onIdx_2);
delayActiveRatio.rat1071.All = 100*sum(delay_onIdx)/length(delay_onIdx);
delayActiveRatio.rat1071.overlap.on = 2*100*sum((delay_onIdx_1+delay_onIdx_2)==2)/sum(delay_onIdx_1+delay_onIdx_2);

timeCell.rat1071.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx)/length(delay_onIdx);
timeCell.rat1071.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx)/length(delay_onIdx);

time_limited_Cell.rat1071.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & ~endField_1.off10)/length(delay_onIdx);
time_limited_Cell.rat1071.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & ~endField_1.off30)/length(delay_onIdx);

consistentCell.rat1071.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & endField_1.off10)/length(delay_onIdx);
consistentCell.rat1071.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & endField_1.off30)/length(delay_onIdx);


% rat 1072
AnalyzeSes = rat1072;
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
    
    rate_DelayTemp = [Fig8TreadmillArmRate.off10.rateDelay(idx,1)';Fig8TreadmillArmRate.off10.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];
    delay_onIdxTemp = (sum(rate_DelayTemp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
    
end

delay_onIdx_1 = sum(rate_Delay(1,:)>p.avgRateThres,1)>0;
delay_onIdx_2 = sum(rate_Delay(2,:)>p.avgRateThres,1)>0;

delayActiveRatio.rat1072.off10 = 100*sum(delay_onIdx_1)/length(delay_onIdx_1);
delayActiveRatio.rat1072.off30 = 100*sum(delay_onIdx_2)/length(delay_onIdx_2);
delayActiveRatio.rat1072.All = 100*sum(delay_onIdx)/length(delay_onIdx);
delayActiveRatio.rat1072.overlap.on = 2*100*sum((delay_onIdx_1+delay_onIdx_2)==2)/sum(delay_onIdx_1+delay_onIdx_2);

timeCell.rat1072.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx)/length(delay_onIdx);
timeCell.rat1072.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx)/length(delay_onIdx);

time_limited_Cell.rat1072.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & ~endField_1.off10)/length(delay_onIdx);
time_limited_Cell.rat1072.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & ~endField_1.off30)/length(delay_onIdx);

consistentCell.rat1072.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & endField_1.off10)/length(delay_onIdx);
consistentCell.rat1072.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & endField_1.off30)/length(delay_onIdx);



% rat 1076
AnalyzeSes = rat1076;
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
    
    rate_DelayTemp = [Fig8TreadmillArmRate.off10.rateDelay(idx,1)';Fig8TreadmillArmRate.off30.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];
    delay_onIdxTemp = (sum(rate_DelayTemp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
    
end

delay_onIdx_1 = sum(rate_Delay(1,:)>p.avgRateThres,1)>0;
delay_onIdx_2 = sum(rate_Delay(2,:)>p.avgRateThres,1)>0;

delayActiveRatio.rat1076.off10 = 100*sum(delay_onIdx_1)/length(delay_onIdx_1);
delayActiveRatio.rat1076.off30 = 100*sum(delay_onIdx_2)/length(delay_onIdx_2);
delayActiveRatio.rat1076.All = 100*sum(delay_onIdx)/length(delay_onIdx);
delayActiveRatio.rat1076.overlap.on = 2*100*sum((delay_onIdx_1+delay_onIdx_2)==2)/sum(delay_onIdx_1+delay_onIdx_2);

timeCell.rat1076.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx)/length(delay_onIdx);
timeCell.rat1076.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx)/length(delay_onIdx);

time_limited_Cell.rat1076.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & ~endField_1.off10)/length(delay_onIdx);
time_limited_Cell.rat1076.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & ~endField_1.off30)/length(delay_onIdx);

consistentCell.rat1076.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & endField_1.off10)/length(delay_onIdx);
consistentCell.rat1076.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & endField_1.off30)/length(delay_onIdx);



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
    
    rate_DelayTemp = [Fig8TreadmillArmRate.off10.rateDelay(idx,1)';Fig8TreadmillArmRate.off30.rateDelay(idx,1)'];
    rate_Delay = [rate_Delay,rate_DelayTemp];
    delay_onIdxTemp = (sum(rate_DelayTemp>p.avgRateThres,1))>0;
    delay_onIdx = [delay_onIdx,delay_onIdxTemp];
    
end

delay_onIdx_1 = sum(rate_Delay(1,:)>p.avgRateThres,1)>0;
delay_onIdx_2 = sum(rate_Delay(2,:)>p.avgRateThres,1)>0;

delayActiveRatio.ratAll.off10 = 100*sum(delay_onIdx_1)/length(delay_onIdx_1);
delayActiveRatio.ratAll.off30 = 100*sum(delay_onIdx_2)/length(delay_onIdx_2);
delayActiveRatio.ratAll.All = 100*sum(delay_onIdx)/length(delay_onIdx);
delayActiveRatio.ratAll.overlap.on = 2*100*sum((delay_onIdx_1+delay_onIdx_2)==2)/sum(delay_onIdx_1+delay_onIdx_2);

timeCell.ratAll.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx)/length(delay_onIdx);
timeCell.ratAll.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx)/length(delay_onIdx);

time_limited_Cell.ratAll.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & ~endField_1.off10)/length(delay_onIdx);
time_limited_Cell.ratAll.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & ~endField_1.off30)/length(delay_onIdx);

consistentCell.ratAll.off10 = 100*sum(fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx & endField_1.off10)/length(delay_onIdx);
consistentCell.ratAll.off30 = 100*sum(fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx & endField_1.off30)/length(delay_onIdx);