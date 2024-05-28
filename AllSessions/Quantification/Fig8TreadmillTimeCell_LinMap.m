% it needs to be stable in each session
% and stable in 2 session combined
%
function Fig8TreadmillTimeCell_LinMap(inFile,AnalyzeSes)

close all

p.savePlot = 0;
p.writeToFile = 0;
p.avgRateThres = 0.5;

avgRate = [];
fieldLabel.allSession = [];
allClusterNum = 0;

% need to be adjusted to automatic detection later
p.linROI = [33, 46];
p.linChoice = [55, 64];
p.xTickLin = [16 40 50 60 68];
p.xTickLabel = {'Return','Dl','St','Ch','Rw'};

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
% initiate the map
map_All = [];
map_On = [];
map_Off = [];

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
    delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap.mat');
    load(delayFile);
    % load correlation values
    stabilityFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayFireStability_TrialbyTrial_2Session.mat');
    load(stabilityFile);
    % load arm rate file
    armRateFile = fullfile(sessInfo(i).mainDir,'Cell Property', 'Fig8TreadmillArmRate.mat');
    load(armRateFile);
    
    
    clusterNum = length(SpikeProp.AvgRate.Fig8Rate);
    rateLabel = SpikeProp.AvgRate.Fig8Rate>0.1 & SpikeProp.AvgRate.Fig8Rate<5;
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
            post_Map_Def1.(sessDirs{j}) = [post_Map_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).post_mapMean_Def1_2Session{k}(1:20)];
            
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
    
    sessDirs = sessInfo(i).sessDirs;
    
    spikeMap = cell(clusterNum,1);
    occupMap = cell(clusterNum,1);
    
    spikeMap_On = cell(clusterNum,1);
    spikeMap_Off = cell(clusterNum,1);
    
    occupMap_On = cell(clusterNum,1);
    occupMap_Off = cell(clusterNum,1);   
    
    spikeMap_L = cell(clusterNum,1);
    spikeMap_R = cell(clusterNum,1);
    
    occupMap_L = cell(clusterNum,1);
    occupMap_R = cell(clusterNum,1);  
    
    % lin map
    for j = 1:length(sessDirs)
        
        % take map for each delay condition
        if i == AnalyzeSes(1)
            map_Block.(sessDirs{j}) = [];
            map_Block_L.(sessDirs{j}) = [];
            map_Block_R.(sessDirs{j}) = [];
            time_Block_L.(sessDirs{j}) = [];
            time_Block_R.(sessDirs{j}) = [];
        end
        
        % load analyzed map per trial
        load(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'));        

        trialNum = length(ratesByECLR.ECLR);
 
         posLabel = zeros(trialNum,1);
         turnLabel = zeros(trialNum,1);
        for m = 1:trialNum
            if ratesByECLR.ECLR(m) == 1
                posLabel(m) = 2;
                turnLabel(m) = 1;
            elseif ratesByECLR.ECLR(m) == 2
                posLabel(m) = 1;
                turnLabel(m) = 1;
            elseif ratesByECLR.ECLR(m) == 3
                posLabel(m) = 1;
                turnLabel(m) = 0;
            elseif ratesByECLR.ECLR(m) == 4
                posLabel(m) = 2;
                turnLabel(m) = 0;
            else
                error('Turning label ERROR');
            end
        end
        
        % find left / right return arm
        leftReturnIdx = find(posLabel==1); % idx on all trials
        rightReturnIdx = find(posLabel==2); % idx on all trials
        % find Left / Right turn trials
        leftTurnIdx = find(turnLabel==1); % idx on all trials
        rightTurnIdx = find(turnLabel==0); % idx on all trials
                
        spikeMap_Temp = cell(clusterNum,1);
        occupMap_Temp = cell(clusterNum,1);
    
        spikeMap_Temp_L = cell(clusterNum,1);
        occupMap_Temp_L = cell(clusterNum,1);
        
        spikeMap_Temp_R = cell(clusterNum,1);
        occupMap_Temp_R = cell(clusterNum,1);
        
        time_Temp_L = cell(clusterNum,1);
        time_Temp_R = cell(clusterNum,1);
        
        % SIMPLE WAY: put all trials together without considering
        % differences
        % Accurate way: match left to left and right to right (return,
        % base, choice, reward)
        for n = 1:trialNum
            
%             maze_time_Temp = sum(ratesByECLR.occupmap{1,1}([1:32,47:end]));
%             delay_time_Temp = sum(ratesByECLR.occupmap{1,1}([33:46]));
%             trial_time_Temp = trial_time_Temp + sum((ratesByECLR.occupmap{1,1}));
    
            for k = 1:clusterNum
                spikeMap{k} = [spikeMap{k};ratesByECLR.spikemap{k,n}];
                occupMap{k} = [occupMap{k};ratesByECLR.occupmap{k,n}];
                
                if turnLabel(n) == 1
                    spikeMap_L{k} = [spikeMap_L{k};ratesByECLR.spikemap{k,n}];
                    occupMap_L{k} = [occupMap_L{k};ratesByECLR.occupmap{k,n}];
                else
                    spikeMap_R{k} = [spikeMap_R{k};ratesByECLR.spikemap{k,n}];
                    occupMap_R{k} = [occupMap_R{k};ratesByECLR.occupmap{k,n}];
                end
                
    
                spikeMap_Temp{k} = [spikeMap_Temp{k};ratesByECLR.spikemap{k,n}];
                occupMap_Temp{k} = [occupMap_Temp{k};ratesByECLR.occupmap{k,n}];
                
                if turnLabel(n) == 1
                    spikeMap_Temp_L{k} = [spikeMap_Temp_L{k};ratesByECLR.spikemap{k,n}];
                    occupMap_Temp_L{k} = [occupMap_Temp_L{k};ratesByECLR.occupmap{k,n}];
                    time_Temp_L{k} = [time_Temp_L{k};DelayFire.(sessDirs{j}).spikeRate1_Smooth{k}(n,:)];
                    
                else
                    spikeMap_Temp_R{k} = [spikeMap_Temp_R{k};ratesByECLR.spikemap{k,n}];
                    occupMap_Temp_R{k} = [occupMap_Temp_R{k};ratesByECLR.occupmap{k,n}];
                    time_Temp_R{k} = [time_Temp_R{k};DelayFire.(sessDirs{j}).spikeRate1_Smooth{k}(n,:)];
                end
                
                if contains(sessDirs{j},'on')
                    spikeMap_On{k} = [spikeMap_On{k};ratesByECLR.spikemap{k,n}];
                    occupMap_On{k} = [occupMap_On{k};ratesByECLR.occupmap{k,n}]; 
                elseif contains(sessDirs{j},'off')
                    spikeMap_Off{k} = [spikeMap_Off{k};ratesByECLR.spikemap{k,n}];
                    occupMap_Off{k} = [occupMap_Off{k};ratesByECLR.occupmap{k,n}]; 
                else 
                    error('Block type wrong')                    
                end
            end
        end

            
        for k = 1:clusterNum
            if rateLabel(k) == 1
                spikeMap_AllTrial = sum(spikeMap_Temp{k},1);
                occupMap_AllTrial = sum(occupMap_Temp{k},1)+0.0001;
                map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
                map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
                map_Block.(sessDirs{j}) = [map_Block.(sessDirs{j});map_Temp_2];
                
                spikeMap_AllTrial = sum(spikeMap_Temp_L{k},1);
                occupMap_AllTrial = sum(occupMap_Temp_L{k},1)+0.0001;
                map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
                map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
                map_Block_L.(sessDirs{j}) = [map_Block_L.(sessDirs{j});map_Temp_2];
                
                spikeMap_AllTrial = sum(spikeMap_Temp_R{k},1);
                occupMap_AllTrial = sum(occupMap_Temp_R{k},1)+0.0001;
                map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
                map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
                map_Block_R.(sessDirs{j}) = [map_Block_R.(sessDirs{j});map_Temp_2];
                
                time_Block_L.(sessDirs{j}) = [time_Block_L.(sessDirs{j});mean(time_Temp_L{k},1)];
                time_Block_R.(sessDirs{j}) = [time_Block_R.(sessDirs{j});mean(time_Temp_R{k},1)];
            end
        end

    end
    
    for k = 1:clusterNum
        if rateLabel(k) == 1
            spikeMap_AllTrial = sum(spikeMap{k},1);
            occupMap_AllTrial = sum(occupMap{k},1)+0.0001;
            map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
            map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
            map_All = [map_All;map_Temp_2];
            
            % treadmill on
            spikeMap_AllTrial = sum(spikeMap_On{k},1);
            occupMap_AllTrial = sum(occupMap_On{k},1)+0.0001;
            map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
            map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
            map_On = [map_On;map_Temp_2];
            
            % treadmill off
            spikeMap_AllTrial = sum(spikeMap_Off{k},1);
            occupMap_AllTrial = sum(occupMap_Off{k},1)+0.0001;
            map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
            map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
            map_Off = [map_Off;map_Temp_2];
        end
    end 
    
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

%% on10, off10, on30, off30
% time cell delay temporal map and spatial map
h1 = figure;
h1.Position = [100 100 1200 400];
% on10 off10
% on10Idx_Time = (fieldLabel_Def1.on10>0 & delay_onIdx);
on10Idx_Time = (delay_onIdx==1);
cellMapTemp1 = time_Block_L.on10_2(on10Idx_Time,:);

% off10Idx_Time = (fieldLabel_Def1.off10>0 & delay_onIdx & ~endField.off10);
cellMapTemp2 = time_Block_R.on10_2(on10Idx_Time,:);

plotTimeMap2(cellMapTemp1,cellMapTemp2,2,[1,4,1],[1,4,2],{'Delay start at barrier','Time cells on10'},{'-3','0','10','13'},{'Space'},{'-3','0','10','13'});

% off 10
% off10Idx_Time = (fieldLabel_Def1.off10>0 & delay_onIdx);
off10Idx_Time = (delay_onIdx==1);
cellMapTemp1 = time_Block_L.off10_2(off10Idx_Time,:);

% off30Idx_Time = (fieldLabel_Def1.off30>0 & delay_onIdx & ~endField.off30);
cellMapTemp2 = time_Block_R.off10_2(off10Idx_Time,:);
plotTimeMap2(cellMapTemp1,cellMapTemp2,2,[1,4,3],[1,4,4],{'Delay start at barrier','Time cells off10'},{'-3','0','30','33'},{'Space'},{'-3','0','30','33'});


h1 = figure;
h1.Position = [100 100 1200 400];
% on10 off10
% on10Idx_Time = (fieldLabel_Def1.on10>0 & delay_onIdx);
on10Idx_Time = (delay_onIdx==1);
cellMapTemp1 = time_Block_R.on30_2(on10Idx_Time,:);

% off10Idx_Time = (fieldLabel_Def1.off10>0 & delay_onIdx & ~endField.off10);
cellMapTemp2 = map_Block_R.on30_2(on10Idx_Time,[48:end]);

plotTimeMap2(cellMapTemp1,cellMapTemp2,2,[1,4,1],[1,4,2],{'Delay start at barrier','Time cells on10'},{'-3','0','10','13'},{'Space'},{'-3','0','10','13'});

% off 10
% off10Idx_Time = (fieldLabel_Def1.off10>0 & delay_onIdx);
off10Idx_Time = (delay_onIdx==1);
cellMapTemp1 = time_Block_R.off30_2(off10Idx_Time,:);

% off30Idx_Time = (fieldLabel_Def1.off30>0 & delay_onIdx & ~endField.off30);
cellMapTemp2 = map_Block_R.off30_2(off10Idx_Time,[48:end]);
plotTimeMap2(cellMapTemp1,cellMapTemp2,2,[1,4,3],[1,4,4],{'Delay start at barrier','Time cells off10'},{'-3','0','30','33'},{'Space'},{'-3','0','30','33'});

% figure
% Violin(on10_mean,1,'ShowData',false,'ViolinColor',[1,0,0]);
% Violin(off10_mean,2,'ShowData',false,'ViolinColor',[0,0,0]);
% Violin(on30_mean,3,'ShowData',false,'ViolinColor',[1,0,0]);
% Violin(off30_mean,4,'ShowData',false,'ViolinColor',[0,0,0]);
% set(gca, 'XTick', [1,2,3,4], 'XTickLabel', {'on10','off10','on30','off30'});
% 
% [h,pVal] = kstest2(on10_mean,off10_mean)
% [h,pVal] = kstest2(on30_mean,off30_mean)
% [h,pVal] = kstest2(on10_mean,on30_mean)
% [h,pVal] = kstest2(off10_mean,off30_mean)


if p.savePlot == 1
    figName = sprintf('%s%s%s',savedir,'\Onset-Barrier-Time cell in its own session-2Session combined');
    print(figName,'-dpng','-r600');
end

% time cells noemalized by all max
% h2 = figure(2);
% h2.Position = [100 100 1200 400];
%% on10, off10, on30, off30
% all time cells
h2 = figure;
h2.Position = [100 100 1200 400];
% on10 off10
on10Idx_Time = (fieldLabel_Def1.on10>0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
on10_mean = mean(cellMapTemp1,2);

off10Idx_Time = (fieldLabel_Def1.off10>0 & delay_onIdx);
cellMapTemp2 = timeMap_Def1.off10(off10Idx_Time,:);
off10_mean = mean(cellMapTemp2,2);

plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),...
    1,[1,4,1],[1,4,2],{'Delay start at barrier','Stable cells on10'},{'-3','0','10','13'},{'Stable cells off10'},{'-3','0','10','13'});

% on30 and off30
off10Idx_Time = (fieldLabel_Def1.on30>0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.on30(off10Idx_Time,:);
on30_mean = mean(cellMapTemp1,2);

off30Idx_Time = (fieldLabel_Def1.off30>0 & delay_onIdx);
cellMapTemp2 = timeMap_Def1.off30(off30Idx_Time,:);
off30_mean = mean(cellMapTemp2,2);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on30(off10Idx_Time,:),post_Map_Def1.on30(off10Idx_Time,:),pre_Map_Def1.off30(off30Idx_Time,:),post_Map_Def1.off30(off30Idx_Time,:),...
    1,[1,4,3],[1,4,4],{'Delay start at barrier','Stable cells on10'},{'-3','0','30','33'},{'Stable cells off10'},{'-3','0','30','33'});


%%
figure
% ramp up cell
% on10 off10
on10Idx_ramp = (fieldLabel_Def1.on10>0 & delay_onIdx & endField.on10);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_ramp,:);

off10Idx_ramp = (fieldLabel_Def1.off10>0 & delay_onIdx & endField.off10);
cellMapTemp2 = timeMap_Def1.off10(off10Idx_ramp,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_ramp,:),post_Map_Def1.on10(on10Idx_ramp,:),pre_Map_Def1.off10(off10Idx_ramp,:),post_Map_Def1.off10(off10Idx_ramp,:),...
    1,[1,4,1],[1,4,2],{'Delay start at barrier','ramp cells on10'},{'-3','0','10','13'},{'ramp cells off10'},{'-3','0','10','13'});


% on30 and off30
on30Idx_ramp = (fieldLabel_Def1.on30>0 & delay_onIdx & endField.on30);
cellMapTemp1 = timeMap_Def1.on30(on30Idx_ramp,:);

off30Idx_ramp = (fieldLabel_Def1.off30>0 & delay_onIdx & endField.off30);
cellMapTemp2 = timeMap_Def1.off30(off30Idx_ramp,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on30(on30Idx_ramp,:),post_Map_Def1.on30(on30Idx_ramp,:),pre_Map_Def1.off30(off30Idx_ramp,:),post_Map_Def1.off30(off30Idx_ramp,:),...
    1,[1,4,3],[1,4,4],{'Delay start at barrier','ramp cells on10'},{'-3','0','30','33'},{'ramp cells off10'},{'-3','0','30','33'});

if p.savePlot == 1
    figName = sprintf('%s%s%s',savedir,'\Onset-Barrier-Time cell in its own session-2Session combined');
    print(figName,'-dpng','-r600');
end

%% all time cells sorted by on10
h1 = figure;
h1.Position = [100 100 900 800];
% on10 
on10Idx_Time = (fieldLabel_Def1.on10>0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
on10_mean = mean(cellMapTemp1,2);

% on30 
off10Idx_Time = on10Idx_Time;
cellMapTemp2 = timeMap_Def1.on30(off10Idx_Time,:);
on30_mean = mean(cellMapTemp2,2);

plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.on30(off10Idx_Time,:),post_Map_Def1.on30(off10Idx_Time,:),...
    2,[2,4,1],[2,4,2],{'Delay start at barrier','Stable on10'},{'-3','0','10','13'},{'on30 by on10'},{'-3','0','30','33'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));

% off10 
off10Idx_Time = on10Idx_Time;
cellMapTemp2 = timeMap_Def1.off10(off10Idx_Time,:);

plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),...
    2,[2,4,1],[2,4,3],{'Delay start at barrier','Stable on10'},{'-3','0','10','13'},{'off10 by on10'},{'-3','0','10','13'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));

% off30 
off30Idx_Time = on10Idx_Time;
cellMapTemp2 = timeMap_Def1.off30(off30Idx_Time,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.off30(off30Idx_Time,:),post_Map_Def1.off30(off30Idx_Time,:),...
    2,[2,4,1],[2,4,4],{'Delay start at barrier','Stable on10'},{'-3','0','10','13'},{'off30 by on10'},{'-3','0','30','33'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));

% sort by off10
off10Idx_Time = (fieldLabel_Def1.off10>0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.off10(off10Idx_Time,:);
off10_mean = mean(cellMapTemp1,2);

off30Idx_Time = off10Idx_Time;
cellMapTemp2 = timeMap_Def1.off30(off30Idx_Time,:);
off30_mean = mean(cellMapTemp2,2);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),pre_Map_Def1.off30(off30Idx_Time,:),post_Map_Def1.off30(off30Idx_Time,:),...
    2,[2,4,5],[2,4,6],{'Delay start at barrier','Stable off10'},{'-3','0','10','13'},{'off30 by off10'},{'-3','0','30','33'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));


% on10 
on10Idx_Time = off10Idx_Time;
cellMapTemp2 = timeMap_Def1.on10(on10Idx_Time,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),...
    2,[2,4,5],[2,4,7],{'Delay start at barrier','Stable off10'},{'-3','0','10','13'},{'on10 by off10'},{'-3','0','10','13'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));


% on30 
off10Idx_Time = off10Idx_Time;
cellMapTemp2 = timeMap_Def1.on30(off10Idx_Time,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),pre_Map_Def1.on30(off10Idx_Time,:),post_Map_Def1.on30(off10Idx_Time,:),...
    2,[2,4,5],[2,4,8],{'Delay start at barrier','Stable off10'},{'-3','0','10','13'},{'on30 by off10'},{'-3','0','30','33'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));


%% all delay on but no stable cells sorted by on10
h1 = figure;
h1.Position = [100 100 900 800];
% on10 
on10Idx_Time = (fieldLabel_Def1.on10==0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
% on30 
off10Idx_Time = on10Idx_Time;
cellMapTemp2 = timeMap_Def1.on30(off10Idx_Time,:);

plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.on30(off10Idx_Time,:),post_Map_Def1.on30(off10Idx_Time,:),...
    2,[2,4,1],[2,4,2],{'Delay start at barrier','Unstable cells on10'},{'-3','0','10','13'},{'on30 by on10'},{'-3','0','30','33'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));


% off10 
off10Idx_Time = on10Idx_Time;
cellMapTemp2 = timeMap_Def1.off10(off10Idx_Time,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),...
    2,[2,4,1],[2,4,3],{'Delay start at barrier','Unstable on10'},{'-3','0','10','13'},{'off10 by on10'},{'-3','0','10','13'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));


% off30 
off30Idx_Time = on10Idx_Time;
cellMapTemp2 = timeMap_Def1.off30(off30Idx_Time,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.off30(off30Idx_Time,:),post_Map_Def1.off30(off30Idx_Time,:),...
    2,[2,4,1],[2,4,4],{'Delay start at barrier','Unstable on10'},{'-3','0','10','13'},{'off30 by on10'},{'-3','0','30','33'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));



% sort by off10
off10Idx_Time = (fieldLabel_Def1.off10==0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.off10(off10Idx_Time,:);

off30Idx_Time = off10Idx_Time;
cellMapTemp2 = timeMap_Def1.off30(off30Idx_Time,:);

plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),pre_Map_Def1.off30(off30Idx_Time,:),post_Map_Def1.off30(off30Idx_Time,:),...
    2,[2,4,5],[2,4,6],{'Delay start at barrier','Unstable cells off10'},{'-3','0','10','13'},{'off30 by off10'},{'-3','0','30','33'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));


% on10 
on10Idx_Time = off10Idx_Time;
cellMapTemp2 = timeMap_Def1.on10(on10Idx_Time,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),...
    2,[2,4,5],[2,4,7],{'Delay start at barrier','Unstable cells off10'},{'-3','0','10','13'},{'on10 by off10'},{'-3','0','10','13'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));


% on30 
off10Idx_Time = off10Idx_Time;
cellMapTemp2 = timeMap_Def1.on30(off10Idx_Time,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),pre_Map_Def1.on30(off10Idx_Time,:),post_Map_Def1.on30(off10Idx_Time,:),...
    2,[2,4,5],[2,4,8],{'Delay start at barrier','Unstable cells off10'},{'-3','0','10','13'},{'on30 by off10'},{'-3','0','30','33'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));

%% all delay on cells sorted by on10
h1 = figure;
h1.Position = [100 100 900 800];
% on10 
on10Idx_Time = (delay_onIdx==1);
on10_mean = rate_Delay(1,on10Idx_Time);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
% on30 
off10Idx_Time = on10Idx_Time;
on30_mean = rate_Delay(3,off10Idx_Time);
cellMapTemp2 = timeMap_Def1.on30(off10Idx_Time,:);

plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.on30(off10Idx_Time,:),post_Map_Def1.on30(off10Idx_Time,:),...
    2,[2,4,1],[2,4,2],{'Delay start at barrier','Delay active cells sort by on10'},{'-3','0','10','13'},{'on30 by on10'},{'-3','0','30','33'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));


% off10 
off10Idx_Time = on10Idx_Time;
off10_mean = rate_Delay(2,off10Idx_Time);
cellMapTemp2 = timeMap_Def1.off10(off10Idx_Time,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),...
    2,[2,4,1],[2,4,3],{'Delay start at barrier','Delay active cells sort by on10'},{'-3','0','10','13'},{'off10 by on10'},{'-3','0','10','13'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));


% off30 
off30Idx_Time = on10Idx_Time;
off30_mean = rate_Delay(4,off10Idx_Time);
cellMapTemp2 = timeMap_Def1.off30(off30Idx_Time,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.off30(off30Idx_Time,:),post_Map_Def1.off30(off30Idx_Time,:),...
    2,[2,4,1],[2,4,4],{'Delay start at barrier','Delay active cells sort by on10'},{'-3','0','10','13'},{'off30 by on10'},{'-3','0','30','33'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));


% sort by off10
off10Idx_Time = (delay_onIdx==1);
cellMapTemp1 = timeMap_Def1.off10(off10Idx_Time,:);

off30Idx_Time = off10Idx_Time;
cellMapTemp2 = timeMap_Def1.off30(off30Idx_Time,:);

plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),pre_Map_Def1.off30(off30Idx_Time,:),post_Map_Def1.off30(off30Idx_Time,:),...
    2,[2,4,5],[2,4,6],{'Delay start at barrier','Delay active cells sort by off10'},{'-3','0','10','13'},{'off30 by off10'},{'-3','0','30','33'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));


% on10 
on10Idx_Time = off10Idx_Time;
cellMapTemp2 = timeMap_Def1.on10(on10Idx_Time,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),...
    2,[2,4,5],[2,4,7],{'Delay start at barrier','Delay active cells sort by off10'},{'-3','0','10','13'},{'on10 by off10'},{'-3','0','10','13'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));


% on30 
off10Idx_Time = off10Idx_Time;
cellMapTemp2 = timeMap_Def1.on30(off10Idx_Time,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),pre_Map_Def1.on30(off10Idx_Time,:),post_Map_Def1.on30(off10Idx_Time,:),...
    2,[2,4,5],[2,4,8],{'Delay start at barrier','Delay active cells sort by off10'},{'-3','0','10','13'},{'on30 by off10'},{'-3','0','30','33'});
sizeTemp = size(cellMapTemp1);
corrVal = corr2(cellMapTemp1,cellMapTemp2(1:sizeTemp(1),1:sizeTemp(2)));
xlabel(sprintf('corr = %1.3f',corrVal));

% figure
figure
plot(on10_mean,on30_mean,'ro')
hold on
plot(off10_mean,off30_mean,'k^')
figure
plot(on10_mean,off10_mean,'m*')
hold on
plot(on30_mean,off30_mean,'b*')

%% all delay-on cell 
figure
% on 10
% on10 off10
on10Idx_Time = delay_onIdx_1;
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);

off10Idx_Time = delay_onIdx_2;
cellMapTemp2 = timeMap_Def1.off10(off10Idx_Time,:);

plotTimeMap2(cellMapTemp1,cellMapTemp2,1,[1,4,1],[1,4,2],...
    {'Delay start at barrier','Delay on cells on10'},{'0','10'},{'Delay on cells off10'},{'0','10'});

% on30 and off30
off10Idx_Time = delay_onIdx_3;
cellMapTemp1 = timeMap_Def1.on30(off10Idx_Time,:);

off30Idx_Time = delay_onIdx_4;
cellMapTemp2 = timeMap_Def1.off30(off30Idx_Time,:);
plotTimeMap2(cellMapTemp1,cellMapTemp2,1,[1,4,3],[1,4,4],...
    {'Delay on cells on30'},{'0','30'},{'Delay on cells off30'},{'0','30'});

%% all non-time cell
figure
% on 10
on10Idx = (~on10Idx_Time & delay_onIdx);
cellMapTemp1 = timeMap_Def1.on10(on10Idx,:);
on10_mean = mean(rate_Delay(1,on10Idx),1);

[~,mean_Sort] = sort(on10_mean);
cellMapTemp1_Sort = cellMapTemp1(mean_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);
subplot(1,4,1)
imagesc(cellMapTemp1_Sort)
% colormap(jet)
set(gca, 'xtick', [0 size(cellMapTemp1_Sort,2)]);
set(gca, 'xticklabels', [0 10]);
caxis([0 5])


% off 10
off10Idx = (~off10Idx_Time & delay_onIdx);
cellMapTemp1 = timeMap_Def1.off10(off10Idx,:);
off10_mean = mean(rate_Delay(2,off10Idx),1);

[~,mean_Sort] = sort(off10_mean);
cellMapTemp1_Sort = cellMapTemp1(mean_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);
subplot(1,4,2)
imagesc(cellMapTemp1_Sort)
colormap(jet)
set(gca, 'xtick', [0 size(cellMapTemp1_Sort,2)]);
set(gca, 'xticklabels', [0 10]);
caxis([0 5])


% on30
on30Idx = (~off10Idx_Time & delay_onIdx);
cellMapTemp1 = timeMap_Def1.on30(on30Idx,:);
on30_mean = mean(rate_Delay(3,on30Idx),1);

[~,mean_Sort] = sort(on30_mean);
cellMapTemp1_Sort = cellMapTemp1(mean_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);
subplot(1,4,3)
imagesc(cellMapTemp1_Sort)
colormap(jet)
set(gca, 'xtick', [0 size(cellMapTemp1_Sort,2)]);
set(gca, 'xticklabels', [0 30]);
caxis([0 5])


% off 30
off30Idx = (~off30Idx_Time & delay_onIdx);
cellMapTemp1 = timeMap_Def1.off30(off30Idx,:);
off30_mean = mean(rate_Delay(4,off30Idx),1);

[~,mean_Sort] = sort(off30_mean);
cellMapTemp1_Sort = cellMapTemp1(mean_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);
subplot(1,4,4)
imagesc(cellMapTemp1_Sort)
colormap(jet)
set(gca, 'xtick', [0 size(cellMapTemp1_Sort,2)]);
set(gca, 'xticklabels', [0 30]);
caxis([0 5])

figure
Violin(on10_mean,1,'ShowData',false,'ViolinColor',[1,0,0]);
Violin(off10_mean,2,'ShowData',false,'ViolinColor',[0,0,0]);
Violin(on30_mean,3,'ShowData',false,'ViolinColor',[1,0,0]);
Violin(off30_mean,4,'ShowData',false,'ViolinColor',[0,0,0]);
set(gca, 'XTick', [1,2,3,4], 'XTickLabel', {'on10','off10','on30','off30'});


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