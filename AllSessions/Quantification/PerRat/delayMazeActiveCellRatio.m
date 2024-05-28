% Plot time % being used in each delay conditions
clear all
inFile = 'Fig8Treadmill_OnOff.xlsx';

close all

p.avgRateThres = 0.5;

% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'

% Read in input information
sessInfo = SessInfoImport(inFile);

% specifiy each rat session
rat1043 = 1:4;
rat1044 = 5:8;
rat1046 = 9:12;
rat1058 = 13:15;
rat1079 = 16:18;
ratAll = 1:18;

% rat 1043
AnalyzeSes = rat1043;
mazeRatio.on10 = []; mazeRatio.off10 = []; mazeRatio.on30 = []; mazeRatio.off30 = []; mazeRatio.All = [];
delayRatio.on10 = []; delayRatio.off10 = []; delayRatio.on30 = []; delayRatio.off30 = []; delayRatio.All = [];
delayRatio.overlap.on = []; delayRatio.overlap.off = []; delayRatio.overlap.length_10 = []; delayRatio.overlap.length_30 = [];

for i = AnalyzeSes(1:end)
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    rateCluster = SpikeProp.AvgRate.Fig8Rate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    for j = 1:length(sessDirs)
        
            maze_time.(sessDirs{j}) = [];
            delay_time.(sessDirs{j}) = [];
            block_time.(sessDirs{j}) = [];
            
            mazeSpike.(sessDirs{j}) = [];
            delaySpike.(sessDirs{j}) = [];
        
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);       
  
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
        %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
        %         delayTend2 = Fig8DelayZonePos.delayPos2.endT;
        
        trialNum = size(delayTstart1,2);
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
            delay_timeTemp = 10*trialNum;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
            delay_timeTemp = 30*trialNum;
        else
            error('Delay time is wrong')
        end
  
        block_timeTemp = pathData.t(end) - pathData.t(1);
        maze_timeTemp = pathData.t(end) - pathData.t(1) - delay_timeTemp;
        
         for k = 1:clusterNum             % get each spike time, change unit to msec from sec             % ts unit: ms             tSp = Spike_Session.(sessDirs{j}){k};             if rateLabel(k) == 1                 if ~isempty(tSp)                     % only consider pupative pyramidal cells                      ts_DelayCount = 0;                     for m = 1:trialNum                                                 % delay definition 1 start from barrier                         ts_DelayCount = ts_DelayCount + sum(tSp>delayTstart1(m) & tSp<delayTend1_2(m));                     end                                          delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),ts_DelayCount];                     mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),length(tSp)-ts_DelayCount];                     block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];                     maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];                     delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];                 else                     delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),0];                     mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),0];                     block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];                     maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];                     delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];                 end             end          end
            % get each spike time, change unit to msec from sec
            % ts unit: ms
            tSp = Spike_Session.(sessDirs{j}){k};
            if rateLabel(k) == 1
                if ~isempty(tSp)
                    % only consider pupative pyramidal cells

                    ts_DelayCount = 0;
                    for m = 1:trialNum                        
                        % delay definition 1 start from barrier
                        ts_DelayCount = ts_DelayCount + sum(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                    end
                    
                    delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),ts_DelayCount];
                    mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),length(tSp)-ts_DelayCount];
                    block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];
                    maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];
                    delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];
                else
                    delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),0];
                    mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),0];
                    block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];
                    maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];
                    delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];
                end
            end
         end 
    end
  
    % time used in 10 or 30 sec condition
    mazeRate.on10 = (mazeSpike.on10_1 + mazeSpike.on10_2)./(maze_time.on10_1 + maze_time.on10_2);
    mazeRate.off10 = (mazeSpike.off10_1 + mazeSpike.off10_2)./(maze_time.off10_1 + maze_time.off10_2);
    mazeRate.on30 = (mazeSpike.on30_1 + mazeSpike.on30_2)./(maze_time.on30_1 + maze_time.on30_2);
    mazeRate.off30 = (mazeSpike.off30_1 + mazeSpike.off30_2)./(maze_time.off30_1 + maze_time.off30_2);
    maze_onIdx_on10 = mazeRate.on10>p.avgRateThres;
    maze_onIdx_off10 = mazeRate.off10>p.avgRateThres;
    maze_onIdx_on30 = mazeRate.on30>p.avgRateThres;
    maze_onIdx_off30 = mazeRate.off30>p.avgRateThres;
    maze_onIdx = (maze_onIdx_on10+maze_onIdx_off10+maze_onIdx_on30+maze_onIdx_off30)>0;
    
    delayRate.on10 = (delaySpike.on10_1 + delaySpike.on10_2)./(delay_time.on10_1 + delay_time.on10_2);
    delayRate.off10 = (delaySpike.off10_1 + delaySpike.off10_2)./(delay_time.off10_1 + delay_time.off10_2);
    delayRate.on30 = (delaySpike.on30_1 + delaySpike.on30_2)./(delay_time.on30_1 + delay_time.on30_2);
    delayRate.off30 = (delaySpike.off30_1 + delaySpike.off30_2)./(delay_time.off30_1 + delay_time.off30_2);
    delay_onIdx_on10 = delayRate.on10>p.avgRateThres;
    delay_onIdx_off10 = delayRate.off10>p.avgRateThres;
    delay_onIdx_on30 = delayRate.on30>p.avgRateThres;
    delay_onIdx_off30 = delayRate.off30>p.avgRateThres;
    delay_onIdx = (delay_onIdx_on10+delay_onIdx_off10+delay_onIdx_on30+delay_onIdx_off30)>0;
    
    mazeRatio.on10 = [mazeRatio.on10,100*sum(maze_onIdx_on10)/length(maze_onIdx_on10)];
    mazeRatio.off10 = [mazeRatio.off10,100*sum(maze_onIdx_off10)/length(maze_onIdx_off10)];
    mazeRatio.on30 = [mazeRatio.on30,100*sum(maze_onIdx_on30)/length(maze_onIdx_on30)];
    mazeRatio.off30 = [mazeRatio.off30,100*sum(maze_onIdx_off30)/length(maze_onIdx_off30)];
    mazeRatio.All = [mazeRatio.All,100*sum(maze_onIdx)/length(maze_onIdx)];
    
    delayRatio.on10 = [delayRatio.on10,100*sum(delay_onIdx_on10)/length(delay_onIdx_on10)];
    delayRatio.off10 = [delayRatio.off10,100*sum(delay_onIdx_off10)/length(delay_onIdx_off10)];
    delayRatio.on30 = [delayRatio.on30,100*sum(delay_onIdx_on30)/length(delay_onIdx_on30)];
    delayRatio.off30 = [delayRatio.off30,100*sum(delay_onIdx_off30)/length(delay_onIdx_off30)];
    delayRatio.All = [delayRatio.All,100*sum(delay_onIdx)/length(delay_onIdx)];
    
    delayRatio.overlap.on = [delayRatio.overlap.on,2*100*sum((delay_onIdx_on10+delay_onIdx_on30)==2)/sum(delay_onIdx_on10+delay_onIdx_on30)];
    delayRatio.overlap.off = [delayRatio.overlap.off,2*100*sum((delay_onIdx_off10+delay_onIdx_off30)==2)/sum(delay_onIdx_off10+delay_onIdx_off30)];
    delayRatio.overlap.length_10 = [delayRatio.overlap.length_10,2*100*sum((delay_onIdx_on10+delay_onIdx_off10)==2)/sum(delay_onIdx_on10+delay_onIdx_off10)];
    delayRatio.overlap.length_30 = [delayRatio.overlap.length_30,2*100*sum((delay_onIdx_on30+delay_onIdx_off30)==2)/sum(delay_onIdx_on30+delay_onIdx_off30)];
end

mazeActiveRatio.rat1043.on10_nanmean = nanmean(mazeRatio.on10);
mazeActiveRatio.rat1043.off10_nanmean = nanmean(mazeRatio.off10);
mazeActiveRatio.rat1043.on30_nanmean = nanmean(mazeRatio.on30);
mazeActiveRatio.rat1043.off30_nanmean = nanmean(mazeRatio.off30);
mazeActiveRatio.rat1043.All_nanmean = nanmean(mazeRatio.All);

mazeActiveRatio.rat1043.on10_std = std(mazeRatio.on10);
mazeActiveRatio.rat1043.off10_std = std(mazeRatio.off10);
mazeActiveRatio.rat1043.on30_std = std(mazeRatio.on30);
mazeActiveRatio.rat1043.off30_std = std(mazeRatio.off30);
mazeActiveRatio.rat1043.All_std = std(mazeRatio.All);

delayActiveRatio.rat1043.on10_nanmean = nanmean(delayRatio.on10);
delayActiveRatio.rat1043.off10_nanmean = nanmean(delayRatio.off10);
delayActiveRatio.rat1043.on30_nanmean = nanmean(delayRatio.on30);
delayActiveRatio.rat1043.off30_nanmean = nanmean(delayRatio.off30);
delayActiveRatio.rat1043.All_nanmean = nanmean(delayRatio.All);

delayActiveRatio.rat1043.on10_std = std(delayRatio.on10);
delayActiveRatio.rat1043.off10_std = std(delayRatio.off10);
delayActiveRatio.rat1043.on30_std = std(delayRatio.on30);
delayActiveRatio.rat1043.off30_std = std(delayRatio.off30);
delayActiveRatio.rat1043.All_std = std(delayRatio.All);

delayActiveRatio.rat1043.overlap.on_nanmean = nanmean(delayRatio.overlap.on);
delayActiveRatio.rat1043.overlap.off_nanmean = nanmean(delayRatio.overlap.off);
delayActiveRatio.rat1043.overlap.length_10_nanmean = nanmean(delayRatio.overlap.length_10);
delayActiveRatio.rat1043.overlap.length_30_nanmean = nanmean(delayRatio.overlap.length_30);

delayActiveRatio.rat1043.overlap.on_std = std(delayRatio.overlap.on);
delayActiveRatio.rat1043.overlap.off_std = std(delayRatio.overlap.off);
delayActiveRatio.rat1043.overlap.length_10_std = std(delayRatio.overlap.length_10);
delayActiveRatio.rat1043.overlap.length_30_std = std(delayRatio.overlap.length_30);


% rat 1044
AnalyzeSes = rat1044;
mazeRatio.on10 = []; mazeRatio.off10 = []; mazeRatio.on30 = []; mazeRatio.off30 = []; mazeRatio.All = [];
delayRatio.on10 = []; delayRatio.off10 = []; delayRatio.on30 = []; delayRatio.off30 = []; delayRatio.All = [];
delayRatio.overlap.on = []; delayRatio.overlap.off = []; delayRatio.overlap.length_10 = []; delayRatio.overlap.length_30 = [];

for i = AnalyzeSes(1:end)
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    rateCluster = SpikeProp.AvgRate.Fig8Rate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    for j = 1:length(sessDirs)
        
            maze_time.(sessDirs{j}) = [];
            delay_time.(sessDirs{j}) = [];
            block_time.(sessDirs{j}) = [];
            
            mazeSpike.(sessDirs{j}) = [];
            delaySpike.(sessDirs{j}) = [];
        
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);       
  
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
        %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
        %         delayTend2 = Fig8DelayZonePos.delayPos2.endT;
        
        trialNum = size(delayTstart1,2);
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
            delay_timeTemp = 10*trialNum;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
            delay_timeTemp = 30*trialNum;
        else
            error('Delay time is wrong')
        end
  
        block_timeTemp = pathData.t(end) - pathData.t(1);
        maze_timeTemp = pathData.t(end) - pathData.t(1) - delay_timeTemp;
        
         for k = 1:clusterNum             % get each spike time, change unit to msec from sec             % ts unit: ms             tSp = Spike_Session.(sessDirs{j}){k};             if rateLabel(k) == 1                 if ~isempty(tSp)                     % only consider pupative pyramidal cells                      ts_DelayCount = 0;                     for m = 1:trialNum                                                 % delay definition 1 start from barrier                         ts_DelayCount = ts_DelayCount + sum(tSp>delayTstart1(m) & tSp<delayTend1_2(m));                     end                                          delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),ts_DelayCount];                     mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),length(tSp)-ts_DelayCount];                     block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];                     maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];                     delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];                 else                     delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),0];                     mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),0];                     block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];                     maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];                     delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];                 end             end          end
            % get each spike time, change unit to msec from sec
            % ts unit: ms
            tSp = Spike_Session.(sessDirs{j}){k};
            if rateLabel(k) == 1
                if ~isempty(tSp)
                    % only consider pupative pyramidal cells

                    ts_DelayCount = 0;
                    for m = 1:trialNum                        
                        % delay definition 1 start from barrier
                        ts_DelayCount = ts_DelayCount + sum(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                    end
                    
                    delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),ts_DelayCount];
                    mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),length(tSp)-ts_DelayCount];
                    block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];
                    maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];
                    delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];
                else
                    delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),0];
                    mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),0];
                    block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];
                    maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];
                    delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];
                end
            end
         end 
    end
  
    % time used in 10 or 30 sec condition
    mazeRate.on10 = (mazeSpike.on10_1 + mazeSpike.on10_2)./(maze_time.on10_1 + maze_time.on10_2);
    mazeRate.off10 = (mazeSpike.off10_1 + mazeSpike.off10_2)./(maze_time.off10_1 + maze_time.off10_2);
    mazeRate.on30 = (mazeSpike.on30_1 + mazeSpike.on30_2)./(maze_time.on30_1 + maze_time.on30_2);
    mazeRate.off30 = (mazeSpike.off30_1 + mazeSpike.off30_2)./(maze_time.off30_1 + maze_time.off30_2);
    maze_onIdx_on10 = mazeRate.on10>p.avgRateThres;
    maze_onIdx_off10 = mazeRate.off10>p.avgRateThres;
    maze_onIdx_on30 = mazeRate.on30>p.avgRateThres;
    maze_onIdx_off30 = mazeRate.off30>p.avgRateThres;
    maze_onIdx = (maze_onIdx_on10+maze_onIdx_off10+maze_onIdx_on30+maze_onIdx_off30)>0;
    
    delayRate.on10 = (delaySpike.on10_1 + delaySpike.on10_2)./(delay_time.on10_1 + delay_time.on10_2);
    delayRate.off10 = (delaySpike.off10_1 + delaySpike.off10_2)./(delay_time.off10_1 + delay_time.off10_2);
    delayRate.on30 = (delaySpike.on30_1 + delaySpike.on30_2)./(delay_time.on30_1 + delay_time.on30_2);
    delayRate.off30 = (delaySpike.off30_1 + delaySpike.off30_2)./(delay_time.off30_1 + delay_time.off30_2);
    delay_onIdx_on10 = delayRate.on10>p.avgRateThres;
    delay_onIdx_off10 = delayRate.off10>p.avgRateThres;
    delay_onIdx_on30 = delayRate.on30>p.avgRateThres;
    delay_onIdx_off30 = delayRate.off30>p.avgRateThres;
    delay_onIdx = (delay_onIdx_on10+delay_onIdx_off10+delay_onIdx_on30+delay_onIdx_off30)>0;
    
    mazeRatio.on10 = [mazeRatio.on10,100*sum(maze_onIdx_on10)/length(maze_onIdx_on10)];
    mazeRatio.off10 = [mazeRatio.off10,100*sum(maze_onIdx_off10)/length(maze_onIdx_off10)];
    mazeRatio.on30 = [mazeRatio.on30,100*sum(maze_onIdx_on30)/length(maze_onIdx_on30)];
    mazeRatio.off30 = [mazeRatio.off30,100*sum(maze_onIdx_off30)/length(maze_onIdx_off30)];
    mazeRatio.All = [mazeRatio.All,100*sum(maze_onIdx)/length(maze_onIdx)];
    
    delayRatio.on10 = [delayRatio.on10,100*sum(delay_onIdx_on10)/length(delay_onIdx_on10)];
    delayRatio.off10 = [delayRatio.off10,100*sum(delay_onIdx_off10)/length(delay_onIdx_off10)];
    delayRatio.on30 = [delayRatio.on30,100*sum(delay_onIdx_on30)/length(delay_onIdx_on30)];
    delayRatio.off30 = [delayRatio.off30,100*sum(delay_onIdx_off30)/length(delay_onIdx_off30)];
    delayRatio.All = [delayRatio.All,100*sum(delay_onIdx)/length(delay_onIdx)];
    
    delayRatio.overlap.on = [delayRatio.overlap.on,2*100*sum((delay_onIdx_on10+delay_onIdx_on30)==2)/sum(delay_onIdx_on10+delay_onIdx_on30)];
    delayRatio.overlap.off = [delayRatio.overlap.off,2*100*sum((delay_onIdx_off10+delay_onIdx_off30)==2)/sum(delay_onIdx_off10+delay_onIdx_off30)];
    delayRatio.overlap.length_10 = [delayRatio.overlap.length_10,2*100*sum((delay_onIdx_on10+delay_onIdx_off10)==2)/sum(delay_onIdx_on10+delay_onIdx_off10)];
    delayRatio.overlap.length_30 = [delayRatio.overlap.length_30,2*100*sum((delay_onIdx_on30+delay_onIdx_off30)==2)/sum(delay_onIdx_on30+delay_onIdx_off30)];
end

mazeActiveRatio.rat1044.on10_nanmean = nanmean(mazeRatio.on10);
mazeActiveRatio.rat1044.off10_nanmean = nanmean(mazeRatio.off10);
mazeActiveRatio.rat1044.on30_nanmean = nanmean(mazeRatio.on30);
mazeActiveRatio.rat1044.off30_nanmean = nanmean(mazeRatio.off30);
mazeActiveRatio.rat1044.All_nanmean = nanmean(mazeRatio.All);

mazeActiveRatio.rat1044.on10_std = std(mazeRatio.on10);
mazeActiveRatio.rat1044.off10_std = std(mazeRatio.off10);
mazeActiveRatio.rat1044.on30_std = std(mazeRatio.on30);
mazeActiveRatio.rat1044.off30_std = std(mazeRatio.off30);
mazeActiveRatio.rat1044.All_std = std(mazeRatio.All);

delayActiveRatio.rat1044.on10_nanmean = nanmean(delayRatio.on10);
delayActiveRatio.rat1044.off10_nanmean = nanmean(delayRatio.off10);
delayActiveRatio.rat1044.on30_nanmean = nanmean(delayRatio.on30);
delayActiveRatio.rat1044.off30_nanmean = nanmean(delayRatio.off30);
delayActiveRatio.rat1044.All_nanmean = nanmean(delayRatio.All);

delayActiveRatio.rat1044.on10_std = std(delayRatio.on10);
delayActiveRatio.rat1044.off10_std = std(delayRatio.off10);
delayActiveRatio.rat1044.on30_std = std(delayRatio.on30);
delayActiveRatio.rat1044.off30_std = std(delayRatio.off30);
delayActiveRatio.rat1044.All_std = std(delayRatio.All);

delayActiveRatio.rat1044.overlap.on_nanmean = nanmean(delayRatio.overlap.on);
delayActiveRatio.rat1044.overlap.off_nanmean = nanmean(delayRatio.overlap.off);
delayActiveRatio.rat1044.overlap.length_10_nanmean = nanmean(delayRatio.overlap.length_10);
delayActiveRatio.rat1044.overlap.length_30_nanmean = nanmean(delayRatio.overlap.length_30);

delayActiveRatio.rat1044.overlap.on_std = std(delayRatio.overlap.on);
delayActiveRatio.rat1044.overlap.off_std = std(delayRatio.overlap.off);
delayActiveRatio.rat1044.overlap.length_10_std = std(delayRatio.overlap.length_10);
delayActiveRatio.rat1044.overlap.length_30_std = std(delayRatio.overlap.length_30);


% rat 1046
AnalyzeSes = rat1046;
mazeRatio.on10 = []; mazeRatio.off10 = []; mazeRatio.on30 = []; mazeRatio.off30 = []; mazeRatio.All = [];
delayRatio.on10 = []; delayRatio.off10 = []; delayRatio.on30 = []; delayRatio.off30 = []; delayRatio.All = [];
delayRatio.overlap.on = []; delayRatio.overlap.off = []; delayRatio.overlap.length_10 = []; delayRatio.overlap.length_30 = [];

for i = AnalyzeSes(1:end)
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    rateCluster = SpikeProp.AvgRate.Fig8Rate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    for j = 1:length(sessDirs)
        
            maze_time.(sessDirs{j}) = [];
            delay_time.(sessDirs{j}) = [];
            block_time.(sessDirs{j}) = [];
            
            mazeSpike.(sessDirs{j}) = [];
            delaySpike.(sessDirs{j}) = [];
        
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);       
  
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
        %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
        %         delayTend2 = Fig8DelayZonePos.delayPos2.endT;
        
        trialNum = size(delayTstart1,2);
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
            delay_timeTemp = 10*trialNum;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
            delay_timeTemp = 30*trialNum;
        else
            error('Delay time is wrong')
        end
  
        block_timeTemp = pathData.t(end) - pathData.t(1);
        maze_timeTemp = pathData.t(end) - pathData.t(1) - delay_timeTemp;
        
         for k = 1:clusterNum             % get each spike time, change unit to msec from sec             % ts unit: ms             tSp = Spike_Session.(sessDirs{j}){k};             if rateLabel(k) == 1                 if ~isempty(tSp)                     % only consider pupative pyramidal cells                      ts_DelayCount = 0;                     for m = 1:trialNum                                                 % delay definition 1 start from barrier                         ts_DelayCount = ts_DelayCount + sum(tSp>delayTstart1(m) & tSp<delayTend1_2(m));                     end                                          delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),ts_DelayCount];                     mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),length(tSp)-ts_DelayCount];                     block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];                     maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];                     delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];                 else                     delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),0];                     mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),0];                     block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];                     maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];                     delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];                 end             end          end
            % get each spike time, change unit to msec from sec
            % ts unit: ms
            tSp = Spike_Session.(sessDirs{j}){k};
            if rateLabel(k) == 1
                if ~isempty(tSp)
                    % only consider pupative pyramidal cells

                    ts_DelayCount = 0;
                    for m = 1:trialNum                        
                        % delay definition 1 start from barrier
                        ts_DelayCount = ts_DelayCount + sum(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                    end
                    
                    delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),ts_DelayCount];
                    mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),length(tSp)-ts_DelayCount];
                    block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];
                    maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];
                    delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];
                else
                    delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),0];
                    mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),0];
                    block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];
                    maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];
                    delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];
                end
            end
         end 
    end
  
    % time used in 10 or 30 sec condition
    mazeRate.on10 = (mazeSpike.on10_1 + mazeSpike.on10_2)./(maze_time.on10_1 + maze_time.on10_2);
    mazeRate.off10 = (mazeSpike.off10_1 + mazeSpike.off10_2)./(maze_time.off10_1 + maze_time.off10_2);
    mazeRate.on30 = (mazeSpike.on30_1 + mazeSpike.on30_2)./(maze_time.on30_1 + maze_time.on30_2);
    mazeRate.off30 = (mazeSpike.off30_1 + mazeSpike.off30_2)./(maze_time.off30_1 + maze_time.off30_2);
    maze_onIdx_on10 = mazeRate.on10>p.avgRateThres;
    maze_onIdx_off10 = mazeRate.off10>p.avgRateThres;
    maze_onIdx_on30 = mazeRate.on30>p.avgRateThres;
    maze_onIdx_off30 = mazeRate.off30>p.avgRateThres;
    maze_onIdx = (maze_onIdx_on10+maze_onIdx_off10+maze_onIdx_on30+maze_onIdx_off30)>0;
    
    delayRate.on10 = (delaySpike.on10_1 + delaySpike.on10_2)./(delay_time.on10_1 + delay_time.on10_2);
    delayRate.off10 = (delaySpike.off10_1 + delaySpike.off10_2)./(delay_time.off10_1 + delay_time.off10_2);
    delayRate.on30 = (delaySpike.on30_1 + delaySpike.on30_2)./(delay_time.on30_1 + delay_time.on30_2);
    delayRate.off30 = (delaySpike.off30_1 + delaySpike.off30_2)./(delay_time.off30_1 + delay_time.off30_2);
    delay_onIdx_on10 = delayRate.on10>p.avgRateThres;
    delay_onIdx_off10 = delayRate.off10>p.avgRateThres;
    delay_onIdx_on30 = delayRate.on30>p.avgRateThres;
    delay_onIdx_off30 = delayRate.off30>p.avgRateThres;
    delay_onIdx = (delay_onIdx_on10+delay_onIdx_off10+delay_onIdx_on30+delay_onIdx_off30)>0;
    
    mazeRatio.on10 = [mazeRatio.on10,100*sum(maze_onIdx_on10)/length(maze_onIdx_on10)];
    mazeRatio.off10 = [mazeRatio.off10,100*sum(maze_onIdx_off10)/length(maze_onIdx_off10)];
    mazeRatio.on30 = [mazeRatio.on30,100*sum(maze_onIdx_on30)/length(maze_onIdx_on30)];
    mazeRatio.off30 = [mazeRatio.off30,100*sum(maze_onIdx_off30)/length(maze_onIdx_off30)];
    mazeRatio.All = [mazeRatio.All,100*sum(maze_onIdx)/length(maze_onIdx)];
    
    delayRatio.on10 = [delayRatio.on10,100*sum(delay_onIdx_on10)/length(delay_onIdx_on10)];
    delayRatio.off10 = [delayRatio.off10,100*sum(delay_onIdx_off10)/length(delay_onIdx_off10)];
    delayRatio.on30 = [delayRatio.on30,100*sum(delay_onIdx_on30)/length(delay_onIdx_on30)];
    delayRatio.off30 = [delayRatio.off30,100*sum(delay_onIdx_off30)/length(delay_onIdx_off30)];
    delayRatio.All = [delayRatio.All,100*sum(delay_onIdx)/length(delay_onIdx)];
    
    delayRatio.overlap.on = [delayRatio.overlap.on,2*100*sum((delay_onIdx_on10+delay_onIdx_on30)==2)/sum(delay_onIdx_on10+delay_onIdx_on30)];
    delayRatio.overlap.off = [delayRatio.overlap.off,2*100*sum((delay_onIdx_off10+delay_onIdx_off30)==2)/sum(delay_onIdx_off10+delay_onIdx_off30)];
    delayRatio.overlap.length_10 = [delayRatio.overlap.length_10,2*100*sum((delay_onIdx_on10+delay_onIdx_off10)==2)/sum(delay_onIdx_on10+delay_onIdx_off10)];
    delayRatio.overlap.length_30 = [delayRatio.overlap.length_30,2*100*sum((delay_onIdx_on30+delay_onIdx_off30)==2)/sum(delay_onIdx_on30+delay_onIdx_off30)];
end

mazeActiveRatio.rat1046.on10_nanmean = nanmean(mazeRatio.on10);
mazeActiveRatio.rat1046.off10_nanmean = nanmean(mazeRatio.off10);
mazeActiveRatio.rat1046.on30_nanmean = nanmean(mazeRatio.on30);
mazeActiveRatio.rat1046.off30_nanmean = nanmean(mazeRatio.off30);
mazeActiveRatio.rat1046.All_nanmean = nanmean(mazeRatio.All);

mazeActiveRatio.rat1046.on10_std = std(mazeRatio.on10);
mazeActiveRatio.rat1046.off10_std = std(mazeRatio.off10);
mazeActiveRatio.rat1046.on30_std = std(mazeRatio.on30);
mazeActiveRatio.rat1046.off30_std = std(mazeRatio.off30);
mazeActiveRatio.rat1046.All_std = std(mazeRatio.All);

delayActiveRatio.rat1046.on10_nanmean = nanmean(delayRatio.on10);
delayActiveRatio.rat1046.off10_nanmean = nanmean(delayRatio.off10);
delayActiveRatio.rat1046.on30_nanmean = nanmean(delayRatio.on30);
delayActiveRatio.rat1046.off30_nanmean = nanmean(delayRatio.off30);
delayActiveRatio.rat1046.All_nanmean = nanmean(delayRatio.All);

delayActiveRatio.rat1046.on10_std = std(delayRatio.on10);
delayActiveRatio.rat1046.off10_std = std(delayRatio.off10);
delayActiveRatio.rat1046.on30_std = std(delayRatio.on30);
delayActiveRatio.rat1046.off30_std = std(delayRatio.off30);
delayActiveRatio.rat1046.All_std = std(delayRatio.All);

delayActiveRatio.rat1046.overlap.on_nanmean = nanmean(delayRatio.overlap.on);
delayActiveRatio.rat1046.overlap.off_nanmean = nanmean(delayRatio.overlap.off);
delayActiveRatio.rat1046.overlap.length_10_nanmean = nanmean(delayRatio.overlap.length_10);
delayActiveRatio.rat1046.overlap.length_30_nanmean = nanmean(delayRatio.overlap.length_30);

delayActiveRatio.rat1046.overlap.on_std = std(delayRatio.overlap.on);
delayActiveRatio.rat1046.overlap.off_std = std(delayRatio.overlap.off);
delayActiveRatio.rat1046.overlap.length_10_std = std(delayRatio.overlap.length_10);
delayActiveRatio.rat1046.overlap.length_30_std = std(delayRatio.overlap.length_30);


% rat 1058
AnalyzeSes = rat1058;
mazeRatio.on10 = []; mazeRatio.off10 = []; mazeRatio.on30 = []; mazeRatio.off30 = []; mazeRatio.All = [];
delayRatio.on10 = []; delayRatio.off10 = []; delayRatio.on30 = []; delayRatio.off30 = []; delayRatio.All = [];
delayRatio.overlap.on = []; delayRatio.overlap.off = []; delayRatio.overlap.length_10 = []; delayRatio.overlap.length_30 = [];

for i = AnalyzeSes(1:end)
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    rateCluster = SpikeProp.AvgRate.Fig8Rate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    for j = 1:length(sessDirs)
        
            maze_time.(sessDirs{j}) = [];
            delay_time.(sessDirs{j}) = [];
            block_time.(sessDirs{j}) = [];
            
            mazeSpike.(sessDirs{j}) = [];
            delaySpike.(sessDirs{j}) = [];
        
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);       
  
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
        %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
        %         delayTend2 = Fig8DelayZonePos.delayPos2.endT;
        
        trialNum = size(delayTstart1,2);
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
            delay_timeTemp = 10*trialNum;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
            delay_timeTemp = 30*trialNum;
        else
            error('Delay time is wrong')
        end
  
        block_timeTemp = pathData.t(end) - pathData.t(1);
        maze_timeTemp = pathData.t(end) - pathData.t(1) - delay_timeTemp;
        
         for k = 1:clusterNum             % get each spike time, change unit to msec from sec             % ts unit: ms             tSp = Spike_Session.(sessDirs{j}){k};             if rateLabel(k) == 1                 if ~isempty(tSp)                     % only consider pupative pyramidal cells                      ts_DelayCount = 0;                     for m = 1:trialNum                                                 % delay definition 1 start from barrier                         ts_DelayCount = ts_DelayCount + sum(tSp>delayTstart1(m) & tSp<delayTend1_2(m));                     end                                          delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),ts_DelayCount];                     mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),length(tSp)-ts_DelayCount];                     block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];                     maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];                     delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];                 else                     delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),0];                     mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),0];                     block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];                     maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];                     delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];                 end             end          end
            % get each spike time, change unit to msec from sec
            % ts unit: ms
            tSp = Spike_Session.(sessDirs{j}){k};
            if rateLabel(k) == 1
                if ~isempty(tSp)
                    % only consider pupative pyramidal cells

                    ts_DelayCount = 0;
                    for m = 1:trialNum                        
                        % delay definition 1 start from barrier
                        ts_DelayCount = ts_DelayCount + sum(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                    end
                    
                    delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),ts_DelayCount];
                    mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),length(tSp)-ts_DelayCount];
                    block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];
                    maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];
                    delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];
                else
                    delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),0];
                    mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),0];
                    block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];
                    maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];
                    delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];
                end
            end
         end 
    end
  
    % time used in 10 or 30 sec condition
    mazeRate.on10 = (mazeSpike.on10_1 + mazeSpike.on10_2)./(maze_time.on10_1 + maze_time.on10_2);
    mazeRate.off10 = (mazeSpike.off10_1 + mazeSpike.off10_2)./(maze_time.off10_1 + maze_time.off10_2);
    mazeRate.on30 = (mazeSpike.on30_1 + mazeSpike.on30_2)./(maze_time.on30_1 + maze_time.on30_2);
    mazeRate.off30 = (mazeSpike.off30_1 + mazeSpike.off30_2)./(maze_time.off30_1 + maze_time.off30_2);
    maze_onIdx_on10 = mazeRate.on10>p.avgRateThres;
    maze_onIdx_off10 = mazeRate.off10>p.avgRateThres;
    maze_onIdx_on30 = mazeRate.on30>p.avgRateThres;
    maze_onIdx_off30 = mazeRate.off30>p.avgRateThres;
    maze_onIdx = (maze_onIdx_on10+maze_onIdx_off10+maze_onIdx_on30+maze_onIdx_off30)>0;
    
    delayRate.on10 = (delaySpike.on10_1 + delaySpike.on10_2)./(delay_time.on10_1 + delay_time.on10_2);
    delayRate.off10 = (delaySpike.off10_1 + delaySpike.off10_2)./(delay_time.off10_1 + delay_time.off10_2);
    delayRate.on30 = (delaySpike.on30_1 + delaySpike.on30_2)./(delay_time.on30_1 + delay_time.on30_2);
    delayRate.off30 = (delaySpike.off30_1 + delaySpike.off30_2)./(delay_time.off30_1 + delay_time.off30_2);
    delay_onIdx_on10 = delayRate.on10>p.avgRateThres;
    delay_onIdx_off10 = delayRate.off10>p.avgRateThres;
    delay_onIdx_on30 = delayRate.on30>p.avgRateThres;
    delay_onIdx_off30 = delayRate.off30>p.avgRateThres;
    delay_onIdx = (delay_onIdx_on10+delay_onIdx_off10+delay_onIdx_on30+delay_onIdx_off30)>0;
    
    mazeRatio.on10 = [mazeRatio.on10,100*sum(maze_onIdx_on10)/length(maze_onIdx_on10)];
    mazeRatio.off10 = [mazeRatio.off10,100*sum(maze_onIdx_off10)/length(maze_onIdx_off10)];
    mazeRatio.on30 = [mazeRatio.on30,100*sum(maze_onIdx_on30)/length(maze_onIdx_on30)];
    mazeRatio.off30 = [mazeRatio.off30,100*sum(maze_onIdx_off30)/length(maze_onIdx_off30)];
    mazeRatio.All = [mazeRatio.All,100*sum(maze_onIdx)/length(maze_onIdx)];
    
    delayRatio.on10 = [delayRatio.on10,100*sum(delay_onIdx_on10)/length(delay_onIdx_on10)];
    delayRatio.off10 = [delayRatio.off10,100*sum(delay_onIdx_off10)/length(delay_onIdx_off10)];
    delayRatio.on30 = [delayRatio.on30,100*sum(delay_onIdx_on30)/length(delay_onIdx_on30)];
    delayRatio.off30 = [delayRatio.off30,100*sum(delay_onIdx_off30)/length(delay_onIdx_off30)];
    delayRatio.All = [delayRatio.All,100*sum(delay_onIdx)/length(delay_onIdx)];
    
    delayRatio.overlap.on = [delayRatio.overlap.on,2*100*sum((delay_onIdx_on10+delay_onIdx_on30)==2)/sum(delay_onIdx_on10+delay_onIdx_on30)];
    delayRatio.overlap.off = [delayRatio.overlap.off,2*100*sum((delay_onIdx_off10+delay_onIdx_off30)==2)/sum(delay_onIdx_off10+delay_onIdx_off30)];
    delayRatio.overlap.length_10 = [delayRatio.overlap.length_10,2*100*sum((delay_onIdx_on10+delay_onIdx_off10)==2)/sum(delay_onIdx_on10+delay_onIdx_off10)];
    delayRatio.overlap.length_30 = [delayRatio.overlap.length_30,2*100*sum((delay_onIdx_on30+delay_onIdx_off30)==2)/sum(delay_onIdx_on30+delay_onIdx_off30)];
end

mazeActiveRatio.rat1058.on10_nanmean = nanmean(mazeRatio.on10);
mazeActiveRatio.rat1058.off10_nanmean = nanmean(mazeRatio.off10);
mazeActiveRatio.rat1058.on30_nanmean = nanmean(mazeRatio.on30);
mazeActiveRatio.rat1058.off30_nanmean = nanmean(mazeRatio.off30);
mazeActiveRatio.rat1058.All_nanmean = nanmean(mazeRatio.All);

mazeActiveRatio.rat1058.on10_std = std(mazeRatio.on10);
mazeActiveRatio.rat1058.off10_std = std(mazeRatio.off10);
mazeActiveRatio.rat1058.on30_std = std(mazeRatio.on30);
mazeActiveRatio.rat1058.off30_std = std(mazeRatio.off30);
mazeActiveRatio.rat1058.All_std = std(mazeRatio.All);

delayActiveRatio.rat1058.on10_nanmean = nanmean(delayRatio.on10);
delayActiveRatio.rat1058.off10_nanmean = nanmean(delayRatio.off10);
delayActiveRatio.rat1058.on30_nanmean = nanmean(delayRatio.on30);
delayActiveRatio.rat1058.off30_nanmean = nanmean(delayRatio.off30);
delayActiveRatio.rat1058.All_nanmean = nanmean(delayRatio.All);

delayActiveRatio.rat1058.on10_std = std(delayRatio.on10);
delayActiveRatio.rat1058.off10_std = std(delayRatio.off10);
delayActiveRatio.rat1058.on30_std = std(delayRatio.on30);
delayActiveRatio.rat1058.off30_std = std(delayRatio.off30);
delayActiveRatio.rat1058.All_std = std(delayRatio.All);

delayActiveRatio.rat1058.overlap.on_nanmean = nanmean(delayRatio.overlap.on);
delayActiveRatio.rat1058.overlap.off_nanmean = nanmean(delayRatio.overlap.off);
delayActiveRatio.rat1058.overlap.length_10_nanmean = nanmean(delayRatio.overlap.length_10);
delayActiveRatio.rat1058.overlap.length_30_nanmean = nanmean(delayRatio.overlap.length_30);

delayActiveRatio.rat1058.overlap.on_std = std(delayRatio.overlap.on);
delayActiveRatio.rat1058.overlap.off_std = std(delayRatio.overlap.off);
delayActiveRatio.rat1058.overlap.length_10_std = std(delayRatio.overlap.length_10);
delayActiveRatio.rat1058.overlap.length_30_std = std(delayRatio.overlap.length_30);


% rat 1079
AnalyzeSes = rat1079;
mazeRatio.on10 = []; mazeRatio.off10 = []; mazeRatio.on30 = []; mazeRatio.off30 = []; mazeRatio.All = [];
delayRatio.on10 = []; delayRatio.off10 = []; delayRatio.on30 = []; delayRatio.off30 = []; delayRatio.All = [];
delayRatio.overlap.on = []; delayRatio.overlap.off = []; delayRatio.overlap.length_10 = []; delayRatio.overlap.length_30 = [];

for i = AnalyzeSes(1:end)
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    rateCluster = SpikeProp.AvgRate.Fig8Rate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    for j = 1:length(sessDirs)
        
            maze_time.(sessDirs{j}) = [];
            delay_time.(sessDirs{j}) = [];
            block_time.(sessDirs{j}) = [];
            
            mazeSpike.(sessDirs{j}) = [];
            delaySpike.(sessDirs{j}) = [];
        
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);       
  
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
        %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
        %         delayTend2 = Fig8DelayZonePos.delayPos2.endT;
        
        trialNum = size(delayTstart1,2);
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
            delay_timeTemp = 10*trialNum;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
            delay_timeTemp = 30*trialNum;
        else
            error('Delay time is wrong')
        end
  
        block_timeTemp = pathData.t(end) - pathData.t(1);
        maze_timeTemp = pathData.t(end) - pathData.t(1) - delay_timeTemp;
        
         for k = 1:clusterNum             % get each spike time, change unit to msec from sec             % ts unit: ms             tSp = Spike_Session.(sessDirs{j}){k};             if rateLabel(k) == 1                 if ~isempty(tSp)                     % only consider pupative pyramidal cells                      ts_DelayCount = 0;                     for m = 1:trialNum                                                 % delay definition 1 start from barrier                         ts_DelayCount = ts_DelayCount + sum(tSp>delayTstart1(m) & tSp<delayTend1_2(m));                     end                                          delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),ts_DelayCount];                     mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),length(tSp)-ts_DelayCount];                     block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];                     maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];                     delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];                 else                     delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),0];                     mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),0];                     block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];                     maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];                     delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];                 end             end          end
            % get each spike time, change unit to msec from sec
            % ts unit: ms
            tSp = Spike_Session.(sessDirs{j}){k};
            if rateLabel(k) == 1
                if ~isempty(tSp)
                    % only consider pupative pyramidal cells

                    ts_DelayCount = 0;
                    for m = 1:trialNum                        
                        % delay definition 1 start from barrier
                        ts_DelayCount = ts_DelayCount + sum(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                    end
                    
                    delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),ts_DelayCount];
                    mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),length(tSp)-ts_DelayCount];
                    block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];
                    maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];
                    delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];
                else
                    delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),0];
                    mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),0];
                    block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];
                    maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];
                    delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];
                end
            end
         end 
    end
  
    % time used in 10 or 30 sec condition
    mazeRate.on10 = (mazeSpike.on10_1 + mazeSpike.on10_2)./(maze_time.on10_1 + maze_time.on10_2);
    mazeRate.off10 = (mazeSpike.off10_1 + mazeSpike.off10_2)./(maze_time.off10_1 + maze_time.off10_2);
    mazeRate.on30 = (mazeSpike.on30_1 + mazeSpike.on30_2)./(maze_time.on30_1 + maze_time.on30_2);
    mazeRate.off30 = (mazeSpike.off30_1 + mazeSpike.off30_2)./(maze_time.off30_1 + maze_time.off30_2);
    maze_onIdx_on10 = mazeRate.on10>p.avgRateThres;
    maze_onIdx_off10 = mazeRate.off10>p.avgRateThres;
    maze_onIdx_on30 = mazeRate.on30>p.avgRateThres;
    maze_onIdx_off30 = mazeRate.off30>p.avgRateThres;
    maze_onIdx = (maze_onIdx_on10+maze_onIdx_off10+maze_onIdx_on30+maze_onIdx_off30)>0;
    
    delayRate.on10 = (delaySpike.on10_1 + delaySpike.on10_2)./(delay_time.on10_1 + delay_time.on10_2);
    delayRate.off10 = (delaySpike.off10_1 + delaySpike.off10_2)./(delay_time.off10_1 + delay_time.off10_2);
    delayRate.on30 = (delaySpike.on30_1 + delaySpike.on30_2)./(delay_time.on30_1 + delay_time.on30_2);
    delayRate.off30 = (delaySpike.off30_1 + delaySpike.off30_2)./(delay_time.off30_1 + delay_time.off30_2);
    delay_onIdx_on10 = delayRate.on10>p.avgRateThres;
    delay_onIdx_off10 = delayRate.off10>p.avgRateThres;
    delay_onIdx_on30 = delayRate.on30>p.avgRateThres;
    delay_onIdx_off30 = delayRate.off30>p.avgRateThres;
    delay_onIdx = (delay_onIdx_on10+delay_onIdx_off10+delay_onIdx_on30+delay_onIdx_off30)>0;
    
    mazeRatio.on10 = [mazeRatio.on10,100*sum(maze_onIdx_on10)/length(maze_onIdx_on10)];
    mazeRatio.off10 = [mazeRatio.off10,100*sum(maze_onIdx_off10)/length(maze_onIdx_off10)];
    mazeRatio.on30 = [mazeRatio.on30,100*sum(maze_onIdx_on30)/length(maze_onIdx_on30)];
    mazeRatio.off30 = [mazeRatio.off30,100*sum(maze_onIdx_off30)/length(maze_onIdx_off30)];
    mazeRatio.All = [mazeRatio.All,100*sum(maze_onIdx)/length(maze_onIdx)];
    
    delayRatio.on10 = [delayRatio.on10,100*sum(delay_onIdx_on10)/length(delay_onIdx_on10)];
    delayRatio.off10 = [delayRatio.off10,100*sum(delay_onIdx_off10)/length(delay_onIdx_off10)];
    delayRatio.on30 = [delayRatio.on30,100*sum(delay_onIdx_on30)/length(delay_onIdx_on30)];
    delayRatio.off30 = [delayRatio.off30,100*sum(delay_onIdx_off30)/length(delay_onIdx_off30)];
    delayRatio.All = [delayRatio.All,100*sum(delay_onIdx)/length(delay_onIdx)];
    
    delayRatio.overlap.on = [delayRatio.overlap.on,2*100*sum((delay_onIdx_on10+delay_onIdx_on30)==2)/sum(delay_onIdx_on10+delay_onIdx_on30)];
    delayRatio.overlap.off = [delayRatio.overlap.off,2*100*sum((delay_onIdx_off10+delay_onIdx_off30)==2)/sum(delay_onIdx_off10+delay_onIdx_off30)];
    delayRatio.overlap.length_10 = [delayRatio.overlap.length_10,2*100*sum((delay_onIdx_on10+delay_onIdx_off10)==2)/sum(delay_onIdx_on10+delay_onIdx_off10)];
    delayRatio.overlap.length_30 = [delayRatio.overlap.length_30,2*100*sum((delay_onIdx_on30+delay_onIdx_off30)==2)/sum(delay_onIdx_on30+delay_onIdx_off30)];
end

mazeActiveRatio.rat1079.on10_nanmean = nanmean(mazeRatio.on10);
mazeActiveRatio.rat1079.off10_nanmean = nanmean(mazeRatio.off10);
mazeActiveRatio.rat1079.on30_nanmean = nanmean(mazeRatio.on30);
mazeActiveRatio.rat1079.off30_nanmean = nanmean(mazeRatio.off30);
mazeActiveRatio.rat1079.All_nanmean = nanmean(mazeRatio.All);

mazeActiveRatio.rat1079.on10_std = std(mazeRatio.on10);
mazeActiveRatio.rat1079.off10_std = std(mazeRatio.off10);
mazeActiveRatio.rat1079.on30_std = std(mazeRatio.on30);
mazeActiveRatio.rat1079.off30_std = std(mazeRatio.off30);
mazeActiveRatio.rat1079.All_std = std(mazeRatio.All);

delayActiveRatio.rat1079.on10_nanmean = nanmean(delayRatio.on10);
delayActiveRatio.rat1079.off10_nanmean = nanmean(delayRatio.off10);
delayActiveRatio.rat1079.on30_nanmean = nanmean(delayRatio.on30);
delayActiveRatio.rat1079.off30_nanmean = nanmean(delayRatio.off30);
delayActiveRatio.rat1079.All_nanmean = nanmean(delayRatio.All);

delayActiveRatio.rat1079.on10_std = std(delayRatio.on10);
delayActiveRatio.rat1079.off10_std = std(delayRatio.off10);
delayActiveRatio.rat1079.on30_std = std(delayRatio.on30);
delayActiveRatio.rat1079.off30_std = std(delayRatio.off30);
delayActiveRatio.rat1079.All_std = std(delayRatio.All);

delayActiveRatio.rat1079.overlap.on_nanmean = nanmean(delayRatio.overlap.on);
delayActiveRatio.rat1079.overlap.off_nanmean = nanmean(delayRatio.overlap.off);
delayActiveRatio.rat1079.overlap.length_10_nanmean = nanmean(delayRatio.overlap.length_10);
delayActiveRatio.rat1079.overlap.length_30_nanmean = nanmean(delayRatio.overlap.length_30);

delayActiveRatio.rat1079.overlap.on_std = std(delayRatio.overlap.on);
delayActiveRatio.rat1079.overlap.off_std = std(delayRatio.overlap.off);
delayActiveRatio.rat1079.overlap.length_10_std = std(delayRatio.overlap.length_10);
delayActiveRatio.rat1079.overlap.length_30_std = std(delayRatio.overlap.length_30);



% rat all
AnalyzeSes = ratAll;
mazeRatio.on10 = []; mazeRatio.off10 = []; mazeRatio.on30 = []; mazeRatio.off30 = []; mazeRatio.All = [];
delayRatio.on10 = []; delayRatio.off10 = []; delayRatio.on30 = []; delayRatio.off30 = []; delayRatio.All = [];
delayRatio.overlap.on = []; delayRatio.overlap.off = []; delayRatio.overlap.length_10 = []; delayRatio.overlap.length_30 = [];
for i = AnalyzeSes(1:end)
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    rateCluster = SpikeProp.AvgRate.Fig8Rate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    for j = 1:length(sessDirs)
        
            maze_time.(sessDirs{j}) = [];
            delay_time.(sessDirs{j}) = [];
            block_time.(sessDirs{j}) = [];
            
            mazeSpike.(sessDirs{j}) = [];
            delaySpike.(sessDirs{j}) = [];
        
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);       
  
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
        %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
        %         delayTend2 = Fig8DelayZonePos.delayPos2.endT;
        
        trialNum = size(delayTstart1,2);
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
            delay_timeTemp = 10*trialNum;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
            delay_timeTemp = 30*trialNum;
        else
            error('Delay time is wrong')
        end
  
        block_timeTemp = pathData.t(end) - pathData.t(1);
        maze_timeTemp = pathData.t(end) - pathData.t(1) - delay_timeTemp;
        
         for k = 1:clusterNum             % get each spike time, change unit to msec from sec             % ts unit: ms             tSp = Spike_Session.(sessDirs{j}){k};             if rateLabel(k) == 1                 if ~isempty(tSp)                     % only consider pupative pyramidal cells                      ts_DelayCount = 0;                     for m = 1:trialNum                                                 % delay definition 1 start from barrier                         ts_DelayCount = ts_DelayCount + sum(tSp>delayTstart1(m) & tSp<delayTend1_2(m));                     end                                          delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),ts_DelayCount];                     mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),length(tSp)-ts_DelayCount];                     block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];                     maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];                     delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];                 else                     delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),0];                     mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),0];                     block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];                     maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];                     delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];                 end             end          end
            % get each spike time, change unit to msec from sec
            % ts unit: ms
            tSp = Spike_Session.(sessDirs{j}){k};
            if rateLabel(k) == 1
                if ~isempty(tSp)
                    % only consider pupative pyramidal cells

                    ts_DelayCount = 0;
                    for m = 1:trialNum                        
                        % delay definition 1 start from barrier
                        ts_DelayCount = ts_DelayCount + sum(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                    end
                    
                    delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),ts_DelayCount];
                    mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),length(tSp)-ts_DelayCount];
                    block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];
                    maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];
                    delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];
                else
                    delaySpike.(sessDirs{j}) = [delaySpike.(sessDirs{j}),0];
                    mazeSpike.(sessDirs{j}) = [mazeSpike.(sessDirs{j}),0];
                    block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),block_timeTemp];
                    maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),maze_timeTemp];
                    delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),delay_timeTemp];
                end
            end
         end 
    end
  
    % time used in 10 or 30 sec condition
    mazeRate.on10 = (mazeSpike.on10_1 + mazeSpike.on10_2)./(maze_time.on10_1 + maze_time.on10_2);
    mazeRate.off10 = (mazeSpike.off10_1 + mazeSpike.off10_2)./(maze_time.off10_1 + maze_time.off10_2);
    mazeRate.on30 = (mazeSpike.on30_1 + mazeSpike.on30_2)./(maze_time.on30_1 + maze_time.on30_2);
    mazeRate.off30 = (mazeSpike.off30_1 + mazeSpike.off30_2)./(maze_time.off30_1 + maze_time.off30_2);
    maze_onIdx_on10 = mazeRate.on10>p.avgRateThres;
    maze_onIdx_off10 = mazeRate.off10>p.avgRateThres;
    maze_onIdx_on30 = mazeRate.on30>p.avgRateThres;
    maze_onIdx_off30 = mazeRate.off30>p.avgRateThres;
    maze_onIdx = (maze_onIdx_on10+maze_onIdx_off10+maze_onIdx_on30+maze_onIdx_off30)>0;
    
    delayRate.on10 = (delaySpike.on10_1 + delaySpike.on10_2)./(delay_time.on10_1 + delay_time.on10_2);
    delayRate.off10 = (delaySpike.off10_1 + delaySpike.off10_2)./(delay_time.off10_1 + delay_time.off10_2);
    delayRate.on30 = (delaySpike.on30_1 + delaySpike.on30_2)./(delay_time.on30_1 + delay_time.on30_2);
    delayRate.off30 = (delaySpike.off30_1 + delaySpike.off30_2)./(delay_time.off30_1 + delay_time.off30_2);
    delay_onIdx_on10 = delayRate.on10>p.avgRateThres;
    delay_onIdx_off10 = delayRate.off10>p.avgRateThres;
    delay_onIdx_on30 = delayRate.on30>p.avgRateThres;
    delay_onIdx_off30 = delayRate.off30>p.avgRateThres;
    delay_onIdx = (delay_onIdx_on10+delay_onIdx_off10+delay_onIdx_on30+delay_onIdx_off30)>0;
    
    mazeRatio.on10 = [mazeRatio.on10,100*sum(maze_onIdx_on10)/length(maze_onIdx_on10)];
    mazeRatio.off10 = [mazeRatio.off10,100*sum(maze_onIdx_off10)/length(maze_onIdx_off10)];
    mazeRatio.on30 = [mazeRatio.on30,100*sum(maze_onIdx_on30)/length(maze_onIdx_on30)];
    mazeRatio.off30 = [mazeRatio.off30,100*sum(maze_onIdx_off30)/length(maze_onIdx_off30)];
    mazeRatio.All = [mazeRatio.All,100*sum(maze_onIdx)/length(maze_onIdx)];
    
    delayRatio.on10 = [delayRatio.on10,100*sum(delay_onIdx_on10)/length(delay_onIdx_on10)];
    delayRatio.off10 = [delayRatio.off10,100*sum(delay_onIdx_off10)/length(delay_onIdx_off10)];
    delayRatio.on30 = [delayRatio.on30,100*sum(delay_onIdx_on30)/length(delay_onIdx_on30)];
    delayRatio.off30 = [delayRatio.off30,100*sum(delay_onIdx_off30)/length(delay_onIdx_off30)];
    delayRatio.All = [delayRatio.All,100*sum(delay_onIdx)/length(delay_onIdx)];
    
    delayRatio.overlap.on = [delayRatio.overlap.on,2*100*sum((delay_onIdx_on10+delay_onIdx_on30)==2)/sum(delay_onIdx_on10+delay_onIdx_on30)];
    delayRatio.overlap.off = [delayRatio.overlap.off,2*100*sum((delay_onIdx_off10+delay_onIdx_off30)==2)/sum(delay_onIdx_off10+delay_onIdx_off30)];
    delayRatio.overlap.length_10 = [delayRatio.overlap.length_10,2*100*sum((delay_onIdx_on10+delay_onIdx_off10)==2)/sum(delay_onIdx_on10+delay_onIdx_off10)];
    delayRatio.overlap.length_30 = [delayRatio.overlap.length_30,2*100*sum((delay_onIdx_on30+delay_onIdx_off30)==2)/sum(delay_onIdx_on30+delay_onIdx_off30)];
end

mazeActiveRatio.ratAll.on10_nanmean = nanmean(mazeRatio.on10);
mazeActiveRatio.ratAll.off10_nanmean = nanmean(mazeRatio.off10);
mazeActiveRatio.ratAll.on30_nanmean = nanmean(mazeRatio.on30);
mazeActiveRatio.ratAll.off30_nanmean = nanmean(mazeRatio.off30);
mazeActiveRatio.ratAll.All_nanmean = nanmean(mazeRatio.All);

mazeActiveRatio.ratAll.on10_std = std(mazeRatio.on10);
mazeActiveRatio.ratAll.off10_std = std(mazeRatio.off10);
mazeActiveRatio.ratAll.on30_std = std(mazeRatio.on30);
mazeActiveRatio.ratAll.off30_std = std(mazeRatio.off30);
mazeActiveRatio.ratAll.All_std = std(mazeRatio.All);

delayActiveRatio.ratAll.on10_nanmean = nanmean(delayRatio.on10);
delayActiveRatio.ratAll.off10_nanmean = nanmean(delayRatio.off10);
delayActiveRatio.ratAll.on30_nanmean = nanmean(delayRatio.on30);
delayActiveRatio.ratAll.off30_nanmean = nanmean(delayRatio.off30);
delayActiveRatio.ratAll.All_nanmean = nanmean(delayRatio.All);

delayActiveRatio.ratAll.on10_std = std(delayRatio.on10);
delayActiveRatio.ratAll.off10_std = std(delayRatio.off10);
delayActiveRatio.ratAll.on30_std = std(delayRatio.on30);
delayActiveRatio.ratAll.off30_std = std(delayRatio.off30);
delayActiveRatio.ratAll.All_std = std(delayRatio.All);

delayActiveRatio.ratAll.overlap.on_nanmean = nanmean(delayRatio.overlap.on);
delayActiveRatio.ratAll.overlap.off_nanmean = nanmean(delayRatio.overlap.off);
delayActiveRatio.ratAll.overlap.length_10_nanmean = nanmean(delayRatio.overlap.length_10);
delayActiveRatio.ratAll.overlap.length_30_nanmean = nanmean(delayRatio.overlap.length_30);

delayActiveRatio.ratAll.overlap.on_std = std(delayRatio.overlap.on);
delayActiveRatio.ratAll.overlap.off_std = std(delayRatio.overlap.off);
delayActiveRatio.ratAll.overlap.length_10_std = std(delayRatio.overlap.length_10);
delayActiveRatio.ratAll.overlap.length_30_std = std(delayRatio.overlap.length_30);




figure(1)
bar(1,100*sum(mazeRate.on10>0.5)/length(mazeRate.on10));
hold on
bar(2,100*sum(mazeRate.off10>0.5)/length(mazeRate.on10));
bar(3,100*sum(mazeRate.on30>0.5)/length(mazeRate.on10));
bar(4,100*sum(mazeRate.off30>0.5)/length(mazeRate.on10));

bar(6,100*sum(delayRate.on10>0.5)/length(delayRate.on10));
bar(7,100*sum(delayRate.off10>0.5)/length(delayRate.on10));
bar(8,100*sum(delayRate.on30>0.5)/length(delayRate.on10));
bar(9,100*sum(delayRate.off30>0.5)/length(delayRate.on10));

xlim([0 10])
ylim([0 100])
title('Maze active ratio 0.5 Hz')
set(gca, 'XTick', [1,6], 'XTickLabel', {'maze','delay'});

figure(2)
bar(1,100*sum(mazeRate.on10>1)/length(mazeRate.on10));
hold on
bar(2,100*sum(mazeRate.off10>1)/length(mazeRate.on10));
bar(3,100*sum(mazeRate.on30>1)/length(mazeRate.on10));
bar(4,100*sum(mazeRate.off30>1)/length(mazeRate.on10));

bar(6,100*sum(delayRate.on10>1)/length(delayRate.on10));
bar(7,100*sum(delayRate.off10>1)/length(delayRate.on10));
bar(8,100*sum(delayRate.on30>1)/length(delayRate.on10));
bar(9,100*sum(delayRate.off30>1)/length(delayRate.on10));

xlim([0 10])
ylim([0 100])
title('Maze active ratio 1 Hz')
set(gca, 'XTick', [1,6], 'XTickLabel', {'maze','delay'});
