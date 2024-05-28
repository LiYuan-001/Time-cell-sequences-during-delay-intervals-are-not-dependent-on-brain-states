
p.savePlot = 0;
p.writeToFile = 0;

% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'

p.timeBin = 150; % ms;
p.gaussSigma = p.timeBin*2; % sd of the gaussian distribution


% Read in input information
sessInfo = SessInfoImport('V:\LiYuan\Codes\Fig8MazeTreadmill_V2\Fig8Treadmill_OnOff.xlsx');
close all

% use example from 1046-day5
for i = 10
    
    if p.writeToFile
        savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
        if ~exist(savedir2, 'dir')
            mkdir(savedir2);
        end
    end
    
    delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap_2Session.mat');
    load(delayFile); 
    stabilityFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayFireStability_TrialbyTrial_2Session.mat');
    load(stabilityFile);
    
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);
   
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    h1 = figure(1);
    h1.Position = [100 100 900 600];
    h2 = figure(2);
    h1.Position = [100 100 900 600];
    % define where each subplot should be
%     pos1 = [0.1 0.35 0.25 0.6];
%     pos2 = [0.4 0.35 0.25 0.6];
%     pos3 = [0.7 0.35 0.25 0.6];
%     pos4 = [0.4 0.15 0.25 0.1];
            
    for j = 1
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        
 
        % -----------------------------------------------------------------
        % use the mode of delay is always 10 sec or 30 sec
        % rather than identify real time in delay zone
        % to make plot cleaner
        % I can change to delay time in delay zone later
        % Li Yuan, 19-Aug-2020, UCSD
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
            
        trialNum = size(delayTstart1,2);
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT; 
        else
            error('Delay time is wrong')
        end       
                

        for k = 4
            spikeRate_cell_1 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).spikeRate1_Smooth{k};
            spikeRate_Combined_cell_1 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).spikeRate1_Combined_Smooth{k};
            
            pre_spikeRate_cell_1 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).pre_spikeRate1_Smooth{k};
            pre_spikeRate_Combined_cell_1 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).pre_spikeRate1_Combined_Smooth{k};
            
            post_spikeRate_cell_1 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).post_spikeRate1_2session{k}(:,1:20);
            post_spikeRate_Combined_cell_1 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).post_mapMean_Def1_2Session{k}(:,1:20);
            
            
            % get each spike time
            tSp = Spike_Session.(sessDirs{j}){k};
            trialNum_1 = length(delayTstart1);
            
            for m = 1:trialNum_1
                
                % delay definition 1
                spikeTime1 = tSp-delayTstart1(m);
                ts_Delay1 = tSp(tSp>(delayTstart1(m)-3) & tSp<delayTend1_2(m));
                ts_Delay2 = tSp(tSp> delayTend1(m) & tSp<(delayTend1(m)+3)) - delayTend1(m);
                ts_Delay1 = [ts_Delay1-delayTstart1(m);ts_Delay2+maxT];
                
                xPoints = [ts_Delay1';ts_Delay1'];
                yPoints = [m+zeros(size(ts_Delay1'))-0.3;m+zeros(size(ts_Delay1'))+0.3];
       
                figure(1)
                subplot(2,3,1)
                if ~isempty(ts_Delay1)
                    plot(xPoints,yPoints,'k')
                else
                    plot(0,m)
                end
                hold on
            end
        end    
        
        for k = 22
            spikeRate_cell_2 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).spikeRate1_Smooth{k};
            spikeRate_Combined_cell_2 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).spikeRate1_Combined_Smooth{k};
            
            pre_spikeRate_cell_2 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).pre_spikeRate1_Smooth{k};
            pre_spikeRate_Combined_cell_2 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).pre_spikeRate1_Combined_Smooth{k};
            
            post_spikeRate_cell_2 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).post_spikeRate1_2session{k}(:,1:20);
            post_spikeRate_Combined_cell_2 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).post_mapMean_Def1_2Session{k}(:,1:20);
            
            % get each spike time
            tSp = Spike_Session.(sessDirs{j}){k};
            trialNum_1 = length(delayTstart1);
            
            for m = 1:trialNum_1
                % delay definition 1
                spikeTime1 = tSp-delayTstart1(m);
                ts_Delay1 = tSp(tSp>(delayTstart1(m)-3) & tSp<delayTend1_2(m));
                ts_Delay2 = tSp(tSp> delayTend1(m) & tSp<(delayTend1(m)+3)) - delayTend1(m);
                ts_Delay1 = [ts_Delay1-delayTstart1(m);ts_Delay2+maxT];
                
                xPoints = [ts_Delay1';ts_Delay1'];
                yPoints = [m+zeros(size(ts_Delay1'))-0.3;m+zeros(size(ts_Delay1'))+0.3];
       
                figure(2)
                subplot(2,3,1)
                if ~isempty(ts_Delay1)
                    plot(xPoints,yPoints,'k')
                else
                    plot(0,m)
                end
                hold on
            end
        end               
    end
    
    for j = 5
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        
 
        % -----------------------------------------------------------------
        % use the mode of delay is always 10 sec or 30 sec
        % rather than identify real time in delay zone
        % to make plot cleaner
        % I can change to delay time in delay zone later
        % Li Yuan, 19-Aug-2020, UCSD
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
            
        trialNum = size(delayTstart1,2);
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT; 
        else
            error('Delay time is wrong')
        end       
                

        for k = 4
            
            % get each spike time
            tSp = Spike_Session.(sessDirs{j}){k};
            trialNum_2 = length(delayTstart1);
            
            for m = 1:trialNum_2
                spikeTime1 = tSp-delayTstart1(m);
                ts_Delay1 = tSp(tSp>(delayTstart1(m)-3) & tSp<delayTend1_2(m));
                ts_Delay2 = tSp(tSp> delayTend1(m) & tSp<(delayTend1(m)+3)) - delayTend1(m);
                ts_Delay1 = [ts_Delay1-delayTstart1(m);ts_Delay2+maxT];
                
                xPoints = [ts_Delay1';ts_Delay1'];
                yPoints = [trialNum_1+m+zeros(size(ts_Delay1'))-0.3;trialNum_1+m+zeros(size(ts_Delay1'))+0.3];
       
                figure(1)
                subplot(2,3,1)
                if ~isempty(ts_Delay1)
                    plot(xPoints,yPoints,'k')
                else
                    plot(0,m)
                end
                hold on
            end
            subplot(2,3,3)
            shuffleVal = DelayFireStability_TrialbyTrial_2Session.on10.delayCorr_Def1_Trial_Shuffle(:,k);
            sigVal = DelayFireStability_TrialbyTrial_2Session.on10.delayCorr_Def2_Trial_Shuffle95(k);
            originVal = DelayFireStability_TrialbyTrial_2Session.on10.delayCorr_Def2_Trial(k);
            
            Violin(shuffleVal,1,'ShowData',false);
            hold on
            plot([0.7,1.3],ones(1,2).*sigVal,'r-')
            plot(1,originVal,'r*')
            hold off
        end    
           
        for k = 22
           
            % get each spike time
            tSp = Spike_Session.(sessDirs{j}){k};
            trialNum_2 = length(delayTstart1);
            
            for m = 1:trialNum_2
                spikeTime1 = tSp-delayTstart1(m);
                ts_Delay1 = tSp(tSp>(delayTstart1(m)-3) & tSp<delayTend1_2(m));
                ts_Delay2 = tSp(tSp> delayTend1(m) & tSp<(delayTend1(m)+3)) - delayTend1(m);
                ts_Delay1 = [ts_Delay1-delayTstart1(m);ts_Delay2+maxT];
                
                xPoints = [ts_Delay1';ts_Delay1'];
                yPoints = [trialNum_1+m+zeros(size(ts_Delay1'))-0.3;trialNum_1+m+zeros(size(ts_Delay1'))+0.3];
                
                figure(2)
                subplot(2,3,1)
                if ~isempty(ts_Delay1)
                    plot(xPoints,yPoints,'k')
                else
                    plot(0,m)
                end
                hold on
            end
            subplot(2,3,3)
            shuffleVal = DelayFireStability_TrialbyTrial_2Session.on10.delayCorr_Def1_Trial_Shuffle(:,k);
            sigVal = DelayFireStability_TrialbyTrial_2Session.on10.delayCorr_Def2_Trial_Shuffle95(k);
            originVal = DelayFireStability_TrialbyTrial_2Session.on10.delayCorr_Def2_Trial(k);
            
            Violin(shuffleVal,1,'ShowData',false);
            hold on
            plot([0.7,1.3],ones(1,2).*sigVal,'r-')
            plot(1,originVal,'r*')
            hold off
        end       
    end
    
    figure(1)
    subplot(2,3,1)
    TITLE1 = sprintf('%s%d%s%d%s%s%s%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k}(1:end-2),'-',sessDirs{j},'-Spikes');
    title({TITLE1},'Interpreter','None')
    xlabel('Time (s)')
    ylabel('Trials')
    ylim([0.5 trialNum_1+trialNum_2+0.5])
    xlim([-3 maxT+3])
    xTick= [-3 0 maxT maxT+3];
    xTickLabel = [-3 0 maxT maxT+3];
    set (gca,'YDir','reverse')
    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
    
    subplot(2,3,4)
    preRateLine = [pre_spikeRate_Combined_cell_1];
    RateLine = [spikeRate_Combined_cell_1];
    postRateLine = [post_spikeRate_Combined_cell_1];
    plot(movmean([preRateLine,RateLine,postRateLine],10))
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    xlim([0 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)])
    xTick= [0 20 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)];
    xTickLabel = [-3 0 maxT+3];
    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
    
    
    subplot(2,3,2)
    rateMap = [pre_spikeRate_cell_1,spikeRate_cell_1,post_spikeRate_cell_1];
    imagesc(rateMap)
    colormap(jet)
    TITLE1 = 'RateMap - Onset at barrier';    
    title({TITLE1},'Interpreter','None')
    xlabel('Time (s)')
    ylabel('Trials')
    xlim([0 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)])
    xTick= [0 20 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)];
    xTickLabel = [-3 0 maxT+3];
    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
    
    figure(2)
    subplot(2,3,1)
    TITLE1 = sprintf('%s%d%s%d%s%s%s%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k}(1:end-2),'-',sessDirs{j},'-Spikes');
    title({TITLE1},'Interpreter','None')
    xlabel('Time (Sec)')
    ylabel('Trials')
    ylim([0.5 trialNum_1+trialNum_2+0.5])
    xlim([-3 maxT+3])
    xTick= [-3 0 maxT maxT+3];
    xTickLabel = [-3 0 maxT maxT+3];
    set (gca,'YDir','reverse')
    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
    
    subplot(2,3,4)
    preRateLine = [pre_spikeRate_Combined_cell_2];
    RateLine = [spikeRate_Combined_cell_2];
    postRateLine = [post_spikeRate_Combined_cell_2];
    plot(movmean([preRateLine,RateLine,postRateLine],10))
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    xlim([0 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)])
    xTick= [0 20 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)];
    xTickLabel = [-3 0 maxT+3];
    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
    
    
    subplot(2,3,2)
    rateMap = [pre_spikeRate_cell_2,spikeRate_cell_2,post_spikeRate_cell_2];
    imagesc(rateMap)
    colormap(jet)
    TITLE1 = 'RateMap - Onset at barrier';    
    title({TITLE1},'Interpreter','None')
    xlabel('Time (s)')
    ylabel('Trials')
    xlim([0 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)])
    xTick= [0 20 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)];
    xTickLabel = [-3 0 maxT+3];
    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
    
    
    for j = 3
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        
 
        % -----------------------------------------------------------------
        % use the mode of delay is always 10 sec or 30 sec
        % rather than identify real time in delay zone
        % to make plot cleaner
        % I can change to delay time in delay zone later
        % Li Yuan, 19-Aug-2020, UCSD
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
            
        trialNum = size(delayTstart1,2);
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT; 
        else
            error('Delay time is wrong')
        end       
                

        for k = 4
            spikeRate_cell_1 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).spikeRate1_Smooth{k};
            spikeRate_Combined_cell_1 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).spikeRate1_Combined_Smooth{k};
            
            pre_spikeRate_cell_1 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).pre_spikeRate1_Smooth{k};
            pre_spikeRate_Combined_cell_1 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).pre_spikeRate1_Combined_Smooth{k};
            
            post_spikeRate_cell_1 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).post_spikeRate1_2session{k}(:,1:20);
            post_spikeRate_Combined_cell_1 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).post_mapMean_Def1_2Session{k}(:,1:20);
            
            
            % get each spike time
            tSp = Spike_Session.(sessDirs{j}){k};
            trialNum_1 = length(delayTstart1);
            
            for m = 1:trialNum_1
                
                % delay definition 1
                spikeTime1 = tSp-delayTstart1(m);
                ts_Delay1 = tSp(tSp>(delayTstart1(m)-3) & tSp<delayTend1_2(m));
                ts_Delay2 = tSp(tSp> delayTend1(m) & tSp<(delayTend1(m)+3)) - delayTend1(m);
                ts_Delay1 = [ts_Delay1-delayTstart1(m);ts_Delay2+maxT];
                
                xPoints = [ts_Delay1';ts_Delay1'];
                yPoints = [m+zeros(size(ts_Delay1'))-0.3;m+zeros(size(ts_Delay1'))+0.3];
       
                figure(3)
                subplot(2,3,1)
                if ~isempty(ts_Delay1)
                    plot(xPoints,yPoints,'k')
                else
                    plot(0,m)
                end
                hold on
            end
        end    
        
        for k = 22
            spikeRate_cell_2 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).spikeRate1_Smooth{k};
            spikeRate_Combined_cell_2 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).spikeRate1_Combined_Smooth{k};
            
            pre_spikeRate_cell_2 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).pre_spikeRate1_Smooth{k};
            pre_spikeRate_Combined_cell_2 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).pre_spikeRate1_Combined_Smooth{k};
            
            post_spikeRate_cell_2 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).post_spikeRate1_2session{k}(:,1:20);
            post_spikeRate_Combined_cell_2 = Fig8DelayTimeMap_2Session.(sessDirs{j}(1:end-2)).post_mapMean_Def1_2Session{k}(:,1:20);
            
            % get each spike time
            tSp = Spike_Session.(sessDirs{j}){k};
            trialNum_1 = length(delayTstart1);
            
            for m = 1:trialNum_1
                % delay definition 1
                spikeTime1 = tSp-delayTstart1(m);
                ts_Delay1 = tSp(tSp>(delayTstart1(m)-3) & tSp<delayTend1_2(m));
                ts_Delay2 = tSp(tSp> delayTend1(m) & tSp<(delayTend1(m)+3)) - delayTend1(m);
                ts_Delay1 = [ts_Delay1-delayTstart1(m);ts_Delay2+maxT];
                
                xPoints = [ts_Delay1';ts_Delay1'];
                yPoints = [m+zeros(size(ts_Delay1'))-0.3;m+zeros(size(ts_Delay1'))+0.3];
       
                figure(4)
                subplot(2,3,1)
                if ~isempty(ts_Delay1)
                    plot(xPoints,yPoints,'k')
                else
                    plot(0,m)
                end
                hold on
            end
        end 
    end
    
    for j = 7
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        
 
        % -----------------------------------------------------------------
        % use the mode of delay is always 10 sec or 30 sec
        % rather than identify real time in delay zone
        % to make plot cleaner
        % I can change to delay time in delay zone later
        % Li Yuan, 19-Aug-2020, UCSD
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
   
        trialNum = size(delayTstart1,2);
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT; 
        else
            error('Delay time is wrong')
        end       
                

        for k = 4
            
            % get each spike time
            tSp = Spike_Session.(sessDirs{j}){k};
            trialNum_2 = length(delayTstart1);
            
            for m = 1:trialNum_2
                spikeTime1 = tSp-delayTstart1(m);
                ts_Delay1 = tSp(tSp>(delayTstart1(m)-3) & tSp<delayTend1_2(m));
                ts_Delay2 = tSp(tSp> delayTend1(m) & tSp<(delayTend1(m)+3)) - delayTend1(m);
                ts_Delay1 = [ts_Delay1-delayTstart1(m);ts_Delay2+maxT];
                
                xPoints = [ts_Delay1';ts_Delay1'];
                yPoints = [trialNum_1+m+zeros(size(ts_Delay1'))-0.3;trialNum_1+m+zeros(size(ts_Delay1'))+0.3];
       
                figure(3)
                subplot(2,3,1)
                if ~isempty(ts_Delay1)
                    plot(xPoints,yPoints,'k')
                else
                    plot(0,m)
                end
                hold on
            end
            subplot(2,3,3)
            shuffleVal = DelayFireStability_TrialbyTrial_2Session.on10.delayCorr_Def1_Trial_Shuffle(:,k);
            sigVal = DelayFireStability_TrialbyTrial_2Session.on10.delayCorr_Def2_Trial_Shuffle95(k);
            originVal = DelayFireStability_TrialbyTrial_2Session.on10.delayCorr_Def2_Trial(k);
            
            Violin(shuffleVal,1,'ShowData',false);
            hold on
            plot([0.7,1.3],ones(1,2).*sigVal,'r-')
            plot(1,originVal,'r*')
            hold off
        end    
           
        for k = 22
           
            % get each spike time
            tSp = Spike_Session.(sessDirs{j}){k};
            trialNum_2 = length(delayTstart1);
            
            for m = 1:trialNum_2
                spikeTime1 = tSp-delayTstart1(m);
                ts_Delay1 = tSp(tSp>(delayTstart1(m)-3) & tSp<delayTend1_2(m));
                ts_Delay2 = tSp(tSp> delayTend1(m) & tSp<(delayTend1(m)+3)) - delayTend1(m);
                ts_Delay1 = [ts_Delay1-delayTstart1(m);ts_Delay2+maxT];
                
                xPoints = [ts_Delay1';ts_Delay1'];
                yPoints = [trialNum_1+m+zeros(size(ts_Delay1'))-0.3;trialNum_1+m+zeros(size(ts_Delay1'))+0.3];
                
                figure(4)
                subplot(2,3,1)
                if ~isempty(ts_Delay1)
                    plot(xPoints,yPoints,'k')
                else
                    plot(0,m)
                end
                hold on
            end
            subplot(2,3,3)
            shuffleVal = DelayFireStability_TrialbyTrial_2Session.on10.delayCorr_Def1_Trial_Shuffle(:,k);
            sigVal = DelayFireStability_TrialbyTrial_2Session.on10.delayCorr_Def2_Trial_Shuffle95(k);
            originVal = DelayFireStability_TrialbyTrial_2Session.on10.delayCorr_Def2_Trial(k);
            
            Violin(shuffleVal,1,'ShowData',false);
            hold on
            plot([0.7,1.3],ones(1,2).*sigVal,'r-')
            plot(1,originVal,'r*')
            hold off
        end       
    end
    
    figure(3)
    subplot(2,3,1)
    TITLE1 = sprintf('%s%d%s%d%s%s%s%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k}(1:end-2),'-',sessDirs{j},'-Spikes');
    title({TITLE1},'Interpreter','None')
    xlabel('Time (s)')
    ylabel('Trials')
    ylim([0.5 trialNum_1+trialNum_2+0.5])
    xlim([-3 maxT+3])
    xTick= [-3 0 maxT maxT+3];
    xTickLabel = [-3 0 maxT maxT+3];
    set (gca,'YDir','reverse')
    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
    
    subplot(2,3,4)
    preRateLine = [pre_spikeRate_Combined_cell_1];
    RateLine = [spikeRate_Combined_cell_1];
    postRateLine = [post_spikeRate_Combined_cell_1];
    plot(movmean([preRateLine,RateLine,postRateLine],10))
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    xlim([0 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)])
    xTick= [0 20 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)];
    xTickLabel = [-3 0 maxT+3];
    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
    
    
    subplot(2,3,2)
    rateMap = [pre_spikeRate_cell_1,spikeRate_cell_1,post_spikeRate_cell_1];
    imagesc(rateMap)
    colormap(jet)
    TITLE1 = 'RateMap - Onset at barrier';    
    title({TITLE1},'Interpreter','None')
    xlabel('Time (s)')
    ylabel('Trials')
    xlim([0 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)])
    xTick= [0 20 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)];
    xTickLabel = [-3 0 maxT+3];
    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
    
    figure(4)
    subplot(2,3,1)
    TITLE1 = sprintf('%s%d%s%d%s%s%s%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k}(1:end-2),'-',sessDirs{j},'-Spikes');
    title({TITLE1},'Interpreter','None')
    xlabel('Time (Sec)')
    ylabel('Trials')
    ylim([0.5 trialNum_1+trialNum_2+0.5])
    xlim([-3 maxT+3])
    xTick= [-3 0 maxT maxT+3];
    xTickLabel = [-3 0 maxT maxT+3];
    set (gca,'YDir','reverse')
    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
    
    subplot(2,3,4)
    preRateLine = [pre_spikeRate_Combined_cell_2];
    RateLine = [spikeRate_Combined_cell_2];
    postRateLine = [post_spikeRate_Combined_cell_2];
    plot(movmean([preRateLine,RateLine,postRateLine],10))
    xlabel('Time (s)')
    ylabel('Rate (Hz)')
    xlim([0 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)])
    xTick= [0 20 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)];
    xTickLabel = [-3 0 maxT+3];
    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
    
    
    subplot(2,3,2)
    rateMap = [pre_spikeRate_cell_2,spikeRate_cell_2,post_spikeRate_cell_2];
    imagesc(rateMap)
    colormap(jet)
    TITLE1 = 'RateMap - Onset at barrier';    
    title({TITLE1},'Interpreter','None')
    xlabel('Time (s)')
    ylabel('Trials')
    xlim([0 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)])
    xTick= [0 20 size(pre_spikeRate_cell_1,2)+size(spikeRate_cell_1,2)+size(post_spikeRate_cell_1,2)];
    xTickLabel = [-3 0 maxT+3];
    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
    

    clear DelayFire
    fprintf('Finished position analysis for session %d\n',i);
end

