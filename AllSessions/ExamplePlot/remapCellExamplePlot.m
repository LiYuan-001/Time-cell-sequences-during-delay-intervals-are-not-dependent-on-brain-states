
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
for i = 4
    
    
    % load delayFire map
    delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap_2Session.mat');
    load(delayFile);
        
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
%     [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
%     clusterNum = length(TList);

    TList = Fig8DelayTimeMap_2Session.tList;
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    
    % get each phase names (no delay etc)
    sessDirs = {'on10','off10','on30','off30'};
    
    h1 = figure;
    h1.Position = [100 100 900 300];
    
            
    for j = 1:4
       
        subplot(1,4,j);
        % tt5_04;        
        for k = 6
            
            timeMap_Def1 = Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate1_Smooth{k};
            pre_Map_Def1 = Fig8DelayTimeMap_2Session.(sessDirs{j}).pre_spikeRate1_Smooth{k};
            post_Map_Def1 = Fig8DelayTimeMap_2Session.(sessDirs{j}).post_spikeRate1_2session{k}(:,1:20);
            
            imagesc([pre_Map_Def1,timeMap_Def1,post_Map_Def1])
            colormap(jet)
            xlabel(sprintf('%s%1.2f','Mean rate = ',mean(mean(timeMap_Def1))));
            if j == 1
                TITLE1 = sprintf('%s%d%s%d%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k});
                TITLE2 = sessDirs{j};
                title({TITLE1;TITLE2},'Interpreter','None')
            else
                TITLE1 = [];
                TITLE2 = sessDirs{j};
                title({TITLE1;TITLE2},'Interpreter','None')
            end
        end               
    end
    
    clear Fig8DelayTimeMap_2Session
end

for i = 6
    
    
    % load delayFire map
    delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap_2Session.mat');
    load(delayFile);
        
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
%     [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
%     clusterNum = length(TList);

    TList = Fig8DelayTimeMap_2Session.tList;
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    
    % get each phase names (no delay etc)
    sessDirs = {'on10','off10','on30','off30'};
    
    h1 = figure;
    h1.Position = [100 100 900 300];
    
            
    for j = 1:4
       
        subplot(1,4,j);
        % tt5_04;        
        for k = 24
            
            timeMap_Def1 = Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate1_Smooth{k};
            pre_Map_Def1 = Fig8DelayTimeMap_2Session.(sessDirs{j}).pre_spikeRate1_Smooth{k};
            post_Map_Def1 = Fig8DelayTimeMap_2Session.(sessDirs{j}).post_spikeRate1_2session{k}(:,1:20);
            
            imagesc([pre_Map_Def1,timeMap_Def1,post_Map_Def1])
            colormap(jet)
            xlabel(sprintf('%s%1.2f','Mean rate = ',mean(mean(timeMap_Def1))));
            if j == 1
                TITLE1 = sprintf('%s%d%s%d%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k});
                TITLE2 = sessDirs{j};
                title({TITLE1;TITLE2},'Interpreter','None')
            else
                TITLE1 = [];
                TITLE2 = sessDirs{j};
                title({TITLE1;TITLE2},'Interpreter','None')
            end
        end               
    end
    
    clear Fig8DelayTimeMap_2Session
end

for i = 9
    
    
    % load delayFire map
    delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap_2Session.mat');
    load(delayFile);
        
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
%     [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
%     clusterNum = length(TList);

    TList = Fig8DelayTimeMap_2Session.tList;
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    
    % get each phase names (no delay etc)
    sessDirs = {'on10','off10','on30','off30'};
    
    h1 = figure;
    h1.Position = [100 100 900 300];
    
            
    for j = 1:4
       
        subplot(1,4,j);
        % tt5_04;        
        for k = 10
            
            timeMap_Def1 = Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate1_Smooth{k};
            pre_Map_Def1 = Fig8DelayTimeMap_2Session.(sessDirs{j}).pre_spikeRate1_Smooth{k};
            post_Map_Def1 = Fig8DelayTimeMap_2Session.(sessDirs{j}).post_spikeRate1_2session{k}(:,1:20);
            
            imagesc([pre_Map_Def1,timeMap_Def1,post_Map_Def1])
            colormap(jet)
            xlabel(sprintf('%s%1.2f','Mean rate = ',mean(mean(timeMap_Def1))));
            if j == 1
                TITLE1 = sprintf('%s%d%s%d%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k});
                TITLE2 = sessDirs{j};
                title({TITLE1;TITLE2},'Interpreter','None')
            else
                TITLE1 = [];
                TITLE2 = sessDirs{j};
                title({TITLE1;TITLE2},'Interpreter','None')
            end
        end               
    end
    
    clear Fig8DelayTimeMap_2Session
end