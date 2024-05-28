% Plot time % being used in each delay conditions
function Fig8Treadmill_MazeRate_Quant(inFile,AnalyzeSes)
close all

p.savePlot = 0;
p.writeToFile = 0;
% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'

% Read in input information
sessInfo = SessInfoImport(inFile);

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
        
        if i == AnalyzeSes(1)
            maze_time.(sessDirs{j}) = [];
            delay_time.(sessDirs{j}) = [];
            block_time.(sessDirs{j}) = [];
            
            mazeSpike.(sessDirs{j}) = [];
            delaySpike.(sessDirs{j}) = [];
        end
        
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
        
         for k = 1:clusterNum
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
    fprintf('Finished position analysis for session %d\n',i);
end

% time used in 10 or 30 sec condition
mazeRate.on10 = (mazeSpike.on10_1 + mazeSpike.on10_2)./(maze_time.on10_1 + maze_time.on10_2);
mazeRate.off10 = (mazeSpike.off10_1 + mazeSpike.off10_2)./(maze_time.off10_1 + maze_time.off10_2);
mazeRate.on30 = (mazeSpike.on30_1 + mazeSpike.on30_2)./(maze_time.on30_1 + maze_time.on30_2);
mazeRate.off30 = (mazeSpike.off30_1 + mazeSpike.off30_2)./(maze_time.off30_1 + maze_time.off30_2);

delayRate.on10 = (delaySpike.on10_1 + delaySpike.on10_2)./(delay_time.on10_1 + delay_time.on10_2);
delayRate.off10 = (delaySpike.off10_1 + delaySpike.off10_2)./(delay_time.off10_1 + delay_time.off10_2);
delayRate.on30 = (delaySpike.on30_1 + delaySpike.on30_2)./(delay_time.on30_1 + delay_time.on30_2);
delayRate.off30 = (delaySpike.off30_1 + delaySpike.off30_2)./(delay_time.off30_1 + delay_time.off30_2);

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

end