% Plot time % being used in each delay conditions
function Fig8Treadmill_MazeTime_Quant(inFile,AnalyzeSes)
close all

p.savePlot = 0;
p.writeToFile = 0;


% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    
    for j = 1:length(sessDirs)
        
        if i == AnalyzeSes(1)
            maze_time.(sessDirs{j}) = [];
            delay_time.(sessDirs{j}) = [];
            block_time.(sessDirs{j}) = [];
        end
        
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        % load analyzed map per trial
        load(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'));        

        trialNum = length(ratesByECLR.ECLR);
        if contains(sessDirs{j},'10')
            delayTemp = 10*trialNum;
            delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),10*trialNum];
        else
            delayTemp = 30*trialNum;
            delay_time.(sessDirs{j}) = [delay_time.(sessDirs{j}),30*trialNum];
        end
        
        block_time.(sessDirs{j}) = [block_time.(sessDirs{j}),pathData.t(end) - pathData.t(1)];
        maze_time.(sessDirs{j}) = [maze_time.(sessDirs{j}),pathData.t(end) - pathData.t(1) - delayTemp];
    end
    fprintf('Finished position analysis for session %d\n',i);
end

% time used in 10 or 30 sec condition
maze.block10 = 100*(maze_time.on10_1+maze_time.on10_2+maze_time.off10_1+maze_time.off10_2)./...
    (block_time.on10_1+block_time.on10_2+block_time.off10_1+block_time.off10_2);
maze.block30 = 100*(maze_time.on30_1+maze_time.on30_2+maze_time.off30_1+maze_time.off30_2)./...
    (block_time.on30_1+block_time.on30_2+block_time.off30_1+block_time.off30_2);
delay.block10 = 100*(delay_time.on10_1+delay_time.on10_2+delay_time.off10_1+delay_time.off10_2)./...
    (block_time.on10_1+block_time.on10_2+block_time.off10_1+block_time.off10_2);
delay.block30 = 100*(delay_time.on30_1+delay_time.on30_2+delay_time.off30_1+delay_time.off30_2)./...
    (block_time.on30_1+block_time.on30_2+block_time.off30_1+block_time.off30_2);

% time used in 10 or 30 sec condition abslute time
maze_abs.block10 = (maze_time.on10_1+maze_time.on10_2+maze_time.off10_1+maze_time.off10_2);
maze_abs.block30 = (maze_time.on30_1+maze_time.on30_2+maze_time.off30_1+maze_time.off30_2);
delay_abs.block10 = (delay_time.on10_1+delay_time.on10_2+delay_time.off10_1+delay_time.off10_2);
delay_abs.block30 = (delay_time.on30_1+delay_time.on30_2+delay_time.off30_1+delay_time.off30_2);

% time %
figure(1)
h1 = subplot(1,2,1);
plot([maze.block10;delay.block10],'d-','Color',[0.8,0.8,0.8])
hold on
boxplot(maze.block10,'Position',1);
boxplot(delay.block10,'Position',2);
xlim([0 3])
ylim([0 100])
title('10 sec maze time')
set(gca, 'XTick', [1,2], 'XTickLabel', {'maze','delay'});

h2 = subplot(1,2,2);
plot([maze.block30;delay.block30],'d-','Color',[0.8,0.8,0.8])
hold on
boxplot(maze.block30,'Position',1);
boxplot(delay.block30,'Position',2);
xlim([0 3])
ylim([0 100])
title('10 sec maze time')
set(gca, 'XTick', [1,2], 'XTickLabel', {'maze','delay'});

% time absolute
figure(2)
h1 = subplot(1,2,1);
plot([maze_abs.block10;delay_abs.block10],'d-','Color',[0.8,0.8,0.8])
hold on
boxplot(maze_abs.block10,'Position',1);
boxplot(delay_abs.block10,'Position',2);
xlim([0 3])
ylim([0 100])
title('10 sec maze time')
set(gca, 'XTick', [1,2], 'XTickLabel', {'maze','delay'});

h2 = subplot(1,2,2);
plot([maze_abs.block30;delay_abs.block30],'d-','Color',[0.8,0.8,0.8])
hold on
boxplot(maze_abs.block30,'Position',1);
boxplot(delay_abs.block30,'Position',2);
xlim([0 3])
ylim([0 2000])
title('10 sec maze time')
set(gca, 'XTick', [1,2], 'XTickLabel', {'maze','delay'});


% maze length difference
figure(3)
bar(1,(400-75)/400*100)
hold on
bar(2,(75)/400*100)
xlim([0 3])
ylim([0 100])
title('maze length ratio')
set(gca, 'XTick', [1,2], 'XTickLabel', {'maze','delay'});
end