% Plot linear maps for all putative pyr cells on the maze
function Fig8Treadmill_AllLinearMap_Quant(inFile,AnalyzeSes)
close all

p.savePlot = 0;
p.writeToFile = 0;

% need to be adjusted to automatic detection later
p.linROI = [33, 46];
p.linChoice = [55, 64];
p.xTickLin = [16 40 50 60 68];
p.xTickLabel = {'Return','Dl','St','Ch','Rw'};

% Read in input information
sessInfo = SessInfoImport(inFile);
% initiate the map
map_All = [];
map_On = [];
map_Off = [];

for i = AnalyzeSes(1:end)
    
%     if p.savePlot
%         % directory for plot figures
%         % generate a folder for each rat eah day under the current folder
%         savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\Spatial map');
%         if ~exist(savedir, 'dir')
%             mkdir(savedir);
%         end
%     end
    
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);

    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    
    % get valid cell ind
    % time cell / non-time cell
    rateCluster = SpikeProp.AvgRate.Fig8Rate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    clusterNum = length(SpikeProp.AvgRate.Fig8Rate);
    
    
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    spikeMap = cell(clusterNum,1);
    occupMap = cell(clusterNum,1);
    
    spikeMap_On = cell(clusterNum,1);
    spikeMap_Off = cell(clusterNum,1);
    
    occupMap_On = cell(clusterNum,1);
    occupMap_Off = cell(clusterNum,1);   
    
    for j = 1:length(sessDirs)
        
        % take map for each delay condition
        if i == AnalyzeSes(1)
            map_Block.(sessDirs{j}) = [];
        end
        
        % load analyzed map per trial
        load(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'));        

        trialNum = length(ratesByECLR.ECLR);
        
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
                spikeMap_AllTrial = sum(spikeMap{k},1);
                occupMap_AllTrial = sum(occupMap{k},1);
                map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
                map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
                map_Block.(sessDirs{j}) = [map_Block.(sessDirs{j});map_Temp_2];
            end
        end

    end
    
    for k = 1:clusterNum
        if rateLabel(k) == 1
            spikeMap_AllTrial = sum(spikeMap{k},1);
            occupMap_AllTrial = sum(occupMap{k},1);
            map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
            map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
            map_All = [map_All;map_Temp_2];
            
            % treadmill on
            spikeMap_AllTrial = sum(spikeMap_On{k},1);
            occupMap_AllTrial = sum(occupMap_On{k},1);
            map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
            map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
            map_On = [map_On;map_Temp_2];
            
            % treadmill off
            spikeMap_AllTrial = sum(spikeMap_Off{k},1);
            occupMap_AllTrial = sum(occupMap_Off{k},1);
            map_Temp = spikeMap_AllTrial./occupMap_AllTrial;
            map_Temp_2 = gaussfilt(1:length(map_Temp),map_Temp,1);
            map_Off = [map_Off;map_Temp_2];
        end
    end 
    fprintf('Finished position analysis for session %d\n',i);
end

% unnormalized map
figure(1)
[~,peakIdx] = max(map_All,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = map_All(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);

imagesc(cellMapTemp1_Sort)
colormap(jet)
set(gca, 'xtick', p.xTickLin);
set(gca, 'xticklabels', p.xTickLabel);
hold on; vertplot(p.linROI, 0, size(cellMapTemp1_Sort,1), 'r--');
vertplot(p.linChoice, 0, size(cellMapTemp1_Sort,1), 'r--');


figure(2)
imagesc(cellMapTemp1_Sort_Norm)
colormap(jet)
set(gca, 'xtick', p.xTickLin);
set(gca, 'xticklabels', p.xTickLabel);
hold on; vertplot(p.linROI, 0, size(cellMapTemp1_Sort,1), 'r--');
vertplot(p.linChoice, 0, size(cellMapTemp1_Sort,1), 'r--');
% title(TITLEText1,'Interpreter','None')
% axis on
% set(gca, 'xtick', [0 size(cellMapTemp1_Sort_Norm,2)]);
% set(gca, 'xticklabels', TICKText1);

% each delay condition
figure(3)
for j = 1:length(sessDirs)
    subplot(2,4,j)
    cellMapTemp1_Sort = map_Block.(sessDirs{j})(peak_Sort,:);
    cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);
    imagesc(cellMapTemp1_Sort_Norm)
    colormap(jet)
    set(gca, 'xtick', p.xTickLin);
    set(gca, 'xticklabels', p.xTickLabel);
    hold on; vertplot(p.linROI, 0, size(cellMapTemp1_Sort,1), 'r--');
    vertplot(p.linChoice, 0, size(cellMapTemp1_Sort,1), 'r--');
    title(sessDirs{j},'Interpreter','None');   
end

figure(5)
[~,peakIdx] = max(map_On,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = map_All(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);

imagesc(cellMapTemp1_Sort_Norm)
colormap(jet)
set(gca, 'xtick', p.xTickLin);
set(gca, 'xticklabels', p.xTickLabel);
hold on; vertplot(p.linROI, 0, size(cellMapTemp1_Sort,1), 'r--');
vertplot(p.linChoice, 0, size(cellMapTemp1_Sort,1), 'r--');
title('Treadmill on')


figure(6)
[~,peakIdx] = max(map_Off,[],2);
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = map_All(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);

imagesc(cellMapTemp1_Sort_Norm)
colormap(jet)
set(gca, 'xtick', p.xTickLin);
set(gca, 'xticklabels', p.xTickLabel);
hold on; vertplot(p.linROI, 0, size(cellMapTemp1_Sort,1), 'r--');
vertplot(p.linChoice, 0, size(cellMapTemp1_Sort,1), 'r--');
title('Treadmill off')



% get delay active cells
delay_onIdx = (sum(delayMeanRate_Def1>p.avgRateThres,2))>0;

% get maze active cells
mazePeakRate = max(map_All(:,[1:p.linROI(1)-2,p.linROI(2)+2:end]),[],2);
mazeMeanRate = mean(map_All(:,[1:p.linROI(1)-2,p.linROI(2)+2:end]),2);
mazeMeanRate_SD = std(map_All(:,[1:p.linROI(1)-2,p.linROI(2)+2:end]),0,2);

maze_onIdx_1 = mazePeakRate > 3;
maze_onIdx_2 = mazePeakRate > (mazeMeanRate+2*mazeMeanRate_SD);
maze_onIdx = maze_onIdx_1 & maze_onIdx_2;


figure(3)
maze_On_map = map_Pyr(maze_onIdx,:);
[~,peakIdx] = max(maze_On_map,[],2);
% consider cells outside reward area
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = maze_On_map(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);

imagesc(cellMapTemp1_Sort)
colormap(jet)
set(gca, 'xtick', p.xTickLin);
set(gca, 'xticklabels', p.xTickLabel);
hold on; vertplot(p.linROI, 0, size(cellMapTemp1_Sort,1), 'r--');
vertplot(p.linChoice, 0, size(cellMapTemp1_Sort,1), 'r--');
title('All maze on cells');
caxis([0 max(max(cellMapTemp1_Sort))/4])

figure(4)
delay_On_map = map_Pyr(delay_onIdx,:);
[~,peakIdx] = max(delay_On_map,[],2);
% consider cells outside reward area
[~,peak_Sort] = sort(peakIdx);
cellMapTemp1_Sort = delay_On_map(peak_Sort,:);
cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);

imagesc(cellMapTemp1_Sort)
colormap(jet)
set(gca, 'xtick', p.xTickLin);
set(gca, 'xticklabels', p.xTickLabel);
hold on; vertplot(p.linROI, 0, size(cellMapTemp1_Sort,1), 'r--');
vertplot(p.linChoice, 0, size(cellMapTemp1_Sort,1), 'r--');
title('All Delay on cells');
caxis([0 max(max(cellMapTemp1_Sort))/4])

end