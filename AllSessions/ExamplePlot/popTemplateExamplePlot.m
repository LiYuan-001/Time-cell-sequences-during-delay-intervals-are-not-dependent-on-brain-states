% example for population activity in delay area for manuscript
% example for popilation activity in delay and reward area
close all

p.savePlot = 0;
p.writeToFile = 0;

% Read in input information
sessInfo = SessInfoImport('W:\Li Yuan\Codes\Fig8MazeTreadmill\Fig8Treadmill_OnOff.xlsx');

for i = 6

%     if p.savePlot
%         % directory for plot figures
%         % generate a folder for each rat eah day under the current folder
%         savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures',sessInfo(i).animal,'-day',sessInfo(i).day,'\Cell pair cofire');
%         if ~exist(savedir, 'dir')
%             mkdir(savedir);
%         end
%     end
    
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);

    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    % get valid cell ind
    % rate [0.1 10] hz on fig 8 maze
    % time cell / non-time cell
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    timeCellLabel = zeros(1,clusterNum);
    
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs; 
    
    %% 
    % delay and reward area compare
    % delay area
    
    % load reactivation delay specific file
    reactDelayFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_Delay.mat');
    load(reactDelayFile);
    
    for j = 4 % session
        eventStartTs = DelayPopFire_Delay.(sessDirs{j}).eventStartTs;
        eventEndTs = DelayPopFire_Delay.(sessDirs{j}).eventEndTs;
        
        meanTs_1 = nanmean(DelayPopFire_Delay.(sessDirs{j}).eventTsp2Start(1:2:end,:));
        [meanTsSort,cellOrder] = sort(meanTs_1);
        eventTspSort = DelayPopFire_Delay.(sessDirs{j}).eventTsp2Start(1:2:end,cellOrder);
        
        for k = 1:length(cellOrder)
            xPoints = [meanTsSort(k)*3;meanTsSort(k)*3];
            yPoints = [k-0.35;k+0.35];
            plot(xPoints,yPoints,'k')
            hold on
        end
        
        meanTs_2 = nanmean(DelayPopFire_Delay.(sessDirs{j}).eventTsp2Start(2:2:end-1,:));
        meanTs_2_sort = meanTs_2(cellOrder);
        for k = 1:length(cellOrder)
            xPoints = [2+meanTs_2_sort(k)*3;2+meanTs_2_sort(k)*3];
            yPoints = [k-0.35;k+0.35];
            plot(xPoints,yPoints,'k')
            hold on
        end
        
%         for m = 7:11
%             for k = 1:length(cellOrder)
%                 xPoints = [m-6+eventTspSort(m,k)*3;m-6+eventTspSort(m,k)*3];
%                 yPoints = [k-0.35;k+0.35];
%                 plot(xPoints,yPoints,'b')
%                 hold on
%             end
%         end
    end
    
    
    
    % load reactivation file
    reactFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_WholeSes.mat');
    load(reactFile);
    % reward area
    for j = 3
        
        validInd = DelayPopFire_WholeSes.(sessDirs{j}).eventLoc == 5;
        eventStartTs = DelayPopFire_WholeSes.(sessDirs{j}).eventStartTs(validInd);
        eventEndTs = DelayPopFire_WholeSes.(sessDirs{j}).eventEndTs(validInd);
        
        validTrial = find(validInd==1);
        
        meanTs_1 = nanmean(DelayPopFire_WholeSes.(sessDirs{j}).eventTsp2Start(validTrial(1:2:end),:));
        [meanTsSort,cellOrder] = sort(meanTs_1);
        eventTspSort = DelayPopFire_WholeSes.(sessDirs{j}).eventTsp2Start(validTrial(1:2:end),cellOrder);
        
        figure
        for k = 1:length(cellOrder)
            xPoints = [meanTsSort(k)*3;meanTsSort(k)*3];
            yPoints = [k-0.35;k+0.35];
            plot(xPoints,yPoints,'k')
            hold on
        end
        
        meanTs_2 = nanmean(DelayPopFire_WholeSes.(sessDirs{j}).eventTsp2Start(validTrial(2:2:end-1),:));
        meanTs_2_sort = meanTs_2(cellOrder);
        for k = 1:length(cellOrder)
            xPoints = [2+meanTs_2_sort(k)*3;2+meanTs_2_sort(k)*3];
            yPoints = [k-0.35;k+0.35];
            plot(xPoints,yPoints,'k')
            hold on
        end
        
%         for m = 5:9
%             for k = 1:length(cellOrder)
%                 xPoints = [m-3+eventTspSort(m,k)*3;m-3+eventTspSort(m,k)*3];
%                 yPoints = [k-0.35;k+0.35];
%                 plot(xPoints,yPoints,'b')
%                 hold on
%             end
%         end
    end
end
