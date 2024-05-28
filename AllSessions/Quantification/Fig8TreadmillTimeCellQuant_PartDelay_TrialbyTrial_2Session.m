function Fig8TreadmillTimeCellQuant_PartDelay_TrialbyTrial_2Session(inFile,AnalyzeSes)

close all

p.savePlot = 0;
p.writeToFile = 0;

avgRate = [];
fieldLabel_Def2.allSession = [];
sameTypeStability.on10Corr = [];
sameTypeStability.off10Corr = [];
sameTypeStability.on30Corr = [];
sameTypeStability.off30Corr = [];

SecName = {'first','second','third'};
% % Read in input information
sessInfo = SessInfoImport(inFile);

if p.savePlot
    % directory for plot figures
    % generate a folder for each rat eah day under the current folder
    savedir = sprintf('%s%s',cd,'\Figures\Quantification Figures\Time Cell Num and stability quantification');
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end
end
    
for i = AnalyzeSes(1:end)
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load stability file
    stabilityFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayFireStability30to10_Shuffle.mat');
    load(stabilityFile);
    % load time field file
    timeFieldFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayTimeField_PartDelay.mat');
    load(timeFieldFile);
    % load delayFire map
    delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap.mat');
    load(delayFile);
    
    idx = find(SpikeProp.AvgRate.Fig8Rate>0.1 & SpikeProp.AvgRate.Fig8Rate<10);
    
    % get each phase names (no delay etc)
    SessDirs = sessInfo(i).sessDirs;
    avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate(idx)];
    for section = 1:3
        fieldLabelTemp.(SecName{section}) = zeros(1,length(idx));
        if i == AnalyzeSes(1)
            fieldLabel_Def2.allSession.(SecName{section}) = [];
        end
    end
    for j = 1:length(SessDirs)
        if contains(SessDirs{j},'30')
            % get time map length for this session
            spikeRate2_Smooth = DelayFire.(SessDirs{j}).spikeRate2_Smooth{1};
            binNum = floor(size(spikeRate2_Smooth,2)/3);
            for section = 1:3
                if i == AnalyzeSes(1)
                    fieldLabel_Def1.(SessDirs{j}).(SecName{section}) = [];
                    fieldLabel_Def2.(SessDirs{j}).(SecName{section}) = [];          
                    timeMap_Def1.(SessDirs{j}).(SecName{section}) = [];
                    timeMap_Def2.(SessDirs{j}).(SecName{section}) = [];
                end
                fieldLabel_Def1.(SessDirs{j}).(SecName{section}) = [fieldLabel_Def1.(SessDirs{j}).(SecName{section}),DelayField_PartDelay.(SessDirs{j}).(SecName{section}).timeField_1(idx)];
                fieldLabel_Def2.(SessDirs{j}).(SecName{section}) = [fieldLabel_Def2.(SessDirs{j}).(SecName{section}),DelayField_PartDelay.(SessDirs{j}).(SecName{section}).timeField_2(idx)];
                % get map for each section
                % map = [cell number, 10 sec map]
                for k = idx
                    
                    spikeRate1_Combined_Smooth = DelayFire.(SessDirs{j}).spikeRate1_Combined_Smooth{k};
                    spikeRate1_Combined_Smooth_Temp = spikeRate1_Combined_Smooth(:,(section-1)*binNum+1:section*binNum);
                    timeMap_Def1.(SessDirs{j}).(SecName{section}) = [timeMap_Def1.(SessDirs{j}).(SecName{section});spikeRate1_Combined_Smooth_Temp];
                    
                    spikeRate2_Combined_Smooth = DelayFire.(SessDirs{j}).spikeRate2_Combined_Smooth{k};
                    spikeRate2_Combined_Smooth_Temp = spikeRate2_Combined_Smooth(:,(section-1)*binNum+1:section*binNum);                 
                    timeMap_Def2.(SessDirs{j}).(SecName{section}) = [timeMap_Def2.(SessDirs{j}).(SecName{section});spikeRate2_Combined_Smooth_Temp];
                end               
            end
        end
    end
end

close all
%% calculate % of time cells
% stability among sessions
% timecellQuant_PartDelay.Num.All = length(fieldLabel.on10_1);
% timecellQuant_PartDelay.Num.AllTimeCells = sum(fieldLabel.allSession);
% % on10_1 and on10_2
% timecellQuant_PartDelay.Num.on10_1 = sum(fieldLabel.on10_1);
% timecellQuant_PartDelay.Num.on10_2 = sum(fieldLabel.on10_2);
% timecellQuant_PartDelay.Num.on10 = mean(timecellQuant_PartDelay.Num.on10_1,timecellQuant_PartDelay.Num.on10_2);
% timecellQuant_PartDelay.Num.on10_Common = sum((fieldLabel.on10_1+fieldLabel.on10_2)==2);
% 

% on30_1 and on30_2 in its own session
h1 = figure(1);
h1.Position = [100 100 1200 600];

% on 30 session 0-10 sec
cellMapTemp1 = timeMap_Def1.on30_1.first(fieldLabel_Def1.on30_1.first==1,:);
cellMapTemp2 = timeMap_Def1.on30_2.first(fieldLabel_Def1.on30_2.first==1,:);
% sort by peak from on30_1 first block
plotTimeMap2(cellMapTemp1,cellMapTemp2,1,[2,6,1],[2,6,7],...
    {'Delay start at barrier','Time cells on30_1 0-10 sec'},{'0','10'},{'Time cells on30_2 0-10 sec'},{'0','10'});
ymax_1 = size(cellMapTemp1,2);
ymax_2 = size(cellMapTemp2,2);

% on 30 session 10-20 sec
cellMapTemp1 = timeMap_Def1.on30_1.second(fieldLabel_Def1.on30_1.second==1,:);
cellMapTemp2 = timeMap_Def1.on30_2.second(fieldLabel_Def1.on30_2.second==1,:);

plotTimeMap2(cellMapTemp1,cellMapTemp2,1,[2,6,2],[2,6,8],...
    {'Delay start at barrier','Time cells on30_1 10-20 sec'},{'0','10'},{'Time cells on30_2 10-20 sec'},{'0','10'});
subplot(2,6,2)
ylim([0 ymax_1])
subplot(2,6,8)
ylim([0 ymax_2])


% on 30 session 20-30 sec
cellMapTemp1 = timeMap_Def1.on30_1.third(fieldLabel_Def1.on30_1.third==1,:);
cellMapTemp2 = timeMap_Def1.on30_2.third(fieldLabel_Def1.on30_2.third==1,:);
% sort by peak from on30_1 third block
plotTimeMap2(cellMapTemp1,cellMapTemp2,1,[2,6,3],[2,6,9],...
    {'Delay start at barrier','Time cells on30_1 20-30 sec'},{'0','10'},{'Time cells on30_2 20-30 sec'},{'0','10'});
subplot(2,6,3)
ylim([0 ymax_1])
subplot(2,6,9)
ylim([0 ymax_2])


% off30 session 0-10 sec
cellMapTemp1 = timeMap_Def1.off30_1.first(fieldLabel_Def1.off30_1.first==1,:);
cellMapTemp2 = timeMap_Def1.off30_2.first(fieldLabel_Def1.off30_2.first==1,:);
plotTimeMap2(cellMapTemp1,cellMapTemp2,1,[2,6,4],[2,6,10],...
    {'Delay start at barrier','Time cells off30_1 0-10 sec'},{'0','10'},{'Time cells off30_2 0-10 sec'},{'0','10'});
ymax_1 = size(cellMapTemp1,2);
ymax_2 = size(cellMapTemp2,2);

% off30 session 10-20 sec
cellMapTemp1 = timeMap_Def1.off30_1.second(fieldLabel_Def1.off30_1.second==1,:);
cellMapTemp2 = timeMap_Def1.off30_2.second(fieldLabel_Def1.off30_2.second==1,:);

plotTimeMap2(cellMapTemp1,cellMapTemp2,1,[2,6,5],[2,6,11],...
    {'Delay start at barrier','Time cells off30_1 10-20 sec'},{'0','10'},{'Time cells off30_2 10-20 sec'},{'0','10'});
subplot(2,6,5)
ylim([0 ymax_1])
subplot(2,6,11)
ylim([0 ymax_2])


% off30 session 20-30 sec
cellMapTemp1 = timeMap_Def1.off30_1.third(fieldLabel_Def1.off30_1.third==1,:);
cellMapTemp2 = timeMap_Def1.off30_2.third(fieldLabel_Def1.off30_2.third==1,:);

plotTimeMap2(cellMapTemp1,cellMapTemp2,1,[2,6,6],[2,6,12],...
    {'Delay start at barrier','Time cells off30_1 20-30 sec'},{'0','10'},{'Time cells off30_2 20-30 sec'},{'0','10'});
subplot(2,6,6)
ylim([0 ymax_1])
subplot(2,6,12)
ylim([0 ymax_2])



% on30_1 and on30_2 sorted by on30_1 time cells etc
h2 = figure(2);
h2.Position = [100 100 1200 600];

% on 30 session 0-10 sec
cellMapTemp1 = timeMap_Def1.on30_1.first(fieldLabel_Def1.on30_1.first==1,:);
cellMapTemp2 = timeMap_Def1.on30_2.first(fieldLabel_Def1.on30_1.first==1,:);
% sort by peak from on30_1 first block
plotTimeMap2(cellMapTemp1,cellMapTemp2,2,[2,6,1],[2,6,7],...
    {'Delay start at barrier','Time cells on30_1 0-10 sec'},{'0','10'},{'Time cells on30_2 0-10 sec sorted by on30_1'},{'0','10'});
ymax_1 = size(cellMapTemp1,2);


% on 30 session 10-20 sec
cellMapTemp1 = timeMap_Def1.on30_1.second(fieldLabel_Def1.on30_1.second==1,:);
cellMapTemp2 = timeMap_Def1.on30_2.second(fieldLabel_Def1.on30_1.second==1,:);

plotTimeMap2(cellMapTemp1,cellMapTemp2,2,[2,6,2],[2,6,8],...
    {'Delay start at barrier','Time cells on30_1 10-20 sec'},{'0','10'},{'Time cells on30_2 10-20 sec sorted by on30_1'},{'0','10'});
subplot(2,6,2)
ylim([0 ymax_1])
subplot(2,6,8)
ylim([0 ymax_1])


% on 30 session 20-30 sec
cellMapTemp1 = timeMap_Def1.on30_1.third(fieldLabel_Def1.on30_1.third==1,:);
cellMapTemp2 = timeMap_Def1.on30_2.third(fieldLabel_Def1.on30_1.third==1,:);
% sort by peak from on30_1 third block
plotTimeMap2(cellMapTemp1,cellMapTemp2,2,[2,6,3],[2,6,9],...
    {'Delay start at barrier','Time cells on30_1 20-30 sec'},{'0','10'},{'Time cells on30_2 20-30 sec sorted by on30_1'},{'0','10'});
subplot(2,6,3)
ylim([0 ymax_1])
subplot(2,6,9)
ylim([0 ymax_1])


% off30 session 0-10 sec
cellMapTemp1 = timeMap_Def1.off30_1.first(fieldLabel_Def1.off30_1.first==1,:);
cellMapTemp2 = timeMap_Def1.off30_2.first(fieldLabel_Def1.off30_1.first==1,:);
plotTimeMap2(cellMapTemp1,cellMapTemp2,2,[2,6,4],[2,6,10],...
    {'Delay start at barrier','Time cells off30_1 0-10 sec'},{'0','10'},{'Time cells off30_2 0-10 sec sorted by off30_1'},{'0','10'});
ymax_1 = size(cellMapTemp1,2);


% off30 session 10-20 sec
cellMapTemp1 = timeMap_Def1.off30_1.second(fieldLabel_Def1.off30_1.second==1,:);
cellMapTemp2 = timeMap_Def1.off30_2.second(fieldLabel_Def1.off30_1.second==1,:);

plotTimeMap2(cellMapTemp1,cellMapTemp2,2,[2,6,5],[2,6,11],...
    {'Delay start at barrier','Time cells off30_1 10-20 sec'},{'0','10'},{'Time cells off30_2 10-20 sec sorted by off30_1'},{'0','10'});
subplot(2,6,5)
ylim([0 ymax_1])
subplot(2,6,11)
ylim([0 ymax_1])


% off30 session 20-30 sec
cellMapTemp1 = timeMap_Def1.off30_1.third(fieldLabel_Def1.off30_1.third==1,:);
cellMapTemp2 = timeMap_Def1.off30_2.third(fieldLabel_Def1.off30_1.third==1,:);

plotTimeMap2(cellMapTemp1,cellMapTemp2,2,[2,6,6],[2,6,12],...
    {'Delay start at barrier','Time cells off30_1 20-30 sec'},{'0','10'},{'Time cells off30_2 20-30 sec sorted by off30_1'},{'0','10'});
subplot(2,6,6)
ylim([0 ymax_1])
subplot(2,6,12)
ylim([0 ymax_1])


if p.writeToFile == 1
    save(fullfile(savedir,'timecellQuant_PartDelay.mat'), 'timecellQuant_PartDelay');
end
clear timecellQuant_PartDelay

end