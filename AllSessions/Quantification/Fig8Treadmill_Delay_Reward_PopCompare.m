% Plot linear maps for all putative pyr cells on the maze
function Fig8Treadmill_Delay_Reward_PopCompare(inFile,AnalyzeSes)
close all

p.savePlot = 0;
p.writeToFile = 0;

% 
p.cellNumThres = 20;

% Read in input information
sessInfo = SessInfoImport(inFile);

rtFile = 'W:\Li Yuan\Codes\Fig8MazeTreadmill\Cell Property\ReactivationCode-Leibold\matlab code\Results\rt.mat';
load(rtFile);

dayCount = 0;

for i = AnalyzeSes(1:end)
    
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');

    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);

    % load reactivation file
    reactFile = fullfile(sessInfo(i).mainDir,'Cell Property','RewardPopFire_Delay.mat');
    load(reactFile);
    % load reactivation delay specific file
    reactDelayFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_Delay.mat');
    load(reactDelayFile);
    

    % get valid cell ind
    % time cell / non-time cell
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    clusterNum = length(rateCluster);
    pyrNum = sum(rateLabel);
    
    if pyrNum >= p.cellNumThres
        
        if p.savePlot
            % directory for plot figures
            % generate a folder for each rat eah day under the current folder
            savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\Replay Delay Reward Compare');
            if ~exist(savedir, 'dir')
                mkdir(savedir);
            end
        end
        
        % get processed data from each subfolder
        mainDir = sessInfo(i).mainDir;
        % get each phase names (no delay etc)
        sessDirs = sessInfo(i).sessDirs;
        
        % only consider off sessions
        % delay all events
        delay_Num = DelayPopFire_Delay.off10_1.eventNum + DelayPopFire_Delay.off10_2.eventNum + ...
            DelayPopFire_Delay.off30_1.eventNum + DelayPopFire_Delay.off30_2.eventNum;
        delay_eventTsp2Start = [DelayPopFire_Delay.off10_1.eventTsp2Start;DelayPopFire_Delay.off10_2.eventTsp2Start;...
            DelayPopFire_Delay.off30_1.eventTsp2Start;DelayPopFire_Delay.off30_2.eventTsp2Start];
        delay_DirectionLabel = [DelayPopFire_Delay.off10_1.eventTurnDir;DelayPopFire_Delay.off10_2.eventTurnDir;...
            DelayPopFire_Delay.off30_1.eventTurnDir;DelayPopFire_Delay.off30_2.eventTurnDir];
        
%         reward_Num = RewardPopFire_Delay.on10_1.eventNum + RewardPopFire_Delay.on10_2.eventNum + ...
%             RewardPopFire_Delay.on30_1.eventNum + RewardPopFire_Delay.on30_2.eventNum + ...
%             RewardPopFire_Delay.off10_1.eventNum + RewardPopFire_Delay.off10_2.eventNum + ...
%             RewardPopFire_Delay.off30_1.eventNum + RewardPopFire_Delay.off30_2.eventNum;
        
%         reward_eventTsp2Start = [RewardPopFire_Delay.on10_1.eventTsp2Start;RewardPopFire_Delay.on10_2.eventTsp2Start;...
%             RewardPopFire_Delay.on30_1.eventTsp2Start;RewardPopFire_Delay.on30_2.eventTsp2Start; ...
%             RewardPopFire_Delay.off10_1.eventTsp2Start;RewardPopFire_Delay.off10_2.eventTsp2Start;...
%             RewardPopFire_Delay.off30_1.eventTsp2Start;RewardPopFire_Delay.off30_2.eventTsp2Start];
%         
%         reward_DirectionLabel = [RewardPopFire_Delay.on10_1.eventTurnDir;RewardPopFire_Delay.on10_2.eventTurnDir;...
%             RewardPopFire_Delay.on30_1.eventTurnDir;RewardPopFire_Delay.on30_2.eventTurnDir; ...
%             RewardPopFire_Delay.off10_1.eventTurnDir;RewardPopFire_Delay.off10_2.eventTurnDir;...
%             RewardPopFire_Delay.off30_1.eventTurnDir;RewardPopFire_Delay.off30_2.eventTurnDir];

        reward_Num = RewardPopFire_Delay.off10_1.eventNum + RewardPopFire_Delay.off10_2.eventNum + ...
                    RewardPopFire_Delay.off30_1.eventNum + RewardPopFire_Delay.off30_2.eventNum;

        reward_eventTsp2Start = [RewardPopFire_Delay.off10_1.eventTsp2Start;RewardPopFire_Delay.off10_2.eventTsp2Start;...
            RewardPopFire_Delay.off30_1.eventTsp2Start;RewardPopFire_Delay.off30_2.eventTsp2Start];
        
        reward_DirectionLabel = [RewardPopFire_Delay.off10_1.eventTurnDir;RewardPopFire_Delay.off10_2.eventTurnDir;...
            RewardPopFire_Delay.off30_1.eventTurnDir;RewardPopFire_Delay.off30_2.eventTurnDir];

        % left one out for test
        % delay area pop event contigency
        cellActiveEventNum = sum(~isnan(delay_eventTsp2Start),1);
        delay_eventTsp2Start2 = delay_eventTsp2Start(:,cellActiveEventNum>=delay_Num/10);        
        [tsTemp1,seqTemp] = sort(mean(delay_eventTsp2Start2,1,'omitnan'));
        delay_Template{1} = seqTemp(~isnan(tsTemp1));
        
        cellActiveEventNum = sum(~isnan(reward_eventTsp2Start),1);
        reward_eventTsp2Start2 = reward_eventTsp2Start(:,cellActiveEventNum>=reward_Num/10);  
        [tsTemp1,seqTemp] = sort(mean(reward_eventTsp2Start2,1,'omitnan'));
        reward_Template{1} = seqTemp(~isnan(tsTemp1));   
        
        
            
        rval_Delay = zeros(100,1);
        len_Delay = zeros(100,1);
        
        rval_Delay2Reward = zeros(100,1);
        len_Delay2Reward = zeros(100,1);
        
        if delay_Num >=20 && reward_Num>= 20
            dayCount = dayCount + 1;
            
            for n = 1:100
                
                templateIdx_1 = randi([1,delay_Num],floor(delay_Num/2),1);
                templateIdx_2 = 1:delay_Num;
                templateIdx_2 = templateIdx_2(~ismember(templateIdx_2,templateIdx_1));
                
                firstSpkTime = delay_eventTsp2Start2(templateIdx_1,:);
                [tsTemp1,seqTemp] = sort(mean(firstSpkTime,1,'omitnan'));
                seq_template = seqTemp(~isnan(tsTemp1));
                
                firstSpkTime = delay_eventTsp2Start2(templateIdx_2,:);
                [tsTemp2,seqTemp] = sort(mean(firstSpkTime,1,'omitnan'));
                seq_Test{1} = seqTemp(~isnan(tsTemp2));
                
                [rval_Delay(n),len_Delay(n)]=checktempseq(seq_Test,seq_template);
                [rval_Delay2Reward(n),len_Delay2Reward(n)]=checktempseq(reward_Template,seq_template);
                
            end
            
            % reward area pop event contigency
            rval_Reward = zeros(100,1);
            len_Reward = zeros(100,1);
            rval_Reward2Delay = zeros(100,1);
            len_Reward2Delay = zeros(100,1);
            for n = 1:100
                
                templateIdx_1 = randi([1,reward_Num],floor(reward_Num/2),1);
                templateIdx_2 = 1:reward_Num;
                templateIdx_2 = templateIdx_2(~ismember(templateIdx_2,templateIdx_1));
                
                firstSpkTime = reward_eventTsp2Start2(templateIdx_1,:);
                [tsTemp1,seqTemp] = sort(mean(firstSpkTime,1,'omitnan'));
                seq_template = seqTemp(~isnan(tsTemp1));
                
                firstSpkTime = reward_eventTsp2Start2(templateIdx_2,:);
                [tsTemp2,seqTemp] = sort(mean(firstSpkTime,1,'omitnan'));
                seq_Test{1} = seqTemp(~isnan(tsTemp2));
                
                [rval_Reward(n),len_Reward(n)]=checktempseq(seq_Test,seq_template);
                [rval_Reward2Delay(n),len_Reward2Delay(n)]=checktempseq(delay_Template,seq_template);
            end
            
            figure
            Violin(rval_Delay,1,'ShowData',false);
            Violin(rval_Delay2Reward,2,'ShowData',false);
            
            Violin(rval_Reward,4,'ShowData',false);
            Violin(rval_Reward2Delay,5,'ShowData',false);
            
            set(gca,'XTick',[1,2,4,5],'XTickLabel',{'Delay','Delay2Reward','Reward','Reward2Delay'})
            
            rval_Delay_Mean(dayCount) = mean(rval_Delay);
            rval_Delay_SE(dayCount) = std(rval_Delay)./sqrt(length(rval_Delay));
            rval_Delay2Reward_Mean(dayCount) = mean(rval_Delay2Reward);
            rval_Delay2Reward_SE(dayCount) = std(rval_Delay2Reward)./sqrt(length(rval_Delay2Reward));
            
            rval_Reward_Mean(dayCount) = mean(rval_Reward);
            rval_Reward_SE(dayCount) = std(rval_Reward)./sqrt(length(rval_Reward));
            rval_Reward2Delay_Mean(dayCount) = mean(rval_Reward2Delay);
            rval_Reward2Delay_SE(dayCount) = std(rval_Reward2Delay)./sqrt(length(rval_Reward2Delay));            
        end
    end

    fprintf('Finished analysis for session %d\n',i);
end
figure
plot([rval_Delay_Mean;rval_Delay2Reward_Mean],'o-','LineWidth',1.5)
% errorbar([rval_Delay_Mean;rval_Delay2Reward_Mean],[rval_Delay2Reward_Mean;rval_Delay2Reward_SE])
xlim([0,3])
title('Delay compare with delay and reward')
ylim([-0.5,1])
    
figure
plot([rval_Reward_Mean;rval_Reward2Delay_Mean],'o-','LineWidth',1.5)
% errorbar([rval_Delay_Mean;rval_Delay2Reward_Mean],[rval_Delay2Reward_Mean;rval_Delay2Reward_SE])
xlim([0,3])
title('Reward compare with reward and delay')
ylim([-0.5,1])
end