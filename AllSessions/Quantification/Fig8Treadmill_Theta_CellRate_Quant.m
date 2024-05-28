% Plot linear maps for all putative pyr cells on the maze
function Fig8Treadmill_Theta_CellRate_Quant(inFile,AnalyzeSes)
close all

p.savePlot = 0;
p.writeToFile = 0;

% Read in input information
sessInfo = SessInfoImport(inFile);

avgRate = [];
avgRate_On = [];
avgRate_Off = [];

% calculate overall delay firing rate
delayRate_All = [];
delayRate_Theta = [];
delayRate_NonTheta = [];

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
    ThetaCellRateFile = fullfile(sessInfo(i).mainDir,'Cell Property','ThetaCellRate.mat');
    load(ThetaCellRateFile);
    
    
    % get valid cell ind
    % time cell / non-time cell
    rateCluster = SpikeProp.AvgRate.Fig8Rate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    clusterNum = length(rateCluster);
    avgRate = [avgRate,rateCluster];
    
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    spike_All = zeros(1,clusterNum);
    spike_Theta = zeros(1,clusterNum);
    spike_NonTheta = zeros(1,clusterNum);
    
    time_All = 0;
    time_Theta = 0;
    time_NonTheta = 0;
    
    for j = 1:length(sessDirs)
        
        % initiate the rate
        if i == AnalyzeSes(1)            
            rate_Block_All.(sessDirs{j}) = [];
            rate_Block_Theta.(sessDirs{j}) = [];
            rate_Block_NonTheta.(sessDirs{j}) = [];
            
            spike_All_Sess.(sessDirs{j}) = [];
            spike_Theta_Sess.(sessDirs{j}) = [];
            spike_NonTheta_Sess.(sessDirs{j}) = [];
            
            time_All_Sess.(sessDirs{j}) = 0;
            time_Theta_Sess.(sessDirs{j}) = 0;
            time_NonTheta_Sess.(sessDirs{j}) = 0;
        end
        
        rate_Block_All.(sessDirs{j}) = [rate_Block_All.(sessDirs{j}),ThetaCellRate.(sessDirs{j}).rate_All];
        rate_Block_Theta.(sessDirs{j}) = [rate_Block_Theta.(sessDirs{j}),ThetaCellRate.(sessDirs{j}).rate_Theta];
        rate_Block_NonTheta.(sessDirs{j}) = [rate_Block_NonTheta.(sessDirs{j}),ThetaCellRate.(sessDirs{j}).rate_NonTheta];
        
        spike_All_Sess.(sessDirs{j}) = [spike_All_Sess.(sessDirs{j}),sum(ThetaCellRate.(sessDirs{j}).delay_SpNum,1)];
        spike_Theta_Sess.(sessDirs{j}) = [spike_Theta_Sess.(sessDirs{j}),sum(ThetaCellRate.(sessDirs{j}).delay_SpNum_Theta,1)];
        spike_NonTheta_Sess.(sessDirs{j}) = [spike_NonTheta_Sess.(sessDirs{j}),sum(ThetaCellRate.(sessDirs{j}).delay_SpNum_NonTheta,1)];
        
        time_All_Sess.(sessDirs{j}) = time_All_Sess.(sessDirs{j}) + sum(ThetaCellRate.(sessDirs{j}).time_All);
        time_Theta_Sess.(sessDirs{j}) = time_Theta_Sess.(sessDirs{j}) + sum(ThetaCellRate.(sessDirs{j}).thetaTime);
        time_NonTheta_Sess.(sessDirs{j}) = time_NonTheta_Sess.(sessDirs{j}) + sum(ThetaCellRate.(sessDirs{j}).non_thetaTime);
    end
    fprintf('Finished position analysis for session %d\n',i);
    
end

% delayRate_All = ;
% 
% delayRate_Theta = ;
% 
% delayRate_NonTheta = ;
    
pyrIdx = avgRate>0.1 & avgRate<5;
% delayRate_All_pyr = delayRate_All(pyrIdx);
% delayRate_Theta_pyr = delayRate_Theta(pyrIdx);
% delayRate_NonTheta_pyr = delayRate_NonTheta(pyrIdx);

% delay_On cells
avgRateBlock = [rate_Block_All.on10_1;rate_Block_All.on10_2;rate_Block_All.on30_1;rate_Block_All.on30_2;...
    rate_Block_All.off10_1;rate_Block_All.off10_2;rate_Block_All.off30_1;rate_Block_All.off30_2];
avgRateBlock = avgRateBlock(:,pyrIdx);
delay_onIdx = (sum(avgRateBlock>1,1))>0;


% block on vs block off
avgRateBlcok_On = (avgRateBlock(1,:) + avgRateBlock(2,:) + avgRateBlock(3,:)*3 + avgRateBlock(4,:)*3)/(1+1+3+3);    
avgRateBlcok_Off = (avgRateBlock(5,:) + avgRateBlock(6,:) + avgRateBlock(7,:)*3 + avgRateBlock(8,:)*3)/(1+1+3+3);

delayOnCell_BlockOn = avgRateBlcok_On(delay_onIdx);
delayOnCell_BlockOff = avgRateBlcok_Off(delay_onIdx);

avgRateBlcok_On_10 = (avgRateBlock(1,:) + avgRateBlock(2,:))/(1+1);
avgRateBlcok_Off_10 = (avgRateBlock(5,:) + avgRateBlock(6,:))/(1+1);

delayOnCell_BlockOn_10 = avgRateBlcok_On_10(delay_onIdx);
delayOnCell_BlockOff_10 = avgRateBlcok_Off_10(delay_onIdx);

avgRateBlcok_On_30 = (avgRateBlock(3,:) + avgRateBlock(4,:))/(1+1);
avgRateBlcok_Off_30 = (avgRateBlock(7,:) + avgRateBlock(8,:))/(1+1);

delayOnCell_BlockOn_30 = avgRateBlcok_On_30(delay_onIdx);
delayOnCell_BlockOff_30 = avgRateBlcok_Off_30(delay_onIdx);


% % block on theta/non-theta vs block off theta/non-theta
% avgRateBlcok_On_theta = 
% 
% avgRateBlcok_Off_theta = 
% 
% avgRateBlcok_On_nontheta = 
% 
% avgRateBlcok_Off_nontheta = 


% delayOnCell_BlockOn_Theta = avgRateBlcok_On_theta(delay_onIdx);
% delayOnCell_BlockOff_Theta = avgRateBlcok_Off_theta(delay_onIdx);
% delayOnCell_BlockOn_NonTheta = avgRateBlcok_On_nontheta(delay_onIdx);
% delayOnCell_BlockOff_NonTheta = avgRateBlcok_Off_nontheta(delay_onIdx);


% delay on cells
rate = 0:0.1:10;
for kkk = 1:length(rate)
on10_dist(kkk) = sum(avgRateBlcok_On_10>rate(kkk));
off10_dist(kkk) = sum(avgRateBlcok_Off_10>rate(kkk));
on30_dist(kkk) = sum(avgRateBlcok_On_30>rate(kkk));
off30_dist(kkk) = sum(avgRateBlcok_Off_30>rate(kkk));

on_dist(kkk) = sum(avgRateBlcok_On>rate(kkk));
off_dist(kkk) = sum(avgRateBlcok_Off>rate(kkk));

on10_delayOn_dist(kkk) = sum(delayOnCell_BlockOn_10>rate(kkk));
off10_delayOn_dist(kkk) = sum(delayOnCell_BlockOff_10>rate(kkk));
on30_delayOn_dist(kkk) = sum(delayOnCell_BlockOn_30>rate(kkk));
off30_delayOn_dist(kkk) = sum(delayOnCell_BlockOff_30>rate(kkk));

on_delayOn_dist(kkk) = sum(delayOnCell_BlockOn>rate(kkk));
off_delayOn_dist(kkk) = sum(delayOnCell_BlockOff>rate(kkk));
end

figure
plot(rate,on_dist,'r-','LineWidth',1.5)
hold on
plot(rate,off_dist,'k-','LineWidth',1.5)

figure
plot(rate,on_delayOn_dist,'r-','LineWidth',1.5)
hold on
plot(rate,off_delayOn_dist,'k-','LineWidth',1.5)

figure
plot(rate,on10_delayOn_dist,'r-','LineWidth',1.5)
hold on
plot(rate,off10_delayOn_dist,'k-','LineWidth',1.5)

figure
plot(rate,on30_delayOn_dist,'r-','LineWidth',1.5)
hold on
plot(rate,off30_delayOn_dist,'k-','LineWidth',1.5)

end