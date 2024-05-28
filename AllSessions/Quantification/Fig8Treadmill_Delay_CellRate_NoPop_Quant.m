% Plot linear maps for all putative pyr cells on the maze
function Fig8Treadmill_Delay_CellRate_NoPop_Quant(inFile,AnalyzeSes)
close all

p.savePlot = 0;
p.writeToFile = 0;
p.cellNumThres = 20;

% Read in input information
sessInfo = SessInfoImport(inFile);

avgRate = [];

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
    % load reactivation delay specific file
    reactDelayFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_Delay.mat');
    load(reactDelayFile);
    
    % get valid cell ind
    % time cell / non-time cell
    rateCluster = SpikeProp.AvgRate.Fig8Rate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    clusterNum = length(rateCluster);

    
    rateCluster2 = SpikeProp.max_AvgRate;
    rateLabel2 = rateCluster2<5 & rateCluster2>0.1;
    pyrNum = sum(rateLabel2);
    
    % initiate the variables
    sessDirs = sessInfo(i).sessDirs;
    for j = 1:length(sessDirs)
        
        % initiate the rate
        if i == AnalyzeSes(1)
        
            rate_Block_All.(sessDirs{j}) = [];
            rate_Block_Theta.(sessDirs{j}) = [];
            rate_Block_NonTheta.(sessDirs{j}) = [];
            
            rate_Block_NoPop.(sessDirs{j}) = [];
        end
    end
            
    if pyrNum >= p.cellNumThres      
        avgRate = [avgRate,rateCluster];
        % get processed data from each subfolder
        mainDir = sessInfo(i).mainDir;
        % get each phase names (no delay etc)
        sessDirs = sessInfo(i).sessDirs;
        
        for j = 1:length(sessDirs)
            
            rate_Block_All.(sessDirs{j}) = [rate_Block_All.(sessDirs{j}),ThetaCellRate.(sessDirs{j}).rate_All];
            rate_Block_Theta.(sessDirs{j}) = [rate_Block_Theta.(sessDirs{j}),ThetaCellRate.(sessDirs{j}).rate_Theta];
            rate_Block_NonTheta.(sessDirs{j}) = [rate_Block_NonTheta.(sessDirs{j}),ThetaCellRate.(sessDirs{j}).rate_NonTheta];
            
            spike_AllTemp = sum(ThetaCellRate.(sessDirs{j}).delay_SpNum,1);
            time_AllTemp = sum(ThetaCellRate.(sessDirs{j}).time_All);
            
            % row is cell, column is event
            popSpikeNumTemp = cellfun('length',DelayPopFire_Delay.(sessDirs{j}).eventTsp);
            % each cell total spike num in the all pop events
            popSpikeNumTemp_2 = sum(popSpikeNumTemp,1);
            popTimeTemp = sum(DelayPopFire_Delay.(sessDirs{j}).eventEndTs-DelayPopFire_Delay.(sessDirs{j}).eventStartTs);
            
            rate_Block_NoPop.(sessDirs{j}) = [rate_Block_NoPop.(sessDirs{j}),(spike_AllTemp - popSpikeNumTemp_2)./(time_AllTemp - popTimeTemp)];
        end
        fprintf('Finished position analysis for session %d\n',i);
    end
end

% delayRate_All = (spike_All_Sess.on10_1 + spike_All_Sess.on10_2 + spike_All_Sess.off10_1 + spike_All_Sess.off10_2 + ...
%     spike_All_Sess.on30_1 + spike_All_Sess.on30_2 + spike_All_Sess.off30_1 +spike_All_Sess.off30_2)/...
%     (time_All_Sess.on10_1 + time_All_Sess.on10_2 + time_All_Sess.off10_1 +time_All_Sess.off10_2 + ...
%     time_All_Sess.on30_1 + time_All_Sess.on30_2 + time_All_Sess.off30_1 +time_All_Sess.off30_2);

% delayRate_Theta = (spike_Theta_Sess.on10_1 + spike_Theta_Sess.on10_2 + spike_Theta_Sess.off10_1 + spike_Theta_Sess.off10_2 + ...
%     spike_Theta_Sess.on30_1 + spike_Theta_Sess.on30_2 + spike_Theta_Sess.off30_1 +spike_Theta_Sess.off30_2)/...
%     (time_Theta_Sess.on10_1 + time_Theta_Sess.on10_2 + time_Theta_Sess.off10_1 +time_Theta_Sess.off10_2 + ...
%     time_Theta_Sess.on30_1 + time_Theta_Sess.on30_2 + time_Theta_Sess.off30_1 +time_Theta_Sess.off30_2);
% 
% delayRate_NonTheta = (spike_All_Sess.on10_1 + spike_All_Sess.on10_2 + spike_All_Sess.off10_1 + spike_All_Sess.off10_2 + ...
%     spike_All_Sess.on30_1 + spike_All_Sess.on30_2 + spike_All_Sess.off30_1 +spike_All_Sess.off30_2)/...
%     (time_All_Sess.on10_1 + time_All_Sess.on10_2 + time_All_Sess.off10_1 +time_All_Sess.off10_2 + ...
%     time_All_Sess.on30_1 + time_All_Sess.on30_2 + time_All_Sess.off30_1 +time_All_Sess.off30_2);
    
pyrIdx = avgRate>0.1 & avgRate<5;
% delayRate_All_pyr = delayRate_All(pyrIdx);
% delayRate_Theta_pyr = delayRate_Theta(pyrIdx);
% delayRate_NonTheta_pyr = delayRate_NonTheta(pyrIdx);

% delay_On cells
avgRateBlock = [rate_Block_All.on10_1;rate_Block_All.on10_2;rate_Block_All.on30_1;rate_Block_All.on30_2;...
    rate_Block_All.off10_1;rate_Block_All.off10_2;rate_Block_All.off30_1;rate_Block_All.off30_2];
delay_onIdx = (sum(avgRateBlock>1,1))>0;
delay_onIdx = delay_onIdx & pyrIdx;

% block on vs block off whole block
avgRateBlcok_On = (rate_Block_All.on10_1 + rate_Block_All.on10_2 + rate_Block_All.on30_1*3 + rate_Block_All.on30_2*3)/(1+1+3+3);    
avgRateBlcok_Off = (rate_Block_All.off10_1 + rate_Block_All.off10_2 + rate_Block_All.off30_1*3 + rate_Block_All.off30_2*3)/(1+1+3+3);

delayOnCell_BlockOn = avgRateBlcok_On(delay_onIdx);
delayOnCell_BlockOff = avgRateBlcok_Off(delay_onIdx);

avgRateBlcok_On_10 = (rate_Block_All.on10_1 + rate_Block_All.on10_2 )/(1+1);

avgRateBlcok_Off_10 = (rate_Block_All.off10_1 + rate_Block_All.off10_2 )/(1+1);

delayOnCell_BlockOn_10 = avgRateBlcok_On_10(delay_onIdx);
delayOnCell_BlockOff_10 = avgRateBlcok_Off_10(delay_onIdx);

avgRateBlcok_On_30 = (rate_Block_All.on30_1 + rate_Block_All.on30_2 )/(1+1);

avgRateBlcok_Off_30 = (rate_Block_All.off30_1 + rate_Block_All.off30_2 )/(1+1);

delayOnCell_BlockOn_30 = avgRateBlcok_On_30(delay_onIdx);
delayOnCell_BlockOff_30 = avgRateBlcok_Off_30(delay_onIdx);

% block on vs block off without pop event
avgRateBlcok_On_noPop = (rate_Block_NoPop.on10_1 + rate_Block_NoPop.on10_2 + rate_Block_NoPop.on30_1*3 + rate_Block_NoPop.on30_2*3)/(1+1+3+3);    

avgRateBlcok_Off_noPop = (rate_Block_NoPop.off10_1 + rate_Block_NoPop.off10_2 + rate_Block_NoPop.off30_1*3 + rate_Block_NoPop.off30_2*3)/(1+1+3+3);

delayOnCell_BlockOn_noPop = avgRateBlcok_On_noPop(delay_onIdx);
delayOnCell_BlockOff_noPop = avgRateBlcok_Off_noPop(delay_onIdx);


avgRateBlcok_On_10_noPop = (rate_Block_NoPop.on10_1 + rate_Block_NoPop.on10_2)/(1+1);  

avgRateBlcok_Off_10_noPop = (rate_Block_NoPop.off10_1 + rate_Block_NoPop.off10_2)/(1+1);  

delayOnCell_BlockOn_10_noPop = avgRateBlcok_On_10_noPop(delay_onIdx);
delayOnCell_BlockOff_10_noPop = avgRateBlcok_Off_10_noPop(delay_onIdx);


avgRateBlcok_On_30_noPop = (rate_Block_NoPop.on30_1 + rate_Block_NoPop.on30_2)/(1+1); 
avgRateBlcok_Off_30_noPop = (rate_Block_NoPop.off30_1 + rate_Block_NoPop.off30_2)/(1+1);  

delayOnCell_BlockOn_30_noPop = avgRateBlcok_On_30_noPop(delay_onIdx);
delayOnCell_BlockOff_30_noPop = avgRateBlcok_Off_30_noPop(delay_onIdx);



% on blocks thetaperiod rate
% off blocks theta period rate
% on blocks non-theta period rate
% off blocks non-theta period rate
figure(1)
Violin(delayOnCell_BlockOn_10,1)
Violin(delayOnCell_BlockOff_10,2)
Violin(delayOnCell_BlockOn_30,4)
Violin(delayOnCell_BlockOff_30,5)
set(gca, 'xtick', [1 2 4 5]);
set(gca, 'xticklabels', {'On10','Off10','On30','Off30'});

[h,pVal] = kstest2(delayOnCell_BlockOn_10,delayOnCell_BlockOff_10)
[h,pVal] = kstest2(delayOnCell_BlockOn_30,delayOnCell_BlockOff_30)
[h,pVal] = kstest2(delayOnCell_BlockOn_10,delayOnCell_BlockOn_30)
[h,pVal] = kstest2(delayOnCell_BlockOff_10,delayOnCell_BlockOff_30)


figure(2)
Violin(delayOnCell_BlockOn_10_noPop,1)
Violin(delayOnCell_BlockOff_10_noPop,2)
Violin(delayOnCell_BlockOn_30_noPop,4)
Violin(delayOnCell_BlockOff_30_noPop,5)
set(gca, 'xtick', [1 2 4 5]);
set(gca, 'xticklabels', {'On10','Off10','On30','Off30'});

[h,pVal] = kstest2(delayOnCell_BlockOn_10_noPop,delayOnCell_BlockOff_10_noPop)
[h,pVal] = kstest2(delayOnCell_BlockOn_30_noPop,delayOnCell_BlockOff_30_noPop)
[h,pVal] = kstest2(delayOnCell_BlockOn_10_noPop,delayOnCell_BlockOn_30_noPop)
[h,pVal] = kstest2(delayOnCell_BlockOff_10_noPop,delayOnCell_BlockOff_30_noPop)

end