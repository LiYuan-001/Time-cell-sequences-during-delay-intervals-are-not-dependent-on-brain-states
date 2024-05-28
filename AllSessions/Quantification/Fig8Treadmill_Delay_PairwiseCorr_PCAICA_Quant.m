function Fig8Treadmill_Delay_PairwiseCorr_PCAICA_Quant(inFile,AnalyzeSes)

close all
p.savePlot = 0;
p.writeToFile = 0;
p.eventThres = 2; % p.eventThres * sd + mean

% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    sessDirs = sessInfo(i).sessDirs;
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load assembly file
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAseembly_DelayReward.mat');
    load(assemblyFile);
        
    clusterNum = length(SpikeProp.max_AvgRate);
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    for j = 1:length(sessDirs)
        
        if i == AnalyzeSes(1)
            cosSim.Delay.(sessDirs{j}) = [];
            cosDiff.Delay.(sessDirs{j}) = [];
            cosSim.Reward.(sessDirs{j}) = [];
            cosDiff.Reward.(sessDirs{j}) = [];
            
            % mean stregnth of all time bins
            strength.Delay.LL.(sessDirs{j}) = [];
            strength.Delay.LR.(sessDirs{j}) = [];
            strength.Delay.RR.(sessDirs{j}) = [];
            strength.Delay.RL.(sessDirs{j}) = [];
            
            strength.Reward.LL.(sessDirs{j}) = [];
            strength.Reward.LR.(sessDirs{j}) = [];
            strength.Reward.RR.(sessDirs{j}) = [];
            strength.Reward.RL.(sessDirs{j}) = [];
            
            % cell number of each assembly
            cellNum.Delay.L.(sessDirs{j}) = [];
            cellNum.Delay.R.(sessDirs{j}) = [];
            cellNum.Reward.L.(sessDirs{j}) = [];
            cellNum.Reward.R.(sessDirs{j}) = [];
            
            % Event number of each assembly
            eventNum.Delay.L.(sessDirs{j}) = [];
            eventNum.Delay.R.(sessDirs{j}) = [];
            eventNum.Reward.L.(sessDirs{j}) = [];
            eventNum.Reward.R.(sessDirs{j}) = [];
            
            eventRate.Delay.L.(sessDirs{j}) = [];
            eventRate.Delay.R.(sessDirs{j}) = [];
            eventRate.Reward.L.(sessDirs{j}) = [];
            eventRate.Reward.R.(sessDirs{j}) = [];
            
            
            % mean stregnth of all time bins
            event_strength.Delay.LL.(sessDirs{j}) = [];
            event_strength.Delay.LR.(sessDirs{j}) = [];
            event_strength.Delay.RR.(sessDirs{j}) = [];
            event_strength.Delay.RL.(sessDirs{j}) = [];
            
            event_strength.Reward.LL.(sessDirs{j}) = [];
            event_strength.Reward.LR.(sessDirs{j}) = [];
            event_strength.Reward.RR.(sessDirs{j}) = [];
            event_strength.Reward.RL.(sessDirs{j}) = [];
        end
        
        if sum(rateLabel) >= 20
            % delay area
            % load analyzed positions
            delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
            load(delayFile);
            % load path zone time
            pathZoneFile = fullfile(mainDir,sessDirs{j}, 'PathZone.mat');
            load(pathZoneFile);
            % load turning directions and correctness
            tInfoFile = fullfile(mainDir,sessDirs{j}, 'trialInfo.mat');
            tInfo = load(tInfoFile);
            
            % write down the results for the delay area
            %         delay.seqTotal_L = seqTotal_L;
            %         delay.seqTotal_R = seqTotal_R;
            %         delay.AssmblWght_L = AssmblWght_L;
            %         delay.AssmblWght_R = AssmblWght_R;
            %         delay.AssmblPtrnCellIDs_L = AssmblPtrnCellIDs_L;
            %         delay.AssmblPtrnCellIDs_R = AssmblPtrnCellIDs_R;
            %         delay.AssmblStrength_LL = AssmblStrength_LL;
            %         delay.AssmblStrength_LR = AssmblStrength_LR;
            %         delay.AssmblStrength_RR = AssmblStrength_RR;
            %         delay.AssmblStrength_RL = AssmblStrength_RL;
            
            
            delay = CellAseembly.(sessDirs{j}).delay;
            z_seqTotal_L = zscore(delay.seqTotal_L')';
            z_seqTotal_R = zscore(delay.seqTotal_R')';
            
            for kk = 1:size(delay.AssmblWght_L,2)
                pattern_L = delay.AssmblWght_L(:,kk);
                cdTemp = [];
                csTemp = [];
                for mm = 1:size(delay.AssmblWght_R,2)

                    pattern_R = delay.AssmblWght_R(:,mm);
                    [cdTemp(mm), csTemp(mm)] = cosmetric(pattern_L, pattern_R);
                    cosSim.Delay.(sessDirs{j}) = [cosSim.Delay.(sessDirs{j}),max(csTemp)];
                    cosDiff.Delay.(sessDirs{j}) = [cosDiff.Delay.(sessDirs{j}),min(cdTemp)];
                end              
            end
            
            strength.Delay.LL.(sessDirs{j}) = [strength.Delay.LL.(sessDirs{j});mean(delay.AssmblStrength_LL,2)];
            strength.Delay.LR.(sessDirs{j}) = [strength.Delay.LR.(sessDirs{j});mean(delay.AssmblStrength_LR,2)];
            strength.Delay.RR.(sessDirs{j}) = [strength.Delay.RR.(sessDirs{j});mean(delay.AssmblStrength_RR,2)];
            strength.Delay.RL.(sessDirs{j}) = [strength.Delay.RL.(sessDirs{j});mean(delay.AssmblStrength_RL,2)];
            
            for kk = 1:size(delay.AssmblWght_L,2)
                % detect cell number, event num and event strength of this
                % pattern
                strengthTemp = delay.AssmblStrength_LL(kk,:);
                cellNum.Delay.L.(sessDirs{j}) = [cellNum.Delay.L.(sessDirs{j}),length(delay.AssmblPtrnCellIDs_L{kk})];
                thresVal = mean(strengthTemp) + p.eventThres * std(strengthTemp);
                strengthLabel = strengthTemp > thresVal;
                strengthLabelDiff = diff([0,strengthLabel]);
                eventNum.Delay.L.(sessDirs{j}) = [eventNum.Delay.L.(sessDirs{j}),sum(strengthLabelDiff==1)];
                timeAll = length(strengthTemp) * CellAseembly.binWidth;
                eventRate.Delay.L.(sessDirs{j}) = [eventRate.Delay.L.(sessDirs{j}),sum(strengthLabelDiff==1)./timeAll];
                event_strength.Delay.LL.(sessDirs{j}) = [event_strength.Delay.LL.(sessDirs{j}),strengthTemp(strengthLabelDiff==1)];
            end
            for kk = 1:size(delay.AssmblWght_R,2)
                % detect cell number, event num and event strength of this
                % pattern
                strengthTemp = delay.AssmblStrength_RR(kk,:);
                cellNum.Delay.R.(sessDirs{j}) = [cellNum.Delay.R.(sessDirs{j}),length(delay.AssmblPtrnCellIDs_R{kk})];
                thresVal = mean(strengthTemp) + p.eventThres * std(strengthTemp);
                strengthLabel = strengthTemp > thresVal;
                strengthLabelDiff = diff([0,strengthLabel]);
                eventNum.Delay.R.(sessDirs{j}) = [eventNum.Delay.R.(sessDirs{j}),sum(strengthLabelDiff==1)];
                timeAll = length(strengthTemp) * CellAseembly.binWidth;
                eventRate.Delay.R.(sessDirs{j}) = [eventRate.Delay.R.(sessDirs{j}),sum(strengthLabelDiff==1)./timeAll];
                event_strength.Delay.RR.(sessDirs{j}) = [event_strength.Delay.RR.(sessDirs{j}),strengthTemp(strengthLabelDiff==1)];
            end
%             figure(1)
%             assemblyNum = size(delay.AssmblWght_L,2);
%             for kk = 1:assemblyNum
%                 subplot(assemblyNum,8,j+(kk-1)*8)
%                 stem(delay.AssmblWght_L(:,kk))
%                 hold on
%                 stem(delay.AssmblPtrnCellIDs_L{kk},delay.AssmblWght_L(delay.AssmblPtrnCellIDs_L{kk},kk),'r')
%                 if kk == 1
%                     title({sessDirs{j};'L assembly delay'})
%                 end
%             end
%             
%             figure(2)
%             h1 = subplot(4,8,j);
%             imagesc(delay.seqTotal_L)
%             title({sessDirs{j};'L assembly delay'})
%             h2 = subplot(4,8,j+8);
%             plot(delay.AssmblStrength_LL')
%             h3 = subplot(4,8,j+16);
%             imagesc(delay.seqTotal_R)
%             h4 = subplot(4,8,j+24);
%             plot(delay.AssmblStrength_LR')
%             linkaxes([h1,h2,h3,h4],'x');
%             linkaxes([h2,h4],'y');
            %% reward area
            reward = CellAseembly.(sessDirs{j}).reward;
            z_seqTotal_L = zscore(reward.seqTotal_L')';
            z_seqTotal_R = zscore(reward.seqTotal_R')';
            
            for kk = 1:size(reward.AssmblWght_L,2)
                cdTemp = [];
                csTemp = [];
                for mm = 1:size(reward.AssmblWght_R,2)
                    pattern_L = reward.AssmblWght_L(:,kk);
                    pattern_R = reward.AssmblWght_R(:,mm);
                    [cdTemp(mm), csTemp(mm)] = cosmetric(pattern_L, pattern_R);
                    cosSim.Reward.(sessDirs{j}) = [cosSim.Reward.(sessDirs{j}),max(csTemp)];
                    cosDiff.Reward.(sessDirs{j}) = [cosDiff.Reward.(sessDirs{j}),min(cdTemp)];
                end
            end
            strength.Reward.LL.(sessDirs{j}) = [strength.Reward.LL.(sessDirs{j});mean(reward.AssmblStrength_LL,2)];
            strength.Reward.LR.(sessDirs{j}) = [strength.Reward.LR.(sessDirs{j});mean(reward.AssmblStrength_LR,2)];
            strength.Reward.RR.(sessDirs{j}) = [strength.Reward.RR.(sessDirs{j});mean(reward.AssmblStrength_RR,2)];
            strength.Reward.RL.(sessDirs{j}) = [strength.Reward.RL.(sessDirs{j});mean(reward.AssmblStrength_RL,2)];
            
            for kk = 1:size(reward.AssmblWght_L,2)
                % detect cell number, event num and event strength of this
                % pattern
                strengthTemp = reward.AssmblStrength_LL(kk,:);
                cellNum.Reward.L.(sessDirs{j}) = [cellNum.Reward.L.(sessDirs{j}),length(reward.AssmblPtrnCellIDs_L{kk})];
                thresVal = mean(strengthTemp) + p.eventThres * std(strengthTemp);
                strengthLabel = strengthTemp > thresVal;
                strengthLabelDiff = diff([0,strengthLabel]);
                eventNum.Reward.L.(sessDirs{j}) = [eventNum.Reward.L.(sessDirs{j}),sum(strengthLabelDiff==1)];
                timeAll = length(strengthTemp) * CellAseembly.binWidth;
                eventRate.Reward.L.(sessDirs{j}) = [eventRate.Reward.L.(sessDirs{j}),sum(strengthLabelDiff==1)./timeAll];
                event_strength.Reward.LL.(sessDirs{j}) = [event_strength.Reward.LL.(sessDirs{j}),strengthTemp(strengthLabelDiff==1)];
            end
            for kk = 1:size(reward.AssmblWght_R,2)
                % detect cell number, event num and event strength of this
                % pattern
                strengthTemp = reward.AssmblStrength_RR(kk,:);
                cellNum.Reward.R.(sessDirs{j}) = [cellNum.Reward.R.(sessDirs{j}),length(reward.AssmblPtrnCellIDs_R{kk})];
                thresVal = mean(strengthTemp) + p.eventThres * std(strengthTemp);
                strengthLabel = strengthTemp > thresVal;
                strengthLabelDiff = diff([0,strengthLabel]);
                eventNum.Reward.R.(sessDirs{j}) = [eventNum.Reward.R.(sessDirs{j}),sum(strengthLabelDiff==1)];
                timeAll = length(strengthTemp) * CellAseembly.binWidth;
                eventRate.Reward.R.(sessDirs{j}) = [eventRate.Reward.R.(sessDirs{j}),sum(strengthLabelDiff==1)./timeAll];
                event_strength.Reward.RR.(sessDirs{j}) = [event_strength.Reward.RR.(sessDirs{j}),strengthTemp(strengthLabelDiff==1)];
            end
            
%             figure(3)
%             assemblyNum = size(reward.AssmblWght_L,2);
%             for kk = 1:assemblyNum
%                 subplot(assemblyNum,8,j+(kk-1)*8)
%                 stem(reward.AssmblWght_L(:,kk))
%                 hold on
%                 stem(reward.AssmblPtrnCellIDs_L{kk},reward.AssmblWght_L(reward.AssmblPtrnCellIDs_L{kk},kk),'r')
%                 if kk == 1
%                     title({sessDirs{j};'L assembly reward'})
%                 end
%             end
%             
%             figure(4)
%             h1 = subplot(4,8,j);
%             imagesc(reward.seqTotal_L)
%             title({sessDirs{j};'L assembly reward'})
%             h2 = subplot(4,8,j+8);
%             plot(reward.AssmblStrength_LL')
%             h3 = subplot(4,8,j+16);
%             imagesc(reward.seqTotal_R)
%             h4 = subplot(4,8,j+24);
%             plot(reward.AssmblStrength_LR')
%             linkaxes([h1,h2,h3,h4],'x');
%             linkaxes([h2,h4],'y');
            
        end
    end
    fprintf('Finished analysis for session %d\n',i)
end

sess = {'on10','off10','on30','off30'};

for j = 1:length(sess)
%     figure(1)
%     Violin(cosSim.Delay.(sessDirs{j}),j);
%     Violin(cosSim.Reward.(sessDirs{j}),j+9);
%     figure(2)
%     Violin(cosDiff.Delay.(sessDirs{j}),j);
%     Violin(cosDiff.Reward.(sessDirs{j}),j+9);
    sess_1 = strcat(sess{j},'_1');
    sess_2 = strcat(sess{j},'_2');
    figure(6)
    temp1 = [strength.Delay.LL.(sess_1);strength.Delay.LL.(sess_2)]';
    temp2 = [strength.Delay.LR.(sess_1);strength.Delay.LR.(sess_2)]';
    sess_1_mean = nanmean(temp1);
    sess_2_mean = nanmean(temp2);
    sess_1_ste = std(temp1) / sqrt(length(temp1));
    sess_2_ste = std(temp2) / sqrt(length(temp2));
    subplot(1,4,j)
    plot([temp1;temp2],'o-','MarkerSize',2,'Color',[0.8,0.8,0.8])
    hold on
    errorbar(1,sess_1_mean,sess_1_ste,'ko','LineWidth',1)
    errorbar(2,sess_2_mean,sess_2_ste,'ko','LineWidth',1)
    title({'Delay L assembly';sess{j}})
    set(gca, 'XTick', [1,2], 'XTickLabel', {'LL','LR'});
    xlim([0 3])
    ylim([0 2])
    
    sess_1 = strcat(sess{j},'_1');
    sess_2 = strcat(sess{j},'_2');
    figure(7)
    temp1 = [strength.Delay.RR.(sess_1);strength.Delay.RR.(sess_2)]';
    temp2 = [strength.Delay.RL.(sess_1);strength.Delay.RL.(sess_2)]';
    sess_1_mean = nanmean(temp1);
    sess_2_mean = nanmean(temp2);
    sess_1_ste = std(temp1) / sqrt(length(temp1));
    sess_2_ste = std(temp2) / sqrt(length(temp2));
    subplot(1,4,j)
    plot([temp1;temp2],'o-','MarkerSize',2,'Color',[0.8,0.8,0.8])
    hold on
    errorbar(1,sess_1_mean,sess_1_ste,'ko','LineWidth',1)
    errorbar(2,sess_2_mean,sess_2_ste,'ko','LineWidth',1)
    title({'Delay R assembly';sess{j}})
    set(gca, 'XTick', [1,2], 'XTickLabel', {'RR','RL'});
    xlim([0 3])
    ylim([0 2])
    
    sess_1 = strcat(sess{j},'_1');
    sess_2 = strcat(sess{j},'_2');
    figure(8)
    temp1 = [strength.Reward.LL.(sess_1);strength.Reward.LL.(sess_2)]';
    temp2 = [strength.Reward.LR.(sess_1);strength.Reward.LR.(sess_2)]';
    sess_1_mean = nanmean(temp1);
    sess_2_mean = nanmean(temp2);
    sess_1_ste = std(temp1) / sqrt(length(temp1));
    sess_2_ste = std(temp2) / sqrt(length(temp2));
    subplot(1,4,j)
    plot([temp1;temp2],'o-','MarkerSize',2,'Color',[0.8,0.8,0.8])
    hold on
    errorbar(1,sess_1_mean,sess_1_ste,'ko','LineWidth',1)
    errorbar(2,sess_2_mean,sess_2_ste,'ko','LineWidth',1)
    title({'Reward L assembly';sess{j}})
    set(gca, 'XTick', [1,2], 'XTickLabel', {'LL','LR'});
    xlim([0 3])
    ylim([0 2])
    
    sess_1 = strcat(sess{j},'_1');
    sess_2 = strcat(sess{j},'_2');
    figure(9)
    temp1 = [strength.Reward.RR.(sess_1);strength.Reward.RR.(sess_2)]';
    temp2 = [strength.Reward.RL.(sess_1);strength.Reward.RL.(sess_2)]';
    sess_1_mean = nanmean(temp1);
    sess_2_mean = nanmean(temp2);
    sess_1_ste = std(temp1) / sqrt(length(temp1));
    sess_2_ste = std(temp2) / sqrt(length(temp2));
    subplot(1,4,j)
    plot([temp1;temp2],'o-','MarkerSize',2,'Color',[0.8,0.8,0.8])
    hold on
    errorbar(1,sess_1_mean,sess_1_ste,'ko','LineWidth',1)
    errorbar(2,sess_2_mean,sess_2_ste,'ko','LineWidth',1)
    title({'Reward R assembly';sess{j}})
    set(gca, 'XTick', [1,2], 'XTickLabel', {'RR','RL'});
    xlim([0 3])
    ylim([0 2])
    
    % Assembly num
    figure(10)
    sess_1 = strcat(sess{j},'_1');
    sess_2 = strcat(sess{j},'_2');
    h1 = subplot(2,1,1);
    temp1 = [length(cellNum.Delay.L.(sess_1))+length(cellNum.Delay.L.(sess_2))];
    plot(j,temp1,'d');
    hold on
    title('Assembly Num in assembly Delay')
    temp1 = [length(cellNum.Delay.R.(sess_1))+length(cellNum.Delay.R.(sess_2))];
    plot(j+5,temp1,'d');
    xlim([0 10])   
    h2 = subplot(2,1,2);
    temp1 = [length(cellNum.Reward.L.(sess_1))+length(cellNum.Reward.L.(sess_2))];
    plot(j,temp1,'d');
    hold on
    title('Assembly Num in assembly Reward')
    temp1 = [length(cellNum.Reward.R.(sess_1))+length(cellNum.Reward.R.(sess_2))];
    plot(j+5,temp1,'d');
    set(gca, 'XTick', [1,6], 'XTickLabel', {'L','R'});
    xlim([0 10])
    ylim('auto')
    linkaxes([h1,h2],'y')
    
    % cell num
    figure(11)
    sess_1 = strcat(sess{j},'_1');
    sess_2 = strcat(sess{j},'_2');
    h1 = subplot(2,1,1);
    temp1 = [cellNum.Delay.L.(sess_1),cellNum.Delay.L.(sess_2)];
    Violin(temp1,j,'ShowData',false);
    title('Cell Num in assembly Delay')
    temp2 = [cellNum.Delay.R.(sess_1),cellNum.Delay.R.(sess_2)];
    Violin(temp1,j+5,'ShowData',false);
    h2 = subplot(2,1,2);
    temp1 = [cellNum.Reward.L.(sess_1),cellNum.Reward.L.(sess_2)];
    title('Cell Num in assembly Reward')
    Violin(temp1,j,'ShowData',false);
    temp2 = [cellNum.Reward.R.(sess_1),cellNum.Reward.R.(sess_2)];
    Violin(temp1,j+5,'ShowData',false);
    set(gca, 'XTick', [1,6], 'XTickLabel', {'L','R'});
    linkaxes([h1,h2],'y')
    
    % event rate is better
    figure(12)
    sess_1 = strcat(sess{j},'_1');
    sess_2 = strcat(sess{j},'_2');
    h1 = subplot(2,1,1);
    temp1 = [eventRate.Delay.L.(sess_1),eventRate.Delay.L.(sess_2)];
    Violin(temp1,j,'ShowData',false);
    title('Event Rate in assembly Delay')
    temp2 = [eventRate.Delay.R.(sess_1),eventRate.Delay.R.(sess_2)];
    Violin(temp1,j+5,'ShowData',false);
    h2 = subplot(2,1,2);
    temp1 = [eventRate.Reward.L.(sess_1),eventRate.Reward.L.(sess_2)];
    title('Event Rate in assembly Reward')
    Violin(temp1,j,'ShowData',false);
    temp2 = [eventRate.Reward.R.(sess_1),eventRate.Reward.R.(sess_2)];
    Violin(temp1,j+5,'ShowData',false);
    set(gca, 'XTick', [1,6], 'XTickLabel', {'L','R'});
    linkaxes([h1,h2],'y')
    
    
    % event strength
    figure(13)
    sess_1 = strcat(sess{j},'_1');
    sess_2 = strcat(sess{j},'_2');
    h1 = subplot(2,1,1);
    temp1 = log2([event_strength.Delay.LL.(sess_1),event_strength.Delay.LL.(sess_2)]);
    Violin(temp1,j,'ShowData',false);
    title('Event Strength in assembly Delay')
    temp2 = log2([event_strength.Delay.RR.(sess_1),event_strength.Delay.RR.(sess_2)]);
    Violin(temp1,j+5,'ShowData',false);
    ylabel('log2(Strength)')
    h2 = subplot(2,1,2);
    temp1 = log2([event_strength.Reward.LL.(sess_1),event_strength.Reward.LL.(sess_2)]);
    title('Event Strength in assembly Reward')
    Violin(temp1,j,'ShowData',false);
    temp2 = log2([event_strength.Reward.RR.(sess_1),event_strength.Reward.RR.(sess_2)]);
    Violin(temp1,j+5,'ShowData',false);
    set(gca, 'XTick', [1,6], 'XTickLabel', {'L','R'});
    ylabel('log2(Strength)')
    linkaxes([h1,h2],'y')
end

end