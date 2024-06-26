% correlation of delay area time map is calculated 
% median of all 2 trial pearson correlation was used
% 
% Li Yuan, UCSD, Jan-03-2022
function Fig8TreadmillDelayFiringCorrelation_2Session(inFile,AnalyzeSes,shuffleTimes)

close all

p.savePlot = 0;
p.writeToFile = 1;

% Read in input information
sessInfo = SessInfoImport_Dual(inFile);

for i = AnalyzeSes(1:end)   
    
    savedir = sprintf('%s%s',sessInfo(i).NIDQ,'\Processed');
    % load delayFire map
    delayFile = fullfile(savedir,'Fig8DelayTimeMap.mat');
    load(delayFile);
    
    
    % initiate the data
    DelayFireStability_2Session.rat = sessInfo(i).ratID;
    DelayFireStability_2Session.day = sessInfo(i).Date;
    DelayFireStability_2Session.timeBin = DelayFire.timeBin;
    DelayFireStability_2Session.gaussSigma = DelayFire.gaussSigma;
    
    %% imec0
    if isfield(DelayFire,'imec0')
        disp(['About to do session ' sessInfo(i).NIDQ,' Imec0']);
        if p.savePlot
            % directory for plot figures
            % generate a folder for each rat eah day under the current folder
            savedir_1 = sprintf('%s%s%d%s%s%s',cd,'\Figures\',sessInfo(i).ratID,'-day',sessInfo(i).Date,'\DelayFireStability\Imec0');
            if ~exist(savedir_1, 'dir')
                mkdir(savedir_1);
            end
            %         delete(strcat(savedir_1,'\*'));
        end
        
        clusterNum = length(DelayFire.imec0.clusterID);
        
        
        DelayFireStability_2Session.imec0.clusterID = DelayFire.imec0.clusterID;
        DelayFireStability_2Session.imec0.clusterDepth = DelayFire.imec0.clusterDepth;
        DelayFireStability_2Session.imec0.clusterCh = DelayFire.imec0.clusterCh;
        
        % get each phase names (no delay etc)
        sessDirs = sessInfo(i).sessDirs;
        
        % combine maps from session1 and session2
        for j = 1:length(sessDirs)
            sessName = sessDirs{j}(1:end-2);
            trialNumTemp = length(DelayFire.imec0.(sessDirs{j}).delayTstart1);
            if contains(sessDirs{j},'_1')
                trialNum_2Session.(sessName)(1) = trialNumTemp;
                trialInd = [1:trialNum_2Session.(sessName)(1)];
            elseif contains(sessDirs{j},'_2')
                trialNum_2Session.(sessName)(2) = trialNumTemp;
                trialInd = [trialNum_2Session.(sessName)(1)+1:trialNum_2Session.(sessName)(1)+trialNumTemp];
            end
            binCount = size(DelayFire.imec0.(sessDirs{j}).spikeRate1_Smooth{1},2);
            % -----------------------------------------------------------------
            for k = 1:clusterNum
                spikeRate1_Smooth = DelayFire.imec0.(sessDirs{j}).spikeRate1_Smooth{k};
                map_1.(sessName)(trialInd,1:binCount,k) = spikeRate1_Smooth;
                spikeRate2_Smooth = DelayFire.imec0.(sessDirs{j}).spikeRate2_Smooth{k};
                map_2.(sessName)(trialInd,1:binCount,k) = spikeRate2_Smooth;
            end
        end
        clear trialNum_2Session
        
        if length(sessDirs) == 8
            sessName2 = {'on10','off10','on30','off30'};
        elseif contains(sessDirs{1},'on')
            sessName2 = {'on10','on30'};
        else
            sessName2 = {'off10','off30'};
        end
        
        for j = 1:length(sessName2)
            maxT = sessName2{j}(end-1:end);
            delayCorr_Def1_Trial = nan(1,clusterNum);
            delayCorr_Def1_Trial_Shuffle = nan(shuffleTimes,clusterNum);
            delayCorr_Def1_Trial_Shuffle95 = nan(1,clusterNum);
            delayCorr_Def1_Trial_Shuffle99 = nan(1,clusterNum);
            delayCorr_Def1_Trial_Shuffle_Label95 = zeros(1,clusterNum);
            delayCorr_Def1_Trial_Shuffle_Label99 = zeros(1,clusterNum);
            
            delayCorr_Def2_Trial = nan(1,clusterNum);
            delayCorr_Def2_Trial_Shuffle = nan(shuffleTimes,clusterNum);
            delayCorr_Def2_Trial_Shuffle95 = nan(1,clusterNum);
            delayCorr_Def2_Trial_Shuffle99 = nan(1,clusterNum);
            delayCorr_Def2_Trial_Shuffle_Label95 = zeros(1,clusterNum);
            delayCorr_Def2_Trial_Shuffle_Label99 = zeros(1,clusterNum);
            
            
            for k = 1:clusterNum
                % def1, delay starts at barrier
                % get map for 2 sessions combined
                spikeRate1_2session = map_1.(sessName2{j})(:,:,k);
                trialNum = size(spikeRate1_2session,1);
                % calculate correlation in a session trial by trial (only for
                % trials with values)
                % in each session, at least 5 trials should be active
                trialSumRate_1 = sum(spikeRate1_2session~=0,2);
                validTrial_1 = find(trialSumRate_1~=0);
                % count valid trials in first session
                halfTrialNumber = floor(trialNum/2);
                validSum1 = sum(trialSumRate_1(1:halfTrialNumber)~=0)>=5;
                validSum2 = sum(trialSumRate_1(halfTrialNumber+1:end)~=0)>=5;
                validState_1 = validSum1+validSum2 >1;
                if ~validState_1
                    delayCorr_Def1_Trial(k) = NaN;
                else
                    corrTemp_1 = nan(nchoosek(length(validTrial_1),2),1);
                    countTemp_1 = 0;
                    for m = 1:length(validTrial_1)-1
                        for mm = m+1:length(validTrial_1)
                            countTemp_1 = countTemp_1+1;
                            corrTemp_1(countTemp_1) = Xcorrelate(spikeRate1_2session(validTrial_1(m),:),spikeRate1_2session(validTrial_1(mm),:));
                        end
                    end
                    delayCorr_Def1_Trial(k) = median(corrTemp_1,'omitnan');
                end
                
                % def2, delay starts at entrance
                % get map for 2 sessions combined
                spikeRate2_2session = map_2.(sessName2{j})(:,:,k);
                % calculate correlation in a session trial by trial (only for
                % trials with values)
                % in each session, at least 5 trials should be active
                trialSumRate_2 = sum(spikeRate2_2session~=0,2);
                validTrial_2 = find(trialSumRate_2~=0);
                % count valid trials in first session
                halfTrialNumber = floor(trialNum/2);
                validSum1 = sum(trialSumRate_2(1:halfTrialNumber)~=0)>=5;
                validSum2 = sum(trialSumRate_2(halfTrialNumber+1:end)~=0)>=5;
                validState_2 = validSum1+validSum2 >1;
                if ~validState_2
                    delayCorr_Def2_Trial(k) = NaN;
                else
                    corrTemp_2 = nan(nchoosek(length(validTrial_2),2),1);
                    countTemp_2 = 0;
                    for m = 1:length(validTrial_2)-1
                        for mm = m+1:length(validTrial_2)
                            countTemp_2 = countTemp_2+1;
                            corrTemp_2(countTemp_2) = Xcorrelate(spikeRate2_2session(validTrial_2(m),:),spikeRate2_2session(validTrial_2(mm),:));
                        end
                    end
                    delayCorr_Def2_Trial(k) = median(corrTemp_2,'omitnan');
                end
                
                %% shuffle each trial rate bins (move in one direction and loop) and then do delayCorr_Half_Rotate
                binCount = size(spikeRate1_2session,2);
                for n = 1:shuffleTimes
                    % shf from 3rd bin to (max-3)th bin
                    shf_1 = ceil(rand(trialNum,1)*(binCount-6))+3;
                    shf_2 = ceil(rand(trialNum,1)*(binCount-6))+3;
                    spikeRate1_2session_Shuffle = nan(trialNum,binCount);
                    spikeRate2_2session_Shuffle = nan(trialNum,binCount);
                    for m = 1:trialNum
                        if shf_1 < binCount
                            shfOrder_1 = [shf_1(m)+1:binCount,1:shf_1(m)];
                        else
                            shfOrder_1 = [shf_1(m):binCount,1:shf_1(m)-1];
                        end
                        spikeRate1_2session_Shuffle(m,:) = spikeRate1_2session(m,shfOrder_1);
                        if shf_2 < binCount
                            shfOrder_2 = [shf_2(m)+1:binCount,1:shf_2(m)];
                        else
                            shfOrder_2 = [shf_2(m):binCount,1:shf_2(m)-1];
                        end
                        spikeRate2_2session_Shuffle(m,:) = spikeRate2_2session(m,shfOrder_2);
                    end
                    
                    % trial by trial shuffled correlation
                    % def 1
                    if ~validState_1
                        delayCorr_Def1_Trial_Shuffle(n,k) = NaN;
                    else
                        corrTemp_1 = nan(nchoosek(length(validTrial_1),2),1);
                        countTemp_1 = 0;
                        for m = 1:length(validTrial_1)-1
                            for mm = m+1:length(validTrial_1)
                                countTemp_1 = countTemp_1+1;
                                corrTemp_1(countTemp_1) = Xcorrelate(spikeRate1_2session_Shuffle(validTrial_1(m),:),spikeRate1_2session_Shuffle(validTrial_1(mm),:));
                            end
                        end
                        delayCorr_Def1_Trial_Shuffle(n,k) = median(corrTemp_1,'omitnan');
                    end
                    % def2
                    if ~validState_2
                        delayCorr_Def2_Trial_Shuffle(n,k) = NaN;
                    else
                        corrTemp_2 = nan(nchoosek(length(validTrial_2),2),1);
                        countTemp_2 = 0;
                        for m = 1:length(validTrial_2)-1
                            for mm = m+1:length(validTrial_2)
                                countTemp_2 = countTemp_2+1;
                                corrTemp_2(countTemp_2) = Xcorrelate(spikeRate2_2session_Shuffle(validTrial_2(m),:),spikeRate2_2session_Shuffle(validTrial_2(mm),:));
                            end
                        end
                        delayCorr_Def2_Trial_Shuffle(n,k) = median(corrTemp_2,'omitnan');
                    end
                end
                
                % get 95% and 99%
                % def 1
                shuffleTemp_1 = delayCorr_Def1_Trial_Shuffle(:,k);
                tempSort_1 = sort(shuffleTemp_1(~isnan(shuffleTemp_1)));
                if length(tempSort_1) > shuffleTimes/10
                    delayCorr_Def1_Trial_Shuffle95(k) = tempSort_1(ceil(length(tempSort_1)*0.95));
                    delayCorr_Def1_Trial_Shuffle99(k) = tempSort_1(ceil(length(tempSort_1)*0.99));
                    if delayCorr_Def1_Trial(k) > delayCorr_Def1_Trial_Shuffle95(k)
                        delayCorr_Def1_Trial_Shuffle_Label95(k) = 1;
                    end
                    if delayCorr_Def1_Trial(k) > delayCorr_Def1_Trial_Shuffle99(k)
                        delayCorr_Def1_Trial_Shuffle_Label99(k) = 1;
                    end
                else
                    delayCorr_Def1_Trial_Shuffle95(k) = NaN;
                    delayCorr_Def1_Trial_Shuffle99(k) = NaN;
                end
                % def 2
                shuffleTemp_2 = delayCorr_Def2_Trial_Shuffle(:,k);
                tempSort_2 = sort(shuffleTemp_2(~isnan(shuffleTemp_2)));
                if length(tempSort_2) > shuffleTimes/10
                    delayCorr_Def2_Trial_Shuffle95(k) = tempSort_2(ceil(length(tempSort_2)*0.95));
                    delayCorr_Def2_Trial_Shuffle99(k) = tempSort_2(ceil(length(tempSort_2)*0.99));
                    if delayCorr_Def2_Trial(k) > delayCorr_Def2_Trial_Shuffle95(k)
                        delayCorr_Def2_Trial_Shuffle_Label95(k) = 1;
                    end
                    if delayCorr_Def2_Trial(k) > delayCorr_Def2_Trial_Shuffle99(k)
                        delayCorr_Def2_Trial_Shuffle_Label99(k) = 1;
                    end
                else
                    delayCorr_Def2_Trial_Shuffle95(k) = NaN;
                    delayCorr_Def2_Trial_Shuffle99(k) = NaN;
                end
                
                if p.savePlot == 1
                    % plot heat map
                    if j == 1
                        h = figure(k);
                        if length(sessDirs) == 8
                            h.Position = [100 100 1600 600];
                        else
                            h.Position = [100 100 800 600];
                        end
                        h = figure(k+clusterNum);
                        if length(sessDirs) == 8
                            h.Position = [100 100 1600 600];
                        else
                            h.Position = [100 100 800 600];
                        end
                    end
                    
                    % def 1
                    figure(k)
                    subplot(2,length(sessDirs)/2,j)
                    imagesc([0,binCount],[1 size(spikeRate1_2session,1)],spikeRate1_2session/max(max(spikeRate1_2session)))
                    %             hold on
                    %             imagesc([0,binCount],size(spikeRate2_2session,1)+1.5,spikeRate2_Combined_Smooth/max(spikeRate2_Combined_Smooth))
                    colormap(jet)
                    xlabel('Time (Sec)')
                    ylabel('Trials')
                    xlim([0 binCount])
                    xTick= [0 binCount];
                    xTickLabel = [0 str2double(maxT)];
                    %             axis off
                    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
                    
                    if j == 1
                        TITLE1 = sprintf('%s%d%s%s%s%d%s%d%s%s','Imec0-DelayZone1-',sessInfo(i).ratID,'-Day-',sessInfo(i).Date,'-Depth-',DelayFire.imec0.clusterDepth(k),'-ID-',DelayFire.imec0.clusterID(k),'-',sessName2{j});
                        TITLE2 = 'RateMap_NormbyTrial_Onset at barrier';
                    else
                        TITLE1 = sessName2{j};
                        TITLE2 = [];
                    end
                    title({TITLE1;TITLE2},'Interpreter','None')
                    xlabel(sprintf('%s%1.2f','Stability: ',delayCorr_Def1_Trial(k)));
                    % plot shuffle distribution
                    
                    subplot(2,length(sessDirs)/2,j+length(sessDirs)/2)
                    if ~isnan(delayCorr_Def1_Trial_Shuffle95(k))
                        Violin(delayCorr_Def1_Trial_Shuffle(:,k),1,'ShowData',false);
                        hold on
                        plot([0.7,1.3],ones(1,2).*delayCorr_Def1_Trial_Shuffle95(k),'r-')
                        plot(1,delayCorr_Def1_Trial(k),'r*')
                        hold off
                    else
                        text(0.1,0.5,'Comparison not exist','FontSize',8);
                    end
                    %             axis off
                    ylabel('shuffle compare')
                    xlabel('Time shuffle')
                    
                    % def 2
                    figure(k+clusterNum)
                    subplot(2,length(sessDirs)/2,j)
                    imagesc([0,binCount],[1 size(spikeRate2_2session,1)],spikeRate2_2session/max(max(spikeRate2_2session)))
                    %             hold on
                    %             imagesc([0,binCount],size(spikeRate2_2session,1)+1.5,spikeRate2_Combined_Smooth/max(spikeRate2_Combined_Smooth))
                    colormap(jet)
                    xlabel('Time (Sec)')
                    ylabel('Trials')
                    xlim([0 binCount])
                    xTick= [0 binCount];
                    xTickLabel = [0 str2double(maxT)];
                    %             axis off
                    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
                    
                    if j == 1
                        TITLE1 = sprintf('%s%d%s%s%s%d%s%d%s%s','Imec0-DelayZone2-',sessInfo(i).ratID,'-Day-',sessInfo(i).Date,'-Depth-',DelayFire.imec0.clusterDepth(k),'-ID-',DelayFire.imec0.clusterID(k),'-',sessName2{j});
                        TITLE2 = 'RateMap_NormbyTrial_Onset at entrance';
                    else
                        TITLE1 = sessName2{j};
                        TITLE2 = [];
                    end
                    title({TITLE1;TITLE2},'Interpreter','None')
                    xlabel(sprintf('%s%1.2f','Stability: ',delayCorr_Def2_Trial(k)));
                    % plot shuffle distribution
                    
                    subplot(2,length(sessDirs)/2,j+length(sessDirs)/2)
                    if ~isnan(delayCorr_Def2_Trial_Shuffle95(k))
                        Violin(delayCorr_Def2_Trial_Shuffle(:,k),1,'ShowData',false);
                        hold on
                        plot([0.7,1.3],ones(1,2).*delayCorr_Def2_Trial_Shuffle95(k),'r-')
                        plot(1,delayCorr_Def2_Trial(k),'r*')
                        hold off
                    else
                        text(0.1,0.5,'Comparison not exist','FontSize',8);
                    end
                    %             axis off
                    ylabel('shuffle compare')
                    xlabel('Time shuffle')
                end
            end
            
            DelayFireStability_2Session.imec0.(sessName2{j}).delayCorr_Def1_Trial = delayCorr_Def1_Trial;
            %         DelayFireStability_Shuffle.(sessDirs{j}).delayCorr_Trial = delayCorr_Trial;
            DelayFireStability_2Session.imec0.(sessName2{j}).delayCorr_Def1_Trial_Shuffle = delayCorr_Def1_Trial_Shuffle;
            DelayFireStability_2Session.imec0.(sessName2{j}).delayCorr_Def1_Trial_Shuffle95 = delayCorr_Def1_Trial_Shuffle95;
            DelayFireStability_2Session.imec0.(sessName2{j}).delayCorr_Def1_Trial_Shuffle99 = delayCorr_Def1_Trial_Shuffle99;
            DelayFireStability_2Session.imec0.(sessName2{j}).delayCorr_Def1_Trial_Shuffle_Label95 = delayCorr_Def1_Trial_Shuffle_Label95;
            DelayFireStability_2Session.imec0.(sessName2{j}).delayCorr_Def1_Trial_Shuffle_Label99 = delayCorr_Def1_Trial_Shuffle_Label99;
            
            %         DelayFireStability_Shuffle.(sessDirs{j}).delayCorr_Half = delayCorr_Half;
            DelayFireStability_2Session.imec0.(sessName2{j}).delayCorr_Def2_Trial = delayCorr_Def2_Trial;
            %         DelayFireStability_Shuffle.(sessDirs{j}).delayCorr_Trial = delayCorr_Trial;
            DelayFireStability_2Session.imec0.(sessName2{j}).delayCorr_Def2_Trial_Shuffle = delayCorr_Def2_Trial_Shuffle;
            DelayFireStability_2Session.imec0.(sessName2{j}).delayCorr_Def2_Trial_Shuffle95 = delayCorr_Def2_Trial_Shuffle95;
            DelayFireStability_2Session.imec0.(sessName2{j}).delayCorr_Def2_Trial_Shuffle99 = delayCorr_Def2_Trial_Shuffle99;
            DelayFireStability_2Session.imec0.(sessName2{j}).delayCorr_Def2_Trial_Shuffle_Label95 = delayCorr_Def2_Trial_Shuffle_Label95;
            DelayFireStability_2Session.imec0.(sessName2{j}).delayCorr_Def2_Trial_Shuffle_Label99 = delayCorr_Def2_Trial_Shuffle_Label99;
        end
        
        if p.savePlot == 1
            
            for k = 1:clusterNum
                figure(k)
                figName = sprintf('%s%s%d%s%s%s%d%s%d%s',savedir_1,'\Imec0-Rat-',sessInfo(i).ratID,'-Day',sessInfo(i).Date,'-Depth-',DelayFire.imec0.clusterDepth(k),'-ID-',DelayFire.imec0.clusterID(k),'-DelayFireStability_2Session_Onset-barrier');
                print(figName,'-dpng','-r300');
                figure(k+clusterNum)
                figName = sprintf('%s%s%d%s%s%s%d%s%d%s',savedir_1,'\Imec0-Rat-',sessInfo(i).ratID,'-Day',sessInfo(i).Date,'-Depth-',DelayFire.imec0.clusterDepth(k),'-ID-',DelayFire.imec0.clusterID(k),'-DelayFireStability_2Session_Onset-entrance');
                print(figName,'-dpng','-r300');
            end
        end
        
        close all
        clear map_1 map_2
    end
    
    %% imec1
    if isfield(DelayFire,'imec1')
        disp(['About to do session ' sessInfo(i).NIDQ,' Imec1']);
        if p.savePlot
            % directory for plot figures
            % generate a folder for each rat eah day under the current folder
            savedir_1 = sprintf('%s%s%d%s%s%s',cd,'\Figures\',sessInfo(i).ratID,'-day',sessInfo(i).Date,'\DelayFireStability\Imec1');
            if ~exist(savedir_1, 'dir')
                mkdir(savedir_1);
            end
            %         delete(strcat(savedir_1,'\*'));
        end
        
        clusterNum = length(DelayFire.imec1.clusterID);
          
        DelayFireStability_2Session.imec1.clusterID = DelayFire.imec1.clusterID;
        DelayFireStability_2Session.imec1.clusterDepth = DelayFire.imec1.clusterDepth;
        DelayFireStability_2Session.imec1.clusterCh = DelayFire.imec1.clusterCh;
        
        % get each phase names (no delay etc)
        sessDirs = sessInfo(i).sessDirs;
        
        % combine maps from session1 and session2
        for j = 1:length(sessDirs)
            sessName = sessDirs{j}(1:end-2);
            trialNumTemp = length(DelayFire.imec1.(sessDirs{j}).delayTstart1);
            if contains(sessDirs{j},'_1')
                trialNum_2Session.(sessName)(1) = trialNumTemp;
                trialInd = [1:trialNum_2Session.(sessName)(1)];
            elseif contains(sessDirs{j},'_2')
                trialNum_2Session.(sessName)(2) = trialNumTemp;
                trialInd = [trialNum_2Session.(sessName)(1)+1:trialNum_2Session.(sessName)(1)+trialNumTemp];
            end
            binCount = size(DelayFire.imec1.(sessDirs{j}).spikeRate1_Smooth{1},2);
            % -----------------------------------------------------------------
            for k = 1:clusterNum
                spikeRate1_Smooth = DelayFire.imec1.(sessDirs{j}).spikeRate1_Smooth{k};
                map_1.(sessName)(trialInd,1:binCount,k) = spikeRate1_Smooth;
                spikeRate2_Smooth = DelayFire.imec1.(sessDirs{j}).spikeRate2_Smooth{k};
                map_2.(sessName)(trialInd,1:binCount,k) = spikeRate2_Smooth;
            end
        end
        clear trialNum_2Session
        
        if length(sessDirs) == 8
            sessName2 = {'on10','off10','on30','off30'};
        elseif contains(sessDirs{1},'on')
            sessName2 = {'on10','on30'};
        else
            sessName2 = {'off10','off30'};
        end
        
        for j = 1:length(sessName2)
            maxT = sessName2{j}(end-1:end);
            delayCorr_Def1_Trial = nan(1,clusterNum);
            delayCorr_Def1_Trial_Shuffle = nan(shuffleTimes,clusterNum);
            delayCorr_Def1_Trial_Shuffle95 = nan(1,clusterNum);
            delayCorr_Def1_Trial_Shuffle99 = nan(1,clusterNum);
            delayCorr_Def1_Trial_Shuffle_Label95 = zeros(1,clusterNum);
            delayCorr_Def1_Trial_Shuffle_Label99 = zeros(1,clusterNum);
            
            delayCorr_Def2_Trial = nan(1,clusterNum);
            delayCorr_Def2_Trial_Shuffle = nan(shuffleTimes,clusterNum);
            delayCorr_Def2_Trial_Shuffle95 = nan(1,clusterNum);
            delayCorr_Def2_Trial_Shuffle99 = nan(1,clusterNum);
            delayCorr_Def2_Trial_Shuffle_Label95 = zeros(1,clusterNum);
            delayCorr_Def2_Trial_Shuffle_Label99 = zeros(1,clusterNum);
            
            
            for k = 1:clusterNum
                % def1, delay starts at barrier
                % get map for 2 sessions combined
                spikeRate1_2session = map_1.(sessName2{j})(:,:,k);
                trialNum = size(spikeRate1_2session,1);
                % calculate correlation in a session trial by trial (only for
                % trials with values)
                % in each session, at least 5 trials should be active
                trialSumRate_1 = sum(spikeRate1_2session~=0,2);
                validTrial_1 = find(trialSumRate_1~=0);
                % count valid trials in first session
                halfTrialNumber = floor(trialNum/2);
                validSum1 = sum(trialSumRate_1(1:halfTrialNumber)~=0)>=5;
                validSum2 = sum(trialSumRate_1(halfTrialNumber+1:end)~=0)>=5;
                validState_1 = validSum1+validSum2 >1;
                if ~validState_1
                    delayCorr_Def1_Trial(k) = NaN;
                else
                    corrTemp_1 = nan(nchoosek(length(validTrial_1),2),1);
                    countTemp_1 = 0;
                    for m = 1:length(validTrial_1)-1
                        for mm = m+1:length(validTrial_1)
                            countTemp_1 = countTemp_1+1;
                            corrTemp_1(countTemp_1) = Xcorrelate(spikeRate1_2session(validTrial_1(m),:),spikeRate1_2session(validTrial_1(mm),:));
                        end
                    end
                    delayCorr_Def1_Trial(k) = median(corrTemp_1,'omitnan');
                end
                
                % def2, delay starts at entrance
                % get map for 2 sessions combined
                spikeRate2_2session = map_2.(sessName2{j})(:,:,k);
                % calculate correlation in a session trial by trial (only for
                % trials with values)
                % in each session, at least 5 trials should be active
                trialSumRate_2 = sum(spikeRate2_2session~=0,2);
                validTrial_2 = find(trialSumRate_2~=0);
                % count valid trials in first session
                halfTrialNumber = floor(trialNum/2);
                validSum1 = sum(trialSumRate_2(1:halfTrialNumber)~=0)>=5;
                validSum2 = sum(trialSumRate_2(halfTrialNumber+1:end)~=0)>=5;
                validState_2 = validSum1+validSum2 >1;
                if ~validState_2
                    delayCorr_Def2_Trial(k) = NaN;
                else
                    corrTemp_2 = nan(nchoosek(length(validTrial_2),2),1);
                    countTemp_2 = 0;
                    for m = 1:length(validTrial_2)-1
                        for mm = m+1:length(validTrial_2)
                            countTemp_2 = countTemp_2+1;
                            corrTemp_2(countTemp_2) = Xcorrelate(spikeRate2_2session(validTrial_2(m),:),spikeRate2_2session(validTrial_2(mm),:));
                        end
                    end
                    delayCorr_Def2_Trial(k) = median(corrTemp_2,'omitnan');
                end
                
                %% shuffle each trial rate bins (move in one direction and loop) and then do delayCorr_Half_Rotate
                binCount = size(spikeRate1_2session,2);
                for n = 1:shuffleTimes
                    % shf from 3rd bin to (max-3)th bin
                    shf_1 = ceil(rand(trialNum,1)*(binCount-6))+3;
                    shf_2 = ceil(rand(trialNum,1)*(binCount-6))+3;
                    spikeRate1_2session_Shuffle = nan(trialNum,binCount);
                    spikeRate2_2session_Shuffle = nan(trialNum,binCount);
                    for m = 1:trialNum
                        if shf_1 < binCount
                            shfOrder_1 = [shf_1(m)+1:binCount,1:shf_1(m)];
                        else
                            shfOrder_1 = [shf_1(m):binCount,1:shf_1(m)-1];
                        end
                        spikeRate1_2session_Shuffle(m,:) = spikeRate1_2session(m,shfOrder_1);
                        if shf_2 < binCount
                            shfOrder_2 = [shf_2(m)+1:binCount,1:shf_2(m)];
                        else
                            shfOrder_2 = [shf_2(m):binCount,1:shf_2(m)-1];
                        end
                        spikeRate2_2session_Shuffle(m,:) = spikeRate2_2session(m,shfOrder_2);
                    end
                    
                    % trial by trial shuffled correlation
                    % def 1
                    if ~validState_1
                        delayCorr_Def1_Trial_Shuffle(n,k) = NaN;
                    else
                        corrTemp_1 = nan(nchoosek(length(validTrial_1),2),1);
                        countTemp_1 = 0;
                        for m = 1:length(validTrial_1)-1
                            for mm = m+1:length(validTrial_1)
                                countTemp_1 = countTemp_1+1;
                                corrTemp_1(countTemp_1) = Xcorrelate(spikeRate1_2session_Shuffle(validTrial_1(m),:),spikeRate1_2session_Shuffle(validTrial_1(mm),:));
                            end
                        end
                        delayCorr_Def1_Trial_Shuffle(n,k) = median(corrTemp_1,'omitnan');
                    end
                    % def2
                    if ~validState_2
                        delayCorr_Def2_Trial_Shuffle(n,k) = NaN;
                    else
                        corrTemp_2 = nan(nchoosek(length(validTrial_2),2),1);
                        countTemp_2 = 0;
                        for m = 1:length(validTrial_2)-1
                            for mm = m+1:length(validTrial_2)
                                countTemp_2 = countTemp_2+1;
                                corrTemp_2(countTemp_2) = Xcorrelate(spikeRate2_2session_Shuffle(validTrial_2(m),:),spikeRate2_2session_Shuffle(validTrial_2(mm),:));
                            end
                        end
                        delayCorr_Def2_Trial_Shuffle(n,k) = median(corrTemp_2,'omitnan');
                    end
                end
                
                % get 95% and 99%
                % def 1
                shuffleTemp_1 = delayCorr_Def1_Trial_Shuffle(:,k);
                tempSort_1 = sort(shuffleTemp_1(~isnan(shuffleTemp_1)));
                if length(tempSort_1) > shuffleTimes/10
                    delayCorr_Def1_Trial_Shuffle95(k) = tempSort_1(ceil(length(tempSort_1)*0.95));
                    delayCorr_Def1_Trial_Shuffle99(k) = tempSort_1(ceil(length(tempSort_1)*0.99));
                    if delayCorr_Def1_Trial(k) > delayCorr_Def1_Trial_Shuffle95(k)
                        delayCorr_Def1_Trial_Shuffle_Label95(k) = 1;
                    end
                    if delayCorr_Def1_Trial(k) > delayCorr_Def1_Trial_Shuffle99(k)
                        delayCorr_Def1_Trial_Shuffle_Label99(k) = 1;
                    end
                else
                    delayCorr_Def1_Trial_Shuffle95(k) = NaN;
                    delayCorr_Def1_Trial_Shuffle99(k) = NaN;
                end
                % def 2
                shuffleTemp_2 = delayCorr_Def2_Trial_Shuffle(:,k);
                tempSort_2 = sort(shuffleTemp_2(~isnan(shuffleTemp_2)));
                if length(tempSort_2) > shuffleTimes/10
                    delayCorr_Def2_Trial_Shuffle95(k) = tempSort_2(ceil(length(tempSort_2)*0.95));
                    delayCorr_Def2_Trial_Shuffle99(k) = tempSort_2(ceil(length(tempSort_2)*0.99));
                    if delayCorr_Def2_Trial(k) > delayCorr_Def2_Trial_Shuffle95(k)
                        delayCorr_Def2_Trial_Shuffle_Label95(k) = 1;
                    end
                    if delayCorr_Def2_Trial(k) > delayCorr_Def2_Trial_Shuffle99(k)
                        delayCorr_Def2_Trial_Shuffle_Label99(k) = 1;
                    end
                else
                    delayCorr_Def2_Trial_Shuffle95(k) = NaN;
                    delayCorr_Def2_Trial_Shuffle99(k) = NaN;
                end
                
                if p.savePlot == 1
                    % plot heat map
                    if j == 1
                        h = figure(k);
                        if length(sessDirs) == 8
                            h.Position = [100 100 1600 600];
                        else
                            h.Position = [100 100 800 600];
                        end
                        h = figure(k+clusterNum);
                        if length(sessDirs) == 8
                            h.Position = [100 100 1600 600];
                        else
                            h.Position = [100 100 800 600];
                        end
                    end
                    
                    % def 1
                    figure(k)
                    subplot(2,length(sessDirs)/2,j)
                    imagesc([0,binCount],[1 size(spikeRate1_2session,1)],spikeRate1_2session/max(max(spikeRate1_2session)))
                    %             hold on
                    %             imagesc([0,binCount],size(spikeRate2_2session,1)+1.5,spikeRate2_Combined_Smooth/max(spikeRate2_Combined_Smooth))
                    colormap(jet)
                    xlabel('Time (Sec)')
                    ylabel('Trials')
                    xlim([0 binCount])
                    xTick= [0 binCount];
                    xTickLabel = [0 str2double(maxT)];
                    %             axis off
                    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
                    
                    if j == 1
                        TITLE1 = sprintf('%s%d%s%s%s%d%s%d%s%s','Imec1-DelayZone1-',sessInfo(i).ratID,'-Day-',sessInfo(i).Date,'-Depth-',DelayFire.imec1.clusterDepth(k),'-ID-',DelayFire.imec1.clusterID(k),'-',sessName2{j});
                        TITLE2 = 'RateMap_NormbyTrial_Onset at barrier';
                    else
                        TITLE1 = sessName2{j};
                        TITLE2 = [];
                    end
                    title({TITLE1;TITLE2},'Interpreter','None')
                    xlabel(sprintf('%s%1.2f','Stability: ',delayCorr_Def1_Trial(k)));
                    % plot shuffle distribution
                    
                    subplot(2,length(sessDirs)/2,j+length(sessDirs)/2)
                    if ~isnan(delayCorr_Def1_Trial_Shuffle95(k))
                        Violin(delayCorr_Def1_Trial_Shuffle(:,k),1,'ShowData',false);
                        hold on
                        plot([0.7,1.3],ones(1,2).*delayCorr_Def1_Trial_Shuffle95(k),'r-')
                        plot(1,delayCorr_Def1_Trial(k),'r*')
                        hold off
                    else
                        text(0.1,0.5,'Comparison not exist','FontSize',8);
                    end
                    %             axis off
                    ylabel('shuffle compare')
                    xlabel('Time shuffle')
                    
                    % def 2
                    figure(k+clusterNum)
                    subplot(2,length(sessDirs)/2,j)
                    imagesc([0,binCount],[1 size(spikeRate2_2session,1)],spikeRate2_2session/max(max(spikeRate2_2session)))
                    %             hold on
                    %             imagesc([0,binCount],size(spikeRate2_2session,1)+1.5,spikeRate2_Combined_Smooth/max(spikeRate2_Combined_Smooth))
                    colormap(jet)
                    xlabel('Time (Sec)')
                    ylabel('Trials')
                    xlim([0 binCount])
                    xTick= [0 binCount];
                    xTickLabel = [0 str2double(maxT)];
                    %             axis off
                    set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
                    
                    if j == 1
                        TITLE1 = sprintf('%s%d%s%s%s%d%s%d%s%s','Imec1-DelayZone2-',sessInfo(i).ratID,'-Day-',sessInfo(i).Date,'-Depth-',DelayFire.imec1.clusterDepth(k),'-ID-',DelayFire.imec1.clusterID(k),'-',sessName2{j});
                        TITLE2 = 'RateMap_NormbyTrial_Onset at entrance';
                    else
                        TITLE1 = sessName2{j};
                        TITLE2 = [];
                    end
                    title({TITLE1;TITLE2},'Interpreter','None')
                    xlabel(sprintf('%s%1.2f','Stability: ',delayCorr_Def2_Trial(k)));
                    % plot shuffle distribution
                    
                    subplot(2,length(sessDirs)/2,j+length(sessDirs)/2)
                    if ~isnan(delayCorr_Def2_Trial_Shuffle95(k))
                        Violin(delayCorr_Def2_Trial_Shuffle(:,k),1,'ShowData',false);
                        hold on
                        plot([0.7,1.3],ones(1,2).*delayCorr_Def2_Trial_Shuffle95(k),'r-')
                        plot(1,delayCorr_Def2_Trial(k),'r*')
                        hold off
                    else
                        text(0.1,0.5,'Comparison not exist','FontSize',8);
                    end
                    %             axis off
                    ylabel('shuffle compare')
                    xlabel('Time shuffle')
                end
            end
            
            DelayFireStability_2Session.imec1.(sessName2{j}).delayCorr_Def1_Trial = delayCorr_Def1_Trial;
            %         DelayFireStability_Shuffle.(sessDirs{j}).delayCorr_Trial = delayCorr_Trial;
            DelayFireStability_2Session.imec1.(sessName2{j}).delayCorr_Def1_Trial_Shuffle = delayCorr_Def1_Trial_Shuffle;
            DelayFireStability_2Session.imec1.(sessName2{j}).delayCorr_Def1_Trial_Shuffle95 = delayCorr_Def1_Trial_Shuffle95;
            DelayFireStability_2Session.imec1.(sessName2{j}).delayCorr_Def1_Trial_Shuffle99 = delayCorr_Def1_Trial_Shuffle99;
            DelayFireStability_2Session.imec1.(sessName2{j}).delayCorr_Def1_Trial_Shuffle_Label95 = delayCorr_Def1_Trial_Shuffle_Label95;
            DelayFireStability_2Session.imec1.(sessName2{j}).delayCorr_Def1_Trial_Shuffle_Label99 = delayCorr_Def1_Trial_Shuffle_Label99;
            
            %         DelayFireStability_Shuffle.(sessDirs{j}).delayCorr_Half = delayCorr_Half;
            DelayFireStability_2Session.imec1.(sessName2{j}).delayCorr_Def2_Trial = delayCorr_Def2_Trial;
            %         DelayFireStability_Shuffle.(sessDirs{j}).delayCorr_Trial = delayCorr_Trial;
            DelayFireStability_2Session.imec1.(sessName2{j}).delayCorr_Def2_Trial_Shuffle = delayCorr_Def2_Trial_Shuffle;
            DelayFireStability_2Session.imec1.(sessName2{j}).delayCorr_Def2_Trial_Shuffle95 = delayCorr_Def2_Trial_Shuffle95;
            DelayFireStability_2Session.imec1.(sessName2{j}).delayCorr_Def2_Trial_Shuffle99 = delayCorr_Def2_Trial_Shuffle99;
            DelayFireStability_2Session.imec1.(sessName2{j}).delayCorr_Def2_Trial_Shuffle_Label95 = delayCorr_Def2_Trial_Shuffle_Label95;
            DelayFireStability_2Session.imec1.(sessName2{j}).delayCorr_Def2_Trial_Shuffle_Label99 = delayCorr_Def2_Trial_Shuffle_Label99;
        end
        
        if p.savePlot == 1
            
            for k = 1:clusterNum
                figure(k)
                figName = sprintf('%s%s%d%s%s%s%d%s%d%s',savedir_1,'\Imec1-Rat-',sessInfo(i).ratID,'-Day',sessInfo(i).Date,'-Depth-',DelayFire.imec1.clusterDepth(k),'-ID-',DelayFire.imec1.clusterID(k),'-DelayFireStability_2Session_Onset-barrier');
                print(figName,'-dpng','-r300');
                figure(k+clusterNum)
                figName = sprintf('%s%s%d%s%s%s%d%s%d%s',savedir_1,'\Imec1-Rat-',sessInfo(i).ratID,'-Day',sessInfo(i).Date,'-Depth-',DelayFire.imec1.clusterDepth(k),'-ID-',DelayFire.imec1.clusterID(k),'-DelayFireStability_2Session_Onset-entrance');
                print(figName,'-dpng','-r300');
            end
        end
        
        close all
        clear map_1 map_2
    end
    
    
    if p.writeToFile == 1
        save(fullfile(savedir,'DelayFireStability_2Session.mat'), 'DelayFireStability_2Session');
    end
    clear DelayFireStability_2Session 
    fprintf('Finished shuffle analysis for session %d\n',i);
end

end

function medValue = corr_Median(corrValue)
if sum(isnan(corrValue))<=2
    medValue = mean(corrValue(~isnan(corrValue)));
else
    % mean is nan
    medValue = mean(corrValue);
end
end