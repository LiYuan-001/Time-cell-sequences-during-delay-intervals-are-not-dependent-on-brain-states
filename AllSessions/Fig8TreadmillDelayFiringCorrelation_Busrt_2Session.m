% correlation of delay area time map is calculated 
% median of all 2 trial pearson correlation was used
% 
% Li Yuan, UCSD, 14-Sep-2022
function Fig8TreadmillDelayFiringCorrelation_Busrt_2Session(inFile,AnalyzeSes,shuffleTimes)

close all

p.savePlot = 1;
p.writeToFile = 1;

% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\DelayFireStability_Burst');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end

    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    % load delayFire map
    delayFile = fullfile(savedir2, 'Fig8DelayTimeMap_Burst.mat');
    load(delayFile);
    
    % initiate the data
    DelayFireStability_Burst_2Session.rat = sessInfo(i).animal;
    DelayFireStability_Burst_2Session.day = sessInfo(i).day;
    DelayFireStability_Burst_2Session.timeBin = DelayFire_Burst.timeBin;
    DelayFireStability_Burst_2Session.gaussSigma = DelayFire_Burst.gaussSigma;
    
    TList = DelayFire_Burst.tList;
    clusterNum = length(DelayFire_Burst.tList);
    DelayFireStability_Burst_2Session.tList = DelayFire_Burst.tList;
            
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    map_Burst_1 = [];
    map_Burst_2 = [];
    map_SingleSpk_1 = [];
    map_SingleSpk_2 = [];
    
    % combine maps from session1 and session2
    for j = 1:length(sessDirs)
        sessName = sessDirs{j}(1:end-2);
        trialNumTemp = size(DelayFire_Burst.(sessDirs{j}).delayTstart1,2);  
        if contains(sessDirs{j},'_1')
            trialNum_2Session.(sessName)(1) = trialNumTemp;
            trialInd = [1:trialNum_2Session.(sessName)(1)];
        elseif contains(sessDirs{j},'_2')
            trialNum_2Session.(sessName)(2) = trialNumTemp;
            trialInd = [trialNum_2Session.(sessName)(1)+1:trialNum_2Session.(sessName)(1)+trialNumTemp];
        end
        binCount = size(DelayFire_Burst.(sessDirs{j}).Burst.spikeRate1_Smooth{1},2);
        % -----------------------------------------------------------------
        for k = 1:clusterNum           
            spikeRate1_Smooth = DelayFire_Burst.(sessDirs{j}).Burst.spikeRate1_Smooth{k};
            map_Burst_1.(sessName)(trialInd,1:binCount,k) = spikeRate1_Smooth;
            spikeRate2_Smooth = DelayFire_Burst.(sessDirs{j}).Burst.spikeRate2_Smooth{k};
            map_Burst_2.(sessName)(trialInd,1:binCount,k) = spikeRate2_Smooth;
            
            spikeRate1_Smooth = DelayFire_Burst.(sessDirs{j}).SingleSpk.spikeRate1_Smooth{k};
            map_SingleSpk_1.(sessName)(trialInd,1:binCount,k) = spikeRate1_Smooth;
            spikeRate2_Smooth = DelayFire_Burst.(sessDirs{j}).SingleSpk.spikeRate2_Smooth{k};
            map_SingleSpk_2.(sessName)(trialInd,1:binCount,k) = spikeRate2_Smooth;
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
        
        %% burst
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
            spikeRate1_2session = map_Burst_1.(sessName2{j})(:,:,k);
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
            spikeRate2_2session = map_Burst_2.(sessName2{j})(:,:,k);
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
               if abs(delayCorr_Def1_Trial(k)) > abs(delayCorr_Def1_Trial_Shuffle95(k))
                   delayCorr_Def1_Trial_Shuffle_Label95(k) = 1;
               end
               if abs(delayCorr_Def1_Trial(k)) > abs(delayCorr_Def1_Trial_Shuffle99(k))
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
               if abs(delayCorr_Def2_Trial(k)) > abs(delayCorr_Def2_Trial_Shuffle95(k))
                   delayCorr_Def2_Trial_Shuffle_Label95(k) = 1;
               end
               if abs(delayCorr_Def2_Trial(k)) > abs(delayCorr_Def2_Trial_Shuffle99(k))
                   delayCorr_Def2_Trial_Shuffle_Label99(k) = 1;
               end
           else
               delayCorr_Def2_Trial_Shuffle95(k) = NaN;
               delayCorr_Def2_Trial_Shuffle99(k) = NaN;
           end               
      
            % plot heat map
            if j == 1                
                h = figure(k);
                if length(sessDirs) == 8
                    h.Position = [100 100 1600 1200];
                else
                    h.Position = [100 100 800 600];
                end
%                 h = figure(k+clusterNum);
%                 if length(sessDirs) == 8
%                     h.Position = [100 100 1600 1200];
%                 else
%                     h.Position = [100 100 800 600];
%                 end
            end
            
            % def 1
            figure(k)
            subplot(4,length(sessDirs)/2,j)   
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
                TITLE1 = sprintf('%s%d%s%d%s%s%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-Burst-',TList{k},'-',sessName2{j});
                TITLE2 = 'RateMap_NormbyTrial_Onset at barrier';
            else
                TITLE1 = sessName2{j};
                TITLE2 = [];
            end
            title({TITLE1;TITLE2},'Interpreter','None')
            xlabel(sprintf('%s%1.2f','Stability: ',delayCorr_Def1_Trial(k)));
            % plot shuffle distribution
            
            subplot(4,length(sessDirs)/2,j+length(sessDirs)/2)
            if ~isnan(delayCorr_Def1_Trial_Shuffle95(k))
                Violin(delayCorr_Def1_Trial_Shuffle(:,k),1);
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
            
                                 

        end  
        
        DelayFireStability_Burst_2Session.(sessName2{j}).Burst.delayCorr_Def1_Trial = delayCorr_Def1_Trial;
%         DelayFireStability_Shuffle.(sessDirs{j}).delayCorr_Trial = delayCorr_Trial;
        DelayFireStability_Burst_2Session.(sessName2{j}).Burst.delayCorr_Def1_Trial_Shuffle = delayCorr_Def1_Trial_Shuffle;
        DelayFireStability_Burst_2Session.(sessName2{j}).Burst.delayCorr_Def1_Trial_Shuffle95 = delayCorr_Def1_Trial_Shuffle95;
        DelayFireStability_Burst_2Session.(sessName2{j}).Burst.delayCorr_Def1_Trial_Shuffle99 = delayCorr_Def1_Trial_Shuffle99;
        DelayFireStability_Burst_2Session.(sessName2{j}).Burst.delayCorr_Def1_Trial_Shuffle_Label95 = delayCorr_Def1_Trial_Shuffle_Label95;      
        DelayFireStability_Burst_2Session.(sessName2{j}).Burst.delayCorr_Def1_Trial_Shuffle_Label99 = delayCorr_Def1_Trial_Shuffle_Label99;    
        
%         DelayFireStability_Shuffle.(sessDirs{j}).delayCorr_Half = delayCorr_Half;        
        DelayFireStability_Burst_2Session.(sessName2{j}).Burst.delayCorr_Def2_Trial = delayCorr_Def2_Trial;
%         DelayFireStability_Shuffle.(sessDirs{j}).delayCorr_Trial = delayCorr_Trial;
        DelayFireStability_Burst_2Session.(sessName2{j}).Burst.delayCorr_Def2_Trial_Shuffle = delayCorr_Def2_Trial_Shuffle;
        DelayFireStability_Burst_2Session.(sessName2{j}).Burst.delayCorr_Def2_Trial_Shuffle95 = delayCorr_Def2_Trial_Shuffle95;
        DelayFireStability_Burst_2Session.(sessName2{j}).Burst.delayCorr_Def2_Trial_Shuffle99 = delayCorr_Def2_Trial_Shuffle99;
        DelayFireStability_Burst_2Session.(sessName2{j}).Burst.delayCorr_Def2_Trial_Shuffle_Label95 = delayCorr_Def2_Trial_Shuffle_Label95;      
        DelayFireStability_Burst_2Session.(sessName2{j}).Burst.delayCorr_Def2_Trial_Shuffle_Label99 = delayCorr_Def2_Trial_Shuffle_Label99;  
        
        %% single spike
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
            spikeRate1_2session = map_SingleSpk_1.(sessName2{j})(:,:,k);
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
            spikeRate2_2session = map_SingleSpk_2.(sessName2{j})(:,:,k);
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
               if abs(delayCorr_Def1_Trial(k)) > abs(delayCorr_Def1_Trial_Shuffle95(k))
                   delayCorr_Def1_Trial_Shuffle_Label95(k) = 1;
               end
               if abs(delayCorr_Def1_Trial(k)) > abs(delayCorr_Def1_Trial_Shuffle99(k))
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
               if abs(delayCorr_Def2_Trial(k)) > abs(delayCorr_Def2_Trial_Shuffle95(k))
                   delayCorr_Def2_Trial_Shuffle_Label95(k) = 1;
               end
               if abs(delayCorr_Def2_Trial(k)) > abs(delayCorr_Def2_Trial_Shuffle99(k))
                   delayCorr_Def2_Trial_Shuffle_Label99(k) = 1;
               end
           else
               delayCorr_Def2_Trial_Shuffle95(k) = NaN;
               delayCorr_Def2_Trial_Shuffle99(k) = NaN;
           end               
      
%             % plot heat map
%             if j == 1                
%                 h = figure(k);
%                 if length(sessDirs) == 8
%                     h.Position = [100 100 1600 1200];
%                 else
%                     h.Position = [100 100 800 600];
%                 end
% %                 h = figure(k+clusterNum);
% %                 if length(sessDirs) == 8
% %                     h.Position = [100 100 1600 1200];
% %                 else
% %                     h.Position = [100 100 800 600];
% %                 end
%             end
            
            % def 1
            figure(k)
            subplot(4,length(sessDirs)/2,j+length(sessDirs)/2*2)   
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
                TITLE1 = sprintf('%s%d%s%d%s%s%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-SingleSpk-',TList{k},'-',sessName2{j});
                TITLE2 = 'RateMap_NormbyTrial_Onset at barrier';
            else
                TITLE1 = sessName2{j};
                TITLE2 = [];
            end
            title({TITLE1;TITLE2},'Interpreter','None')
            xlabel(sprintf('%s%1.2f','Stability: ',delayCorr_Def1_Trial(k)));
            % plot shuffle distribution
            
            subplot(4,length(sessDirs)/2,j+length(sessDirs)/2*3)
            if ~isnan(delayCorr_Def1_Trial_Shuffle95(k))
                Violin(delayCorr_Def1_Trial_Shuffle(:,k),1);
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
            
                                 

        end  
        
        DelayFireStability_Burst_2Session.(sessName2{j}).SingleSpk.delayCorr_Def1_Trial = delayCorr_Def1_Trial;
%         DelayFireStability_Shuffle.(sessDirs{j}).delayCorr_Trial = delayCorr_Trial;
        DelayFireStability_Burst_2Session.(sessName2{j}).SingleSpk.delayCorr_Def1_Trial_Shuffle = delayCorr_Def1_Trial_Shuffle;
        DelayFireStability_Burst_2Session.(sessName2{j}).SingleSpk.delayCorr_Def1_Trial_Shuffle95 = delayCorr_Def1_Trial_Shuffle95;
        DelayFireStability_Burst_2Session.(sessName2{j}).SingleSpk.delayCorr_Def1_Trial_Shuffle99 = delayCorr_Def1_Trial_Shuffle99;
        DelayFireStability_Burst_2Session.(sessName2{j}).SingleSpk.delayCorr_Def1_Trial_Shuffle_Label95 = delayCorr_Def1_Trial_Shuffle_Label95;      
        DelayFireStability_Burst_2Session.(sessName2{j}).SingleSpk.delayCorr_Def1_Trial_Shuffle_Label99 = delayCorr_Def1_Trial_Shuffle_Label99;    
        
%         DelayFireStability_Shuffle.(sessDirs{j}).delayCorr_Half = delayCorr_Half;        
        DelayFireStability_Burst_2Session.(sessName2{j}).SingleSpk.delayCorr_Def2_Trial = delayCorr_Def2_Trial;
%         DelayFireStability_Shuffle.(sessDirs{j}).delayCorr_Trial = delayCorr_Trial;
        DelayFireStability_Burst_2Session.(sessName2{j}).SingleSpk.delayCorr_Def2_Trial_Shuffle = delayCorr_Def2_Trial_Shuffle;
        DelayFireStability_Burst_2Session.(sessName2{j}).SingleSpk.delayCorr_Def2_Trial_Shuffle95 = delayCorr_Def2_Trial_Shuffle95;
        DelayFireStability_Burst_2Session.(sessName2{j}).SingleSpk.delayCorr_Def2_Trial_Shuffle99 = delayCorr_Def2_Trial_Shuffle99;
        DelayFireStability_Burst_2Session.(sessName2{j}).SingleSpk.delayCorr_Def2_Trial_Shuffle_Label95 = delayCorr_Def2_Trial_Shuffle_Label95;      
        DelayFireStability_Burst_2Session.(sessName2{j}).SingleSpk.delayCorr_Def2_Trial_Shuffle_Label99 = delayCorr_Def2_Trial_Shuffle_Label99; 
    end
     
    if p.savePlot == 1
        for k = 1:clusterNum
            figure(k)
            figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-Burst-',TList{k},'-DelayFireStability_Trial_2Session_Onset-barrier');
            print(figName,'-dpng','-r300');
            figure(k+clusterNum)
            figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',TList{k},'-DelayFireStability_Trial_2Session_Onset-entrance');
            print(figName,'-dpng','-r300');
        end
    end
    if p.writeToFile == 1
        save(fullfile(savedir2,'DelayFireStability_Burst_2Session.mat'), 'DelayFireStability_Burst_2Session');
    end
    clear DelayFireStability_Burst_2Session
    close all
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