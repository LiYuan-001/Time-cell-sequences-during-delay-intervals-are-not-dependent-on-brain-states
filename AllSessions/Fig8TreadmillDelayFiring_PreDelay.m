function Fig8TreadmillDelayFiring_PreDelay(inFile,AnalyzeSes)

p.savePlot = 0;
p.writeToFile = 1;

% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'

p.timeBin = 150; % ms;
p.gaussSigma = p.timeBin*2; % sd of the gaussian distribution
p.preDelayTime = 3000; % unit ms
p.postDelayTime = 10*10^3; % unit ms

% Read in input information
sessInfo = SessInfoImport(inFile);
close all

for i = AnalyzeSes(1:end)
    
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\Delay time map');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end
    
    if p.writeToFile
        savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
        if ~exist(savedir2, 'dir')
            mkdir(savedir2);
        end
    end
    
    % initiate the data
    DelayFire.rat = sessInfo(i).animal;
    DelayFire.day = sessInfo(i).day;
    DelayFire.timeBin = p.timeBin;
    DelayFire.gaussSigma = p.gaussSigma;
    DelayFire.preDelayTime = p.preDelayTime;
    
    
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);
    if length(unique(TList))~=clusterNum
        error('TTList file has repeated clusters')
    end

    for k = 1:clusterNum
        DelayFire.tList{k} = TList{k}(1:end-2);
    end
            
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    
    % get each phase names (no delay etc)
    SessDirs = sessInfo(i).sessDirs;
    
    for j = 1:length(SessDirs)
        % load analyzed positions
        delayFile = fullfile(mainDir,SessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        pathFile = fullfile(mainDir,SessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        pathZoneFile = fullfile(sessInfo(i).mainDir,SessDirs{j}, 'PathZone.mat');
        load(pathZoneFile);
 
        % -----------------------------------------------------------------
        % use the mode of delay is always 10 sec or 30 sec
        % rather than identify real time in delay zone
        % to make plot cleaner
        % I can change to delay time in delay zone later
        % Li Yuan, 19-Aug-2020, UCSD
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT*10^3; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT*10^3;
        delayTstart2 = Fig8DelayZonePos.delayPos2.startT*10^3;
        delayTend2 = Fig8DelayZonePos.delayPos2.endT*10^3;
            
        trialNum = size(delayTstart1,2);
        if contains(SessDirs{j},'10')
            maxT = 10*10^3;
            delayTend1_2 = delayTstart1+maxT;
            delayTend2_2 = delayTstart2+maxT;
            delayBinCount = floor(maxT/p.timeBin);
        elseif contains(SessDirs{j},'30')
            maxT = 30*10^3;
            delayTend1_2 = delayTstart1+maxT;
            delayTend2_2 = delayTstart2+maxT;
            delayBinCount = floor(maxT/p.timeBin);
        else
            error('Delay time is wrong')
        end       
        preBinCount = floor(p.preDelayTime/p.timeBin);
        postBinCount = floor(p.postDelayTime/p.timeBin);
        
        DelayFire.(SessDirs{j}).delayTstart1 = delayTstart1;
        DelayFire.(SessDirs{j}).delayTstart2 = delayTstart2;
        DelayFire.(SessDirs{j}).delayTend1 = delayTend1_2;
        DelayFire.(SessDirs{j}).delayTend2 = delayTend2_2;
         
        % this part was for delay rate bins of positions in delay zone
        % since every delay is slightly different length, thus require a
        % cell variable
        % -----------------------------------------------------------------
%         rateBin1 = cell(trialNum,1);
%         rateBin2 = cell(trialNum,1);
%         for m = 1:trialNum
%             rateBin1{m} = delayTstart1(m):p.timeBin:delayTend1(m);
%             rateBin2{m} = delayTstart2(m):p.timeBin:delayTend2(m);
%         end
%         [~,maxIdx1] = max(delayTend1-delayTstart1);
%         delayBinCount = length(rateBin1{maxIdx1});
%         
%         [~,maxIdx2] = max(delayTend2-delayTstart2);
%         delayBinCount = length(rateBin2{maxIdx2});
        % -----------------------------------------------------------------
        
        for k = 1:clusterNum
            % delay area
            spikeRate1 = zeros(trialNum,delayBinCount);
            spikeRate1_Smooth = zeros(trialNum,delayBinCount);
            spikeRateNorm1 = zeros(trialNum,delayBinCount);
            spikeRateNorm1_Smooth = zeros(trialNum,delayBinCount);
            
            spikeRate2 = zeros(trialNum,delayBinCount);
            spikeRate2_Smooth = zeros(trialNum,delayBinCount);
            spikeRateNorm2 = zeros(trialNum,delayBinCount);
            spikeRateNorm2_Smooth = zeros(trialNum,delayBinCount);
            
            % pre area
            pre_spikeRate1 = zeros(trialNum,preBinCount);
            pre_spikeRate1_Smooth = zeros(trialNum,preBinCount);
            pre_spikeRateNorm1 = zeros(trialNum,preBinCount);
            pre_spikeRateNorm1_Smooth = zeros(trialNum,preBinCount);
            
            pre_spikeRate2 = zeros(trialNum,preBinCount);
            pre_spikeRate2_Smooth = zeros(trialNum,preBinCount);
            pre_spikeRateNorm2 = zeros(trialNum,preBinCount);
            pre_spikeRateNorm2_Smooth = zeros(trialNum,preBinCount);
            
            % post delay
            post_spikeRate1 = zeros(trialNum,postBinCount);
            post_spikeRate1_Smooth = zeros(trialNum,postBinCount);
            post_spikeRateNorm1 = zeros(trialNum,postBinCount);
            post_spikeRateNorm1_Smooth = zeros(trialNum,postBinCount);
            
            
            % pre & delay area 
            % smooth and normalize together
            pre_delay_spikeRate1 = zeros(trialNum,preBinCount+delayBinCount);
            pre_delay_spikeRate1_Smooth = zeros(trialNum,preBinCount+delayBinCount);
            pre_delay_spikeRateNorm1 = zeros(trialNum,preBinCount+delayBinCount);
            pre_delay_spikeRateNorm1_Smooth = zeros(trialNum,preBinCount+delayBinCount);
            
            pre_delay_spikeRate2 = zeros(trialNum,preBinCount+delayBinCount);
            pre_delay_spikeRate2_Smooth = zeros(trialNum,preBinCount+delayBinCount);
            pre_delay_spikeRateNorm2 = zeros(trialNum,preBinCount+delayBinCount);
            pre_delay_spikeRateNorm2_Smooth = zeros(trialNum,preBinCount+delayBinCount);  
            
            pre_delay_post_spikeRate1 = zeros(trialNum,preBinCount+delayBinCount+postBinCount);
            pre_delay_post_spikeRate1_Smooth = zeros(trialNum,preBinCount+delayBinCount+postBinCount);
            pre_delay_post_spikeRateNorm1 = zeros(trialNum,preBinCount+delayBinCount+postBinCount);
            pre_delay_post_spikeRateNorm1_Smooth = zeros(trialNum,preBinCount+delayBinCount+postBinCount);
            
            
            % get each spike time
            tSp = Spike_Session.(SessDirs{j}){k}*10^3;
            
            h1 = figure(1);
            h1.Position = [100 100 1800 800];
            h2 = figure(2);
            h2.Position = [100 100 1800 800];
            % define where each subplot should be
            pos1 = [0.1 0.35 0.25 0.6];
            pos2 = [0.4 0.35 0.25 0.6];
            pos3 = [0.7 0.35 0.25 0.6];
            pos4 = [0.4 0.15 0.25 0.1];
      
            % get spike inside delay zone and calculate rate on time
            for m = 1:trialNum
                
                %% delay area
                rateBin1 = delayTstart1(m):p.timeBin:delayTend1_2(m)-p.timeBin;
                rateBin2 = delayTstart2(m):p.timeBin:delayTend2_2(m)-p.timeBin;
        
                % delay definition 1
                spikeTime1 = tSp-rateBin1(1);
                ts_Delay1 = tSp(tSp>min(rateBin1) & tSp<max(rateBin1));
                ts_Delay1 = ts_Delay1-delayTstart1(m);
                
%                 figure(1)
%                 subplot('Position',pos1)
%                 if ~isempty(ts_Delay1)
%                     plot(ts_Delay1/10^3,m,'rd')
%                 else
%                     plot(0,m)
%                 end
%                 hold on
                
                spikeCount1 = ceil(spikeTime1/p.timeBin);                
                for n = 1:length(rateBin1)
                    spikeRate1(m,n) = 10^3*sum(spikeCount1 == n)/p.timeBin;
                    %                         countNum = countNum + sum(spikeCount == n);
                end
                %                     if countNum ~= spikeNum
                %                         error('Spike count is wrong')
                %                     end
                
                % gaussian smooth the rate map
                spikeRate1_Smooth(m,:) = gaussfilt(rateBin1,spikeRate1(m,:),p.gaussSigma);                   
                spikeRateNorm1(m,:) = spikeRate1(m,:)/max(spikeRate1(m,:));
                spikeRateNorm1_Smooth(m,:) = spikeRate1_Smooth(m,:)/max(spikeRate1_Smooth(m,:));
                
                % delay definition 2
                spikeTime2 = tSp-rateBin2(1);
                ts_Delay2 = tSp(tSp>min(rateBin2) & tSp<max(rateBin2));
                ts_Delay2 = ts_Delay2-delayTstart2(m);
                
%                 figure(2)
%                 subplot('Position',pos1)
%                 if ~isempty(ts_Delay2)
%                     plot(ts_Delay2/10^3,m,'rd')
%                 else
%                     plot(0,m)
%                 end
%                 hold on
                
                spikeCount2 = ceil(spikeTime2/p.timeBin);                
                for n = 1:length(rateBin2)
                    spikeRate2(m,n) = 10^3*sum(spikeCount2 == n)/p.timeBin;
                    %                         countNum = countNum + sum(spikeCount == n);
                end
                %                     if countNum ~= spikeNum
                %                         error('Spike count is wrong')
                %                     end
                
                % gaussian smooth the rate map
                spikeRate2_Smooth(m,:) = gaussfilt(rateBin2,spikeRate2(m,:),p.gaussSigma);             
                spikeRateNorm2(m,:) = spikeRate2(m,:)/max(spikeRate2(m,:));
                spikeRateNorm2_Smooth(m,:) = spikeRate2_Smooth(m,:)/max(spikeRate2_Smooth(m,:));
                
                %% pre-delay area
                rateBin1 = delayTstart1(m)-p.preDelayTime:p.timeBin:delayTstart1(m)-p.timeBin;
                rateBin2 = delayTstart2(m)-p.preDelayTime:p.timeBin:delayTstart2(m)-p.timeBin;
        
                % delay definition 1
                spikeTime1 = tSp-rateBin1(1);
                ts_Delay1 = tSp(tSp>min(rateBin1) & tSp<max(rateBin1));
                ts_Delay1 = ts_Delay1-delayTstart1(m);
                
%                 figure(1)
%                 subplot('Position',pos1)
%                 if ~isempty(ts_Delay1)
%                     plot(ts_Delay1/10^3,m,'rd')
%                 else
%                     plot(0,m)
%                 end
%                 hold on
                
                spikeCount1 = ceil(spikeTime1/p.timeBin);                
                for n = 1:length(rateBin1)
                    pre_spikeRate1(m,n) = 10^3*sum(spikeCount1 == n)/p.timeBin;
                    %                         countNum = countNum + sum(spikeCount == n);
                end
                %                     if countNum ~= spikeNum
                %                         error('Spike count is wrong')
                %                     end
                
                % gaussian smooth the rate map
                pre_spikeRate1_Smooth(m,:) = gaussfilt(rateBin1,pre_spikeRate1(m,:),p.gaussSigma);                   
                pre_spikeRateNorm1(m,:) = pre_spikeRate1(m,:)/max(pre_spikeRate1(m,:));
                pre_spikeRateNorm1_Smooth(m,:) = pre_spikeRate1_Smooth(m,:)./max(pre_spikeRate1_Smooth(m,:));
                
                % delay definition 2
                spikeTime2 = tSp-rateBin2(1);
                ts_Delay2 = tSp(tSp>min(rateBin2) & tSp<max(rateBin2));
                ts_Delay2 = ts_Delay2-delayTstart2(m);
%                 
%                 figure(2)
%                 subplot('Position',pos1)
%                 if ~isempty(ts_Delay2)
%                     plot(ts_Delay2/10^3,m,'rd')
%                 else
%                     plot(0,m)
%                 end
%                 hold on
                
                spikeCount2 = ceil(spikeTime2/p.timeBin);                
                for n = 1:length(rateBin2)
                    pre_spikeRate2(m,n) = 10^3*sum(spikeCount2 == n)/p.timeBin;
                    %                         countNum = countNum + sum(spikeCount == n);
                end
                %                     if countNum ~= spikeNum
                %                         error('Spike count is wrong')
                %                     end
                
                % gaussian smooth the rate map
                pre_spikeRate2_Smooth(m,:) = gaussfilt(rateBin2,pre_spikeRate2(m,:),p.gaussSigma);             
                pre_spikeRateNorm2(m,:) = pre_spikeRate2(m,:)/max(pre_spikeRate2(m,:));
                pre_spikeRateNorm2_Smooth(m,:) = pre_spikeRate2_Smooth(m,:)./max(pre_spikeRate2_Smooth(m,:));
                
                %% post-delay area
                startTime = 10^3*(Fig8DelayZonePos.delayPos1.endT(m)+1/30);
%                 endTime = 10^3*PathZone.posEndT.Center(m);
                rateBin1 = startTime:p.timeBin:p.postDelayTime+startTime-p.timeBin;
        
                % delay definition 1
                spikeTime1 = tSp-rateBin1(1);
                ts_Delay1 = tSp(tSp>min(rateBin1) & tSp<max(rateBin1));
                ts_Delay1 = ts_Delay1-PathZone.posStartT.Center(m);
                

                spikeCount1 = ceil(spikeTime1/p.timeBin);                
                for n = 1:length(rateBin1)
                    post_spikeRate1(m,n) = 10^3*sum(spikeCount1 == n)/p.timeBin;
                    %                         countNum = countNum + sum(spikeCount == n);
                end
                %                     if countNum ~= spikeNum
                %                         error('Spike count is wrong')
                %                     end
                
                % gaussian smooth the rate map
                post_spikeRate1_Smooth(m,:) = gaussfilt(rateBin1,post_spikeRate1(m,:),p.gaussSigma);                   
                post_spikeRateNorm1(m,:) = post_spikeRate1(m,:)/max(post_spikeRate1(m,:));
                post_spikeRateNorm1_Smooth(m,:) = post_spikeRate1_Smooth(m,:)./max(post_spikeRate1_Smooth(m,:));
                  
                
                %% pre & delay area smooth and normalize together
                pre_delay_spikeRate1(m,:) = [pre_spikeRate1(m,:),spikeRate1(m,:)];
                pre_delay_spikeRate1_Smooth(m,:) = gaussfilt(p.timeBin*(1:length(pre_delay_spikeRate1(m,:))),pre_delay_spikeRate1(m,:),p.gaussSigma);
                pre_delay_spikeRateNorm1(m,:) = pre_delay_spikeRate1(m,:)./max(pre_delay_spikeRate1(m,:));
                pre_delay_spikeRateNorm1_Smooth(m,:) = pre_delay_spikeRate1_Smooth(m,:)./max(pre_delay_spikeRate1_Smooth(m,:));
                
                pre_delay_spikeRate2(m,:) = [pre_spikeRate2(m,:),spikeRate2(m,:)];
                pre_delay_spikeRate2_Smooth(m,:) = gaussfilt(p.timeBin*(1:length(pre_delay_spikeRate2(m,:))),pre_delay_spikeRate2(m,:),p.gaussSigma);
                pre_delay_spikeRateNorm2(m,:) = pre_delay_spikeRate2(m,:)./max(pre_delay_spikeRate2(m,:));
                pre_delay_spikeRateNorm2_Smooth(m,:) = pre_delay_spikeRate2_Smooth(m,:)./max(pre_delay_spikeRate2_Smooth(m,:));
                
                %% pre & delay & post area smooth and normalize together
                pre_delay_post_spikeRate1(m,:) = [pre_spikeRate1(m,:),spikeRate1(m,:),post_spikeRate1(m,:)];
                pre_delay_post_spikeRate1_Smooth(m,:) = gaussfilt(p.timeBin*(1:length(pre_delay_post_spikeRate1(m,:))),pre_delay_post_spikeRate1(m,:),p.gaussSigma);
                pre_delay_post_spikeRateNorm1(m,:) = pre_delay_post_spikeRate1(m,:)./max(pre_delay_post_spikeRate1(m,:));
                pre_delay_post_spikeRateNorm1_Smooth(m,:) = pre_delay_post_spikeRate1_Smooth(m,:)./max(pre_delay_post_spikeRate1_Smooth(m,:));

            end

            spikeRate1_Combined = mean(spikeRate1,1);
            spikeRate2_Combined = mean(spikeRate2,1);            
            spikeRate1_Combined_Smooth = gaussfilt(p.timeBin*(1:length(spikeRate1_Combined)),spikeRate1_Combined,p.gaussSigma);
            spikeRate2_Combined_Smooth = gaussfilt(p.timeBin*(1:length(spikeRate1_Combined)),spikeRate2_Combined,p.gaussSigma); 
            
            pre_spikeRate1_Combined = mean(pre_spikeRate1,1);
            pre_spikeRate2_Combined = mean(pre_spikeRate2,1);  
            pre_spikeRate1_Combined_Smooth = gaussfilt(p.timeBin*(1:length(pre_spikeRate1_Combined)),pre_spikeRate1_Combined,p.gaussSigma);
            pre_spikeRate2_Combined_Smooth = gaussfilt(p.timeBin*(1:length(pre_spikeRate1_Combined)),pre_spikeRate2_Combined,p.gaussSigma);   
             
            
            post_spikeRate1_Combined = mean(post_spikeRate1,1);
            post_spikeRate1_Combined_Smooth = gaussfilt(p.timeBin*(1:length(post_spikeRate1_Combined)),post_spikeRate1_Combined,p.gaussSigma);
            
            
            pre_delay_spikeRate1_Combined = mean(pre_delay_spikeRate1,1);
            pre_delay_spikeRate2_Combined = mean(pre_delay_spikeRate2,1);           
            pre_delay_spikeRate1_Combined_Smooth = gaussfilt(p.timeBin*(1:length(pre_delay_spikeRate1_Combined)),pre_delay_spikeRate1_Combined,p.gaussSigma);
            pre_delay_spikeRate2_Combined_Smooth = gaussfilt(p.timeBin*(1:length(pre_delay_spikeRate2_Combined)),pre_delay_spikeRate2_Combined,p.gaussSigma);
            
            pre_delay_post_spikeRate1_Combined = mean(pre_delay_post_spikeRate1,1);           
            pre_delay_post_spikeRate1_Combined_Smooth = gaussfilt(p.timeBin*(1:length(pre_delay_post_spikeRate1_Combined)),pre_delay_post_spikeRate1_Combined,p.gaussSigma);
%             [nFields1,FieldBinX1,fieldLabel1] = timeFieldSearch(spikeRate1_Combined_Smooth,p);
%             [nFields2,FieldBinX2,fieldLabel2] = timeFieldSearch(spikeRate2_Combined_Smooth,p); 
%             
            % consider smoothing and normalization
            % calculate delay predelay seperately and together
            DelayFire.(SessDirs{j}).spikeRate1{k} = spikeRate1;
            DelayFire.(SessDirs{j}).spikeRate1_Smooth{k} = spikeRate1_Smooth;
            DelayFire.(SessDirs{j}).spikeRateNorm1{k} = spikeRateNorm1;
            DelayFire.(SessDirs{j}).spikeRateNorm1_Smooth{k} = spikeRateNorm1_Smooth;
            DelayFire.(SessDirs{j}).spikeRate2{k} = spikeRate2;
            DelayFire.(SessDirs{j}).spikeRate2_Smooth{k} = spikeRate2_Smooth;
            DelayFire.(SessDirs{j}).spikeRateNorm2{k} = spikeRateNorm2;
            DelayFire.(SessDirs{j}).spikeRateNorm2_Smooth{k} = spikeRateNorm2_Smooth;
            DelayFire.(SessDirs{j}).spikeRate1_Combined{k} = spikeRate1_Combined;
            DelayFire.(SessDirs{j}).spikeRate2_Combined{k} = spikeRate2_Combined;
            DelayFire.(SessDirs{j}).spikeRate1_Combined_Smooth{k} = spikeRate1_Combined_Smooth;
            DelayFire.(SessDirs{j}).spikeRate2_Combined_Smooth{k} = spikeRate2_Combined_Smooth;
            
            DelayFire.(SessDirs{j}).pre_spikeRate1{k} = pre_spikeRate1;
            DelayFire.(SessDirs{j}).pre_spikeRate1_Smooth{k} = pre_spikeRate1_Smooth;
            DelayFire.(SessDirs{j}).pre_spikeRateNorm1{k} = pre_spikeRateNorm1;
            DelayFire.(SessDirs{j}).pre_spikeRateNorm1_Smooth{k} = pre_spikeRateNorm1_Smooth;
            DelayFire.(SessDirs{j}).pre_spikeRate2{k} = pre_spikeRate2;
            DelayFire.(SessDirs{j}).pre_spikeRate2_Smooth{k} = pre_spikeRate2_Smooth;
            DelayFire.(SessDirs{j}).pre_spikeRateNorm2{k} = pre_spikeRateNorm2;
            DelayFire.(SessDirs{j}).pre_spikeRateNorm2_Smooth{k} = pre_spikeRateNorm2_Smooth;
            DelayFire.(SessDirs{j}).pre_spikeRate1_Combined{k} = pre_spikeRate1_Combined;
            DelayFire.(SessDirs{j}).pre_spikeRate2_Combined{k} = pre_spikeRate2_Combined;
            DelayFire.(SessDirs{j}).pre_spikeRate1_Combined_Smooth{k} = pre_spikeRate1_Combined_Smooth;
            DelayFire.(SessDirs{j}).pre_spikeRate2_Combined_Smooth{k} = pre_spikeRate2_Combined_Smooth;
            
            DelayFire.(SessDirs{j}).post_spikeRate1{k} = post_spikeRate1;
            DelayFire.(SessDirs{j}).post_spikeRate1_Smooth{k} = post_spikeRate1_Smooth;
            DelayFire.(SessDirs{j}).post_spikeRateNorm1{k} = post_spikeRateNorm1;
            DelayFire.(SessDirs{j}).post_spikeRateNorm1_Smooth{k} = post_spikeRateNorm1_Smooth;
            DelayFire.(SessDirs{j}).post_spikeRate1_Combined{k} = post_spikeRate1_Combined;
            DelayFire.(SessDirs{j}).post_spikeRate1_Combined_Smooth{k} = post_spikeRate1_Combined_Smooth;
            
            
            DelayFire.(SessDirs{j}).pre_delay_spikeRate1{k} = pre_delay_spikeRate1;
            DelayFire.(SessDirs{j}).pre_delay_spikeRate1_Smooth{k} = pre_delay_spikeRate1_Smooth;
            DelayFire.(SessDirs{j}).pre_delay_spikeRateNorm1{k} = pre_delay_spikeRateNorm1;
            DelayFire.(SessDirs{j}).pre_delay_spikeRateNorm1_Smooth{k} = pre_delay_spikeRateNorm1_Smooth;
            DelayFire.(SessDirs{j}).pre_delay_spikeRate2{k} = pre_delay_spikeRate2;
            DelayFire.(SessDirs{j}).pre_delay_spikeRate2_Smooth{k} = pre_delay_spikeRate2_Smooth;
            DelayFire.(SessDirs{j}).pre_delay_spikeRateNorm2{k} = pre_delay_spikeRateNorm2;
            DelayFire.(SessDirs{j}).pre_delay_spikeRateNorm2_Smooth{k} = pre_delay_spikeRateNorm2_Smooth;
            DelayFire.(SessDirs{j}).pre_delay_spikeRate1_Combined{k} = pre_delay_spikeRate1_Combined;
            DelayFire.(SessDirs{j}).pre_delay_spikeRate2_Combined{k} = pre_delay_spikeRate2_Combined;
            DelayFire.(SessDirs{j}).pre_delay_spikeRate1_Combined_Smooth{k} = pre_delay_spikeRate1_Combined_Smooth;
            DelayFire.(SessDirs{j}).pre_delay_spikeRate2_Combined_Smooth{k} = pre_delay_spikeRate2_Combined_Smooth;
            
            DelayFire.(SessDirs{j}).pre_delay_post_post_spikeRate1{k} = pre_delay_post_spikeRate1;
            DelayFire.(SessDirs{j}).pre_delay_post_spikeRate1_Smooth{k} = pre_delay_post_spikeRate1_Smooth;
            DelayFire.(SessDirs{j}).pre_delay_post_spikeRateNorm1{k} = pre_delay_post_spikeRateNorm1;
            DelayFire.(SessDirs{j}).pre_delay_post_spikeRateNorm1_Smooth{k} = pre_delay_post_spikeRateNorm1_Smooth;
            DelayFire.(SessDirs{j}).pre_delay_post_spikeRate1_Combined{k} = pre_delay_post_spikeRate1_Combined;
            DelayFire.(SessDirs{j}).pre_delay_post_spikeRate1_Combined_Smooth{k} = pre_delay_post_spikeRate1_Combined_Smooth;
            
            
%             figure(1)
%             subplot('Position',pos1)
%             TITLE1 = sprintf('%s%d%s%d%s%s%s%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k}(1:end-2),'-',SessDirs{j},'-Spikes');
%             title({TITLE1},'Interpreter','None')
%             xlabel('Time (Sec)')
%             ylabel('Trials')
%             ylim([0.5 trialNum+0.5])
%             xlim([-p.preDelayTime/1000 maxT/1000])
%             xTick= [-p.preDelayTime/1000 0 maxT/1000];
%             xTickLabel = [-p.preDelayTime/1000 0 maxT/1000];
%             set (gca,'YDir','reverse')
%             set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
%             
%             
%             subplot('Position',pos2)
%             imagesc(pre_delay_spikeRate1_Smooth)
%             colormap(jet)
%             TITLE1 = 'RateMap - Onset at barrier';
% %             TITLE2 = sessInfo(i).mainDir;
%             title({TITLE1},'Interpreter','None')
%             xlabel('Time (Sec)')
%             ylabel('Trials')
%             xlim([0 delayBinCount+preBinCount])
%             xTick= [0 preBinCount delayBinCount+preBinCount];
%             xTickLabel = [-p.preDelayTime/1000 0 maxT/1000];    
%             set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
%             
%             subplot('Position',pos3)
%             imagesc(pre_delay_spikeRateNorm1_Smooth)
%             colormap(jet)     
% %             TITLE1 = sprintf('%s%d%s%d%s%s%s%s','DelayZone-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k}(1:end-2),'-',SessDirs{j});
%             TITLE2 = 'Normalized rate';
%             title({TITLE2},'Interpreter','None')
%             xlabel('Time (Sec)')
%             ylabel('Trials')
%             xlim([0 delayBinCount+preBinCount])
%             xTick= [0 preBinCount delayBinCount+preBinCount];
%             xTickLabel = [-p.preDelayTime/1000 0 maxT/1000];         
%             set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
%                                    
%             subplot('Position',pos4)
%             imagesc(pre_delay_spikeRate1_Combined_Smooth)
%             colormap(jet)
%             title('Combined time map')
%             xlim([0 delayBinCount+preBinCount])
%             xTick= [0 preBinCount delayBinCount+preBinCount];
%             xTickLabel = [-p.preDelayTime/1000 0 maxT/1000];          
%             set(gca, 'ytick',[],'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
%             
%             
%             figure(2)
%             subplot('Position',pos1)
%             TITLE1 = sprintf('%s%d%s%d%s%s%s%s%s','DelayZone2-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k}(1:end-2),'-',SessDirs{j},'-Spikes');
%             title({TITLE1},'Interpreter','None')
%             xlabel('Time (Sec)')
%             ylabel('Trials')
%             ylim([0.5 trialNum+0.5])
%             xlim([-p.preDelayTime/1000 maxT/1000])
%             xTick= [-p.preDelayTime/1000 0 maxT/1000];
%             xTickLabel = [-p.preDelayTime/1000 0 maxT/1000];
%             set (gca,'YDir','reverse')
%             set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
%             
%             subplot('Position',pos2)
%             imagesc(pre_delay_spikeRate2_Smooth)
%             colormap(jet)
%             TITLE1 = 'RateMap - Onset at entrance';
% %             TITLE2 = sessInfo(i).mainDir;
%             title({TITLE1},'Interpreter','None')
%             xlabel('Time (Sec)')
%             ylabel('Trials')
%             xlim([0 delayBinCount+preBinCount])
%             xTick= [0 preBinCount delayBinCount+preBinCount];
%             xTickLabel = [-p.preDelayTime/1000 0 maxT/1000];     
%             set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
%             
%             subplot('Position',pos3)
%             imagesc(pre_delay_spikeRateNorm2_Smooth)
%             colormap(jet)
%             
% %             TITLE1 = sprintf('%s%d%s%d%s%s%s%s','DelayZone-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k}(1:end-2),'-',SessDirs{j});
%             TITLE2 = 'Normalized rate';
%             title({TITLE2},'Interpreter','None')
%             xlabel('Time (Sec)')
%             ylabel('Trials')
%             xlim([0 delayBinCount+preBinCount])
%             xTick= [0 preBinCount delayBinCount+preBinCount];
%             xTickLabel = [-p.preDelayTime/1000 0 maxT/1000];          
%             set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);          
%             
%             subplot('Position',pos4)
%             imagesc(pre_delay_spikeRate2_Combined_Smooth)
%             colormap(jet)
%             title('Combined time map')
%             xlim([0 delayBinCount+preBinCount])
%             xTick= [0 preBinCount delayBinCount+preBinCount];
%             xTickLabel = [-p.preDelayTime/1000 0 maxT/1000];        
%             set(gca, 'ytick',[],'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]); 
            
            if p.savePlot == 1
                figure(1)
                figName = sprintf('%s%s%d%s%d%s%s%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',TList{k}(1:end-2),'-',SessDirs{j},'-DelayFiring_def1');
                print(figName,'-dpng','-r300');            
                figure(2)
                figName = sprintf('%s%s%d%s%d%s%s%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',TList{k}(1:end-2),'-',SessDirs{j},'-DelayFiring_def2');
                print(figName,'-dpng','-r300');
            end          
            
            close all
        end
    end
    if p.writeToFile == 1
        save(fullfile(savedir2,'Fig8DelayTimeMap.mat'), 'DelayFire');
    end
    clear DelayFire
    fprintf('Finished position analysis for session %d\n',i);
end

end