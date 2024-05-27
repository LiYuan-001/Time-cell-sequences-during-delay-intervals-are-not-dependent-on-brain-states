function Fig8TreadmillDelayFiring_Assembly(inFile,AnalyzeSes)

p.savePlot = 1;
p.writeToFile = 1;

p.timeBin = 150; % ms;
p.gaussSigma = p.timeBin*2; % sd of the gaussian distribution
p.preDelayTime = 3000; % unit ms
p.postDelayTime = 3000; % unit ms

% Read in input information
sessInfo = SessInfoImport(inFile);
close all

for i = AnalyzeSes(1:end)
    
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\Delay_Assembly_TimeMap_onoff');
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
    DelayFire_Assembly_onoff.rat = sessInfo(i).animal;
    DelayFire_Assembly_onoff.day = sessInfo(i).day;
    DelayFire_Assembly_onoff.timeBin = p.timeBin;
    DelayFire_Assembly_onoff.gaussSigma = p.gaussSigma;
             
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load assembly file
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_DelayOnOff_LR-25ms.mat');
    load(assemblyFile);
    
    % get time stamps for the assembly
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;    
    tList = SpikeProp.tList(rateLabel);
    
    patNum_on = CellAssembly_DelayLR.DelayOn.patNum ;
    AssmblPtrnCellIDs_on = CellAssembly_DelayLR.DelayOn.AssmblPtrnCellIDs;
    AssmblWght_on = CellAssembly_DelayLR.DelayOn.AssmblWght;
    AssmblStrength_on = CellAssembly_DelayLR.DelayOn.AssmblStrength;
    event_Time_on = CellAssembly_DelayLR.DelayOn.event_Time;
    event_strength_on = CellAssembly_DelayLR.DelayOn.event_strength;
    event_Num_on = CellAssembly_DelayLR.DelayOn.event_Num;
    
    patNum_off = CellAssembly_DelayLR.DelayOff.patNum ;
    AssmblPtrnCellIDs_off = CellAssembly_DelayLR.DelayOff.AssmblPtrnCellIDs;
    AssmblWght_off = CellAssembly_DelayLR.DelayOff.AssmblWght;
    AssmblStrength_off = CellAssembly_DelayLR.DelayOff.AssmblStrength;
    event_Time_off = CellAssembly_DelayLR.DelayOff.event_Time;
    event_strength_off = CellAssembly_DelayLR.DelayOff.event_strength;
    event_Num_off = CellAssembly_DelayLR.DelayOff.event_Num;
    
    for j = 1:length(sessDirs)
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'PathZone.mat');
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
        if contains(sessDirs{j},'10')
            maxT = 10*10^3;
            delayTend1_2 = delayTstart1+maxT;
            delayTend2_2 = delayTstart2+maxT;
            delayBinCount = floor(maxT/p.timeBin);
        elseif contains(sessDirs{j},'30')
            maxT = 30*10^3;
            delayTend1_2 = delayTstart1+maxT;
            delayTend2_2 = delayTstart2+maxT;
            delayBinCount = floor(maxT/p.timeBin);
        else
            error('Delay time is wrong')
        end       
        preBinCount = floor(p.preDelayTime/p.timeBin);
        postBinCount = floor(p.postDelayTime/p.timeBin);
        
        DelayFire_Assembly_onoff.(sessDirs{j}).delayTstart1 = delayTstart1;
%         DelayFire_Assembly_onoff.(sessDirs{j}).delayTstart2 = delayTstart2;
        DelayFire_Assembly_onoff.(sessDirs{j}).delayTend1 = delayTend1_2;
%         DelayFire_Assembly_onoff.(sessDirs{j}).delayTend2 = delayTend2_2;
         
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
        
        for k = 1:patNum_on
            % delay area
            % event time as spike
            spikeRate1 = zeros(trialNum,delayBinCount);
            spikeRate1_Smooth = zeros(trialNum,delayBinCount);
            spikeRateNorm1 = zeros(trialNum,delayBinCount);
            spikeRateNorm1_Smooth = zeros(trialNum,delayBinCount);
            
            % event strength
            spikeRate2 = zeros(trialNum,delayBinCount);
            spikeRate2_Smooth = zeros(trialNum,delayBinCount);
            spikeRateNorm2 = zeros(trialNum,delayBinCount);
            spikeRateNorm2_Smooth = zeros(trialNum,delayBinCount);
            
            % pre area
            % event time as spike
            pre_spikeRate1 = zeros(trialNum,preBinCount);
            pre_spikeRate1_Smooth = zeros(trialNum,preBinCount);
            pre_spikeRateNorm1 = zeros(trialNum,preBinCount);
            pre_spikeRateNorm1_Smooth = zeros(trialNum,preBinCount);
            
            % event strength
            pre_spikeRate2 = zeros(trialNum,preBinCount);
            pre_spikeRate2_Smooth = zeros(trialNum,preBinCount);
            pre_spikeRateNorm2 = zeros(trialNum,preBinCount);
            pre_spikeRateNorm2_Smooth = zeros(trialNum,preBinCount);
            
            % post delay
            post_spikeRate1 = zeros(trialNum,postBinCount);
            post_spikeRate1_Smooth = zeros(trialNum,postBinCount);
            post_spikeRateNorm1 = zeros(trialNum,postBinCount);
            post_spikeRateNorm1_Smooth = zeros(trialNum,postBinCount);
            
            post_spikeRate2 = zeros(trialNum,postBinCount);
            post_spikeRate2_Smooth = zeros(trialNum,postBinCount);
            post_spikeRateNorm2 = zeros(trialNum,postBinCount);
            post_spikeRateNorm2_Smooth = zeros(trialNum,postBinCount);
            
            pre_delay_post_spikeRate1 = zeros(trialNum,preBinCount + postBinCount + delayBinCount);
            pre_delay_post_spikeRate1_Smooth = zeros(trialNum,preBinCount + postBinCount + delayBinCount);
            pre_delay_post_spikeRateNorm1 = zeros(trialNum,preBinCount + postBinCount + delayBinCount);
            pre_delay_post_spikeRateNorm1_Smooth = zeros(trialNum,preBinCount + postBinCount + delayBinCount);
            
            pre_delay_post_spikeRate2 = zeros(trialNum,preBinCount + postBinCount + delayBinCount);
            pre_delay_post_spikeRate2_Smooth = zeros(trialNum,preBinCount + postBinCount + delayBinCount);
            pre_delay_post_spikeRateNorm2 = zeros(trialNum,preBinCount + postBinCount + delayBinCount);
            pre_delay_post_spikeRateNorm2_Smooth = zeros(trialNum,preBinCount + postBinCount + delayBinCount);
                
            % get each spike time
            tSp = event_Time_on{k}*10^3;
            tSp_2 = event_strength_on{k};
            
            h1 = figure(1);
            h1.Position = [100 100 900 600];
            % define where each subplot should be
            pos1 = [0.1 0.35 0.3 0.6];
            pos2 = [0.1 0.15 0.3 0.1];
            pos3 = [0.6 0.35 0.3 0.6];
            pos4 = [0.6 0.15 0.3 0.1];
      
            % get spike inside delay zone and calculate rate on time
            for m = 1:trialNum
                rateBin1 = delayTstart1(m):p.timeBin:delayTend1_2(m)-p.timeBin;
%                 rateBin2 = delayTstart2(m):p.timeBin:delayTend2_2(m)-p.timeBin;
        
                % event time map
                spikeTime1 = tSp-rateBin1(1);
                ts_Delay1 = tSp(tSp>min(rateBin1) & tSp<max(rateBin1));
                ts_Delay1 = ts_Delay1-min(rateBin1);
                
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
                spikeRateNorm1_Smooth(m,:) = gaussfilt(rateBin1,spikeRateNorm1(m,:),p.gaussSigma);

                
                % event strength
                for n = 1:length(rateBin1)
                    ts_Delay2 = sum(tSp_2(tSp>rateBin1(n) & tSp<(rateBin1(n)+p.timeBin)));
                    spikeRate2(m,n) = 10^3*ts_Delay2/p.timeBin;
                    %                         countNum = countNum + sum(spikeCount == n);
                end
                %                     if countNum ~= spikeNum
                %                         error('Spike count is wrong')
                %                     end
                
                % gaussian smooth the rate map
                spikeRate2_Smooth(m,:) = gaussfilt(rateBin1,spikeRate2(m,:),p.gaussSigma);             
                spikeRateNorm2(m,:) = spikeRate2(m,:)/max(spikeRate2(m,:));
                spikeRateNorm2_Smooth(m,:) = gaussfilt(rateBin1,spikeRateNorm2(m,:),p.gaussSigma);
                
                                %% pre-delay area
                rateBin1 = delayTstart1(m)-p.preDelayTime:p.timeBin:delayTstart1(m)-p.timeBin;
        
                % delay definition 1
                spikeTime1 = tSp-rateBin1(1);
                
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
                
                % event strength          
                for n = 1:length(rateBin1)
                    ts_Delay2 = sum(tSp_2(tSp>rateBin1(n) & tSp<(rateBin1(n)+p.timeBin)));
                    pre_spikeRate2(m,n) = 10^3*ts_Delay2/p.timeBin;
                end

                % gaussian smooth the rate map
                pre_spikeRate2_Smooth(m,:) = gaussfilt(rateBin1,pre_spikeRate2(m,:),p.gaussSigma);             
                pre_spikeRateNorm2(m,:) = pre_spikeRate2(m,:)/max(pre_spikeRate2(m,:));
                pre_spikeRateNorm2_Smooth(m,:) = pre_spikeRate2_Smooth(m,:)./max(pre_spikeRate2_Smooth(m,:));
                
                % post-delay area
                % event time
                startTime = 10^3*(Fig8DelayZonePos.delayPos1.endT(m)+1/30);
%                 endTime = 10^3*PathZone.posEndT.Center(m);
                rateBin1 = startTime:p.timeBin:p.postDelayTime+startTime-p.timeBin;
        
                % event time
                spikeTime1 = tSp-rateBin1(1);
                ts_Delay1 = tSp(tSp>min(rateBin1) & tSp<max(rateBin1));
                ts_Delay1 = ts_Delay1-PathZone.posStartT.Center(m);

                spikeCount1 = ceil(spikeTime1/p.timeBin);                
                for n = 1:length(rateBin1)
                    post_spikeRate1(m,n) = 10^3*sum(spikeCount1 == n)/p.timeBin;
                end
                
                % gaussian smooth the rate map
                post_spikeRate1_Smooth(m,:) = gaussfilt(rateBin1,post_spikeRate1(m,:),p.gaussSigma);                   
                post_spikeRateNorm1(m,:) = post_spikeRate1(m,:)/max(post_spikeRate1(m,:));
                post_spikeRateNorm1_Smooth(m,:) = post_spikeRate1_Smooth(m,:)./max(post_spikeRate1_Smooth(m,:));
                  
                % event strength          
                for n = 1:length(rateBin1)
                    ts_Delay2 = sum(tSp_2(tSp>rateBin1(n) & tSp<(rateBin1(n)+p.timeBin)));
                    post_spikeRate2(m,n) = 10^3*ts_Delay2/p.timeBin;
                end

                % gaussian smooth the rate map
                post_spikeRate2_Smooth(m,:) = gaussfilt(rateBin1,post_spikeRate2(m,:),p.gaussSigma);             
                post_spikeRateNorm2(m,:) = post_spikeRate2(m,:)/max(post_spikeRate2(m,:));
                post_spikeRateNorm2_Smooth(m,:) = post_spikeRate2_Smooth(m,:)./max(post_spikeRate2_Smooth(m,:));
                
                
                %% pre & delay & post area smooth and normalize together
                pre_delay_post_spikeRate1(m,:) = [pre_spikeRate1(m,:),spikeRate1(m,:),post_spikeRate1(m,:)];
                pre_delay_post_spikeRate1_Smooth(m,:) = gaussfilt(p.timeBin*(1:length(pre_delay_post_spikeRate1(m,:))),pre_delay_post_spikeRate1(m,:),p.gaussSigma);
                pre_delay_post_spikeRateNorm1(m,:) = pre_delay_post_spikeRate1(m,:)./max(pre_delay_post_spikeRate1(m,:));
                pre_delay_post_spikeRateNorm1_Smooth(m,:) = pre_delay_post_spikeRate1_Smooth(m,:)./max(pre_delay_post_spikeRate1_Smooth(m,:));
                
                pre_delay_post_spikeRate2(m,:) = [pre_spikeRate2(m,:),spikeRate2(m,:),post_spikeRate2(m,:)];
                pre_delay_post_spikeRate2_Smooth(m,:) = gaussfilt(p.timeBin*(1:length(pre_delay_post_spikeRate2(m,:))),pre_delay_post_spikeRate2(m,:),p.gaussSigma);
                pre_delay_post_spikeRateNorm2(m,:) = pre_delay_post_spikeRate2(m,:)./max(pre_delay_post_spikeRate2(m,:));
                pre_delay_post_spikeRateNorm2_Smooth(m,:) = pre_delay_post_spikeRate2_Smooth(m,:)./max(pre_delay_post_spikeRate2_Smooth(m,:));
            end

            spikeRate1_Combined = mean(spikeRate1,1);
            spikeRate2_Combined = mean(spikeRate2,1);
            
            spikeRate1_Combined_Smooth = gaussfilt((1:length(spikeRate1_Combined))*p.timeBin,spikeRate1_Combined,p.gaussSigma);
            spikeRate2_Combined_Smooth = gaussfilt((1:length(spikeRate1_Combined))*p.timeBin,spikeRate2_Combined,p.gaussSigma);   
            
            pre_delay_post_spikeRate1_Combined = mean(pre_delay_post_spikeRate1,1);
            pre_delay_post_spikeRate2_Combined = mean(pre_delay_post_spikeRate2,1);
            
            pre_delay_post_spikeRate1_Combined_Smooth = gaussfilt((1:length(pre_delay_post_spikeRate1_Combined))*p.timeBin,pre_delay_post_spikeRate1_Combined,p.gaussSigma);
            pre_delay_post_spikeRate2_Combined_Smooth = gaussfilt((1:length(pre_delay_post_spikeRate1_Combined))*p.timeBin,pre_delay_post_spikeRate2_Combined,p.gaussSigma);   
            
            
            % fix all NaNs to 0;
            % NaNs means the trial with all 0 Hz bins in this code
            spikeRate1(isnan(spikeRate1)) = 0;
            spikeRate1_Smooth(isnan(spikeRate1_Smooth)) = 0;
            spikeRateNorm1(isnan(spikeRateNorm1)) = 0;
            spikeRateNorm1_Smooth(isnan(spikeRateNorm1_Smooth)) = 0;            
            spikeRate2(isnan(spikeRate2)) = 0;
            spikeRate2_Smooth(isnan(spikeRate2_Smooth)) = 0;
            spikeRateNorm2(isnan(spikeRateNorm2)) = 0;
            spikeRateNorm2_Smooth(isnan(spikeRateNorm2_Smooth)) = 0;            
            spikeRate1_Combined(isnan(spikeRate1_Combined)) = 0;
            spikeRate2_Combined(isnan(spikeRate2_Combined)) = 0;
            spikeRate1_Combined_Smooth(isnan(spikeRate1_Combined_Smooth)) = 0;
            spikeRate2_Combined_Smooth(isnan(spikeRate2_Combined_Smooth)) = 0;
            
            DelayFire_Assembly_onoff.on.(sessDirs{j}).spikeRate1{k} = spikeRate1;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).spikeRate1_Smooth{k} = spikeRate1_Smooth;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).spikeRateNorm1{k} = spikeRateNorm1;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).spikeRateNorm1_Smooth{k} = spikeRateNorm1_Smooth;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).spikeRate2{k} = spikeRate2;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).spikeRate2_Smooth{k} = spikeRate2_Smooth;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).spikeRateNorm2{k} = spikeRateNorm2;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).spikeRateNorm2_Smooth{k} = spikeRateNorm2_Smooth;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).spikeRate1_Combined{k} = spikeRate1_Combined;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).spikeRate2_Combined{k} = spikeRate2_Combined;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).spikeRate1_Combined_Smooth{k} = spikeRate1_Combined_Smooth;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).spikeRate2_Combined_Smooth{k} = spikeRate2_Combined_Smooth;
             
            DelayFire_Assembly_onoff.on.(sessDirs{j}).pre_spikeRate1{k} = pre_spikeRate1;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).pre_spikeRate1_Smooth{k} = pre_spikeRate1_Smooth;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).pre_spikeRateNorm1{k} = pre_spikeRateNorm1;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).pre_spikeRateNorm1_Smooth{k} = pre_spikeRateNorm1_Smooth;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).pre_spikeRate2{k} = pre_spikeRate2;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).pre_spikeRate2_Smooth{k} = pre_spikeRate2_Smooth;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).pre_spikeRateNorm2{k} = pre_spikeRateNorm2;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).pre_spikeRateNorm2_Smooth{k} = pre_spikeRateNorm2_Smooth;

            DelayFire_Assembly_onoff.on.(sessDirs{j}).post_spikeRate1{k} = post_spikeRate1;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).post_spikeRate1_Smooth{k} = post_spikeRate1_Smooth;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).post_spikeRateNorm1{k} = post_spikeRateNorm1;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).post_spikeRateNorm1_Smooth{k} = post_spikeRateNorm1_Smooth;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).post_spikeRate2{k} = post_spikeRate2;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).post_spikeRate2_Smooth{k} = post_spikeRate2_Smooth;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).post_spikeRateNorm2{k} = post_spikeRateNorm2;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).post_spikeRateNorm2_Smooth{k} = post_spikeRateNorm2_Smooth;
            
            DelayFire_Assembly_onoff.on.(sessDirs{j}).pre_delay_post_spikeRate1_Combined{k} = pre_delay_post_spikeRate1_Combined;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).pre_delay_post_spikeRate2_Combined{k} = pre_delay_post_spikeRate2_Combined;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).pre_delay_post_spikeRate1_Combined_Smooth{k} = pre_delay_post_spikeRate1_Combined_Smooth;
            DelayFire_Assembly_onoff.on.(sessDirs{j}).pre_delay_post_spikeRate2_Combined_Smooth{k} = pre_delay_post_spikeRate2_Combined_Smooth;
            
            figure(1)            
            subplot('Position',pos1)
            imagesc([pre_spikeRate1_Smooth,spikeRate1_Smooth,post_spikeRate1_Smooth])
            colormap(jet)
            TITLE1 = sprintf('%s%d%s%d%s%d%s%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-on_Assembly-',k,'-',sessDirs{j},'-eventTime');
            title({TITLE1},'Interpreter','None')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
            ylabel('Trials')
%             xlim([0 delayBinCount])
                       
            subplot('Position',pos2)
            imagesc(pre_delay_post_spikeRate1_Combined_Smooth)
            colormap(jet)     
            title('Combined time map')
            xlabel('Time (Sec)')
            ylabel('Trials')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
                                   
            subplot('Position',pos3)
            imagesc([pre_spikeRate2_Smooth,spikeRate2_Smooth,post_spikeRate2_Smooth])
            colormap(jet)
            title('Event strength','Interpreter','None')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
            ylabel('Trials')
%             xlim([0 delayBinCount])
                       
            subplot('Position',pos4)
            imagesc(pre_delay_post_spikeRate2_Combined_Smooth)
            colormap(jet)     
            title('Combined time map')
            xlabel('Time (Sec)')
            ylabel('Trials')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
            
            if p.savePlot == 1
                figure(1)
                figName = sprintf('%s%s%d%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-OnDelayAssembly-',k,'-',sessDirs{j},'-DelayFiring');
                print(figName,'-dpng','-r300');            
            end          
            close all
        end
        
        
        for k = 1:patNum_off
            % delay area
            % event time as spike
            spikeRate1 = zeros(trialNum,delayBinCount);
            spikeRate1_Smooth = zeros(trialNum,delayBinCount);
            spikeRateNorm1 = zeros(trialNum,delayBinCount);
            spikeRateNorm1_Smooth = zeros(trialNum,delayBinCount);
            
            % event strength
            spikeRate2 = zeros(trialNum,delayBinCount);
            spikeRate2_Smooth = zeros(trialNum,delayBinCount);
            spikeRateNorm2 = zeros(trialNum,delayBinCount);
            spikeRateNorm2_Smooth = zeros(trialNum,delayBinCount);
            
            % pre area
            % event time as spike
            pre_spikeRate1 = zeros(trialNum,preBinCount);
            pre_spikeRate1_Smooth = zeros(trialNum,preBinCount);
            pre_spikeRateNorm1 = zeros(trialNum,preBinCount);
            pre_spikeRateNorm1_Smooth = zeros(trialNum,preBinCount);
            
            % event strength
            pre_spikeRate2 = zeros(trialNum,preBinCount);
            pre_spikeRate2_Smooth = zeros(trialNum,preBinCount);
            pre_spikeRateNorm2 = zeros(trialNum,preBinCount);
            pre_spikeRateNorm2_Smooth = zeros(trialNum,preBinCount);
            
            % post delay
            post_spikeRate1 = zeros(trialNum,postBinCount);
            post_spikeRate1_Smooth = zeros(trialNum,postBinCount);
            post_spikeRateNorm1 = zeros(trialNum,postBinCount);
            post_spikeRateNorm1_Smooth = zeros(trialNum,postBinCount);
            
            post_spikeRate2 = zeros(trialNum,postBinCount);
            post_spikeRate2_Smooth = zeros(trialNum,postBinCount);
            post_spikeRateNorm2 = zeros(trialNum,postBinCount);
            post_spikeRateNorm2_Smooth = zeros(trialNum,postBinCount);
            
            % get each spike time
            tSp = event_Time_off{k}*10^3;
            tSp_2 = event_strength_off{k};
            
            h1 = figure(1);
            h1.Position = [100 100 900 600];
            % define where each subplot should be
            pos1 = [0.1 0.35 0.3 0.6];
            pos2 = [0.1 0.15 0.3 0.1];
            pos3 = [0.6 0.35 0.3 0.6];
            pos4 = [0.6 0.15 0.3 0.1];
      
            % get spike inside delay zone and calculate rate on time
            for m = 1:trialNum
                rateBin1 = delayTstart1(m):p.timeBin:delayTend1_2(m)-p.timeBin;
%                 rateBin2 = delayTstart2(m):p.timeBin:delayTend2_2(m)-p.timeBin;
        
                % event time map
                spikeTime1 = tSp-rateBin1(1);
                ts_Delay1 = tSp(tSp>min(rateBin1) & tSp<max(rateBin1));
                ts_Delay1 = ts_Delay1-min(rateBin1);
                
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
                spikeRateNorm1_Smooth(m,:) = gaussfilt(rateBin1,spikeRateNorm1(m,:),p.gaussSigma);

                
                % event strength
                for n = 1:length(rateBin1)
                    ts_Delay2 = sum(tSp_2(tSp>rateBin1(n) & tSp<(rateBin1(n)+p.timeBin)));
                    spikeRate2(m,n) = 10^3*ts_Delay2/p.timeBin;
                    %                         countNum = countNum + sum(spikeCount == n);
                end
                %                     if countNum ~= spikeNum
                %                         error('Spike count is wrong')
                %                     end
                
                % gaussian smooth the rate map
                spikeRate2_Smooth(m,:) = gaussfilt(rateBin1,spikeRate2(m,:),p.gaussSigma);             
                spikeRateNorm2(m,:) = spikeRate2(m,:)/max(spikeRate2(m,:));
                spikeRateNorm2_Smooth(m,:) = gaussfilt(rateBin1,spikeRateNorm2(m,:),p.gaussSigma);
                
                                %% pre-delay area
                rateBin1 = delayTstart1(m)-p.preDelayTime:p.timeBin:delayTstart1(m)-p.timeBin;
        
                % delay definition 1
                spikeTime1 = tSp-rateBin1(1);
                
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
                
                % event strength          
                for n = 1:length(rateBin1)
                    ts_Delay2 = sum(tSp_2(tSp>rateBin1(n) & tSp<(rateBin1(n)+p.timeBin)));
                    pre_spikeRate2(m,n) = 10^3*ts_Delay2/p.timeBin;
                end

                % gaussian smooth the rate map
                pre_spikeRate2_Smooth(m,:) = gaussfilt(rateBin1,pre_spikeRate2(m,:),p.gaussSigma);             
                pre_spikeRateNorm2(m,:) = pre_spikeRate2(m,:)/max(pre_spikeRate2(m,:));
                pre_spikeRateNorm2_Smooth(m,:) = pre_spikeRate2_Smooth(m,:)./max(pre_spikeRate2_Smooth(m,:));
                
                % post-delay area
                % event time
                startTime = 10^3*(Fig8DelayZonePos.delayPos1.endT(m)+1/30);
%                 endTime = 10^3*PathZone.posEndT.Center(m);
                rateBin1 = startTime:p.timeBin:p.postDelayTime+startTime-p.timeBin;
        
                % event time
                spikeTime1 = tSp-rateBin1(1);
                ts_Delay1 = tSp(tSp>min(rateBin1) & tSp<max(rateBin1));
                ts_Delay1 = ts_Delay1-PathZone.posStartT.Center(m);

                spikeCount1 = ceil(spikeTime1/p.timeBin);                
                for n = 1:length(rateBin1)
                    post_spikeRate1(m,n) = 10^3*sum(spikeCount1 == n)/p.timeBin;
                end
                
                % gaussian smooth the rate map
                post_spikeRate1_Smooth(m,:) = gaussfilt(rateBin1,post_spikeRate1(m,:),p.gaussSigma);                   
                post_spikeRateNorm1(m,:) = post_spikeRate1(m,:)/max(post_spikeRate1(m,:));
                post_spikeRateNorm1_Smooth(m,:) = post_spikeRate1_Smooth(m,:)./max(post_spikeRate1_Smooth(m,:));
                  
                % event strength          
                for n = 1:length(rateBin1)
                    ts_Delay2 = sum(tSp_2(tSp>rateBin1(n) & tSp<(rateBin1(n)+p.timeBin)));
                    post_spikeRate2(m,n) = 10^3*ts_Delay2/p.timeBin;
                end

                % gaussian smooth the rate map
                post_spikeRate2_Smooth(m,:) = gaussfilt(rateBin1,post_spikeRate2(m,:),p.gaussSigma);             
                post_spikeRateNorm2(m,:) = post_spikeRate2(m,:)/max(post_spikeRate2(m,:));
                post_spikeRateNorm2_Smooth(m,:) = post_spikeRate2_Smooth(m,:)./max(post_spikeRate2_Smooth(m,:));
                
                
                %% pre & delay & post area smooth and normalize together
                pre_delay_post_spikeRate1(m,:) = [pre_spikeRate1(m,:),spikeRate1(m,:),post_spikeRate1(m,:)];
                pre_delay_post_spikeRate1_Smooth(m,:) = gaussfilt(p.timeBin*(1:length(pre_delay_post_spikeRate1(m,:))),pre_delay_post_spikeRate1(m,:),p.gaussSigma);
                pre_delay_post_spikeRateNorm1(m,:) = pre_delay_post_spikeRate1(m,:)./max(pre_delay_post_spikeRate1(m,:));
                pre_delay_post_spikeRateNorm1_Smooth(m,:) = pre_delay_post_spikeRate1_Smooth(m,:)./max(pre_delay_post_spikeRate1_Smooth(m,:));
                
                pre_delay_post_spikeRate2(m,:) = [pre_spikeRate2(m,:),spikeRate2(m,:),post_spikeRate2(m,:)];
                pre_delay_post_spikeRate2_Smooth(m,:) = gaussfilt(p.timeBin*(1:length(pre_delay_post_spikeRate2(m,:))),pre_delay_post_spikeRate2(m,:),p.gaussSigma);
                pre_delay_post_spikeRateNorm2(m,:) = pre_delay_post_spikeRate2(m,:)./max(pre_delay_post_spikeRate2(m,:));
                pre_delay_post_spikeRateNorm2_Smooth(m,:) = pre_delay_post_spikeRate2_Smooth(m,:)./max(pre_delay_post_spikeRate2_Smooth(m,:));
            end

            spikeRate1_Combined = mean(spikeRate1,1);
            spikeRate2_Combined = mean(spikeRate2,1);
            
            spikeRate1_Combined_Smooth = gaussfilt((1:length(spikeRate1_Combined))*p.timeBin,spikeRate1_Combined,p.gaussSigma);
            spikeRate2_Combined_Smooth = gaussfilt((1:length(spikeRate1_Combined))*p.timeBin,spikeRate2_Combined,p.gaussSigma);   
            
            pre_delay_post_spikeRate1_Combined = mean(pre_delay_post_spikeRate1,1);
            pre_delay_post_spikeRate2_Combined = mean(pre_delay_post_spikeRate2,1);
            
            pre_delay_post_spikeRate1_Combined_Smooth = gaussfilt((1:length(pre_delay_post_spikeRate1_Combined))*p.timeBin,pre_delay_post_spikeRate1_Combined,p.gaussSigma);
            pre_delay_post_spikeRate2_Combined_Smooth = gaussfilt((1:length(pre_delay_post_spikeRate1_Combined))*p.timeBin,pre_delay_post_spikeRate2_Combined,p.gaussSigma);   
            
            
            % fix all NaNs to 0;
            % NaNs means the trial with all 0 Hz bins in this code
            spikeRate1(isnan(spikeRate1)) = 0;
            spikeRate1_Smooth(isnan(spikeRate1_Smooth)) = 0;
            spikeRateNorm1(isnan(spikeRateNorm1)) = 0;
            spikeRateNorm1_Smooth(isnan(spikeRateNorm1_Smooth)) = 0;            
            spikeRate2(isnan(spikeRate2)) = 0;
            spikeRate2_Smooth(isnan(spikeRate2_Smooth)) = 0;
            spikeRateNorm2(isnan(spikeRateNorm2)) = 0;
            spikeRateNorm2_Smooth(isnan(spikeRateNorm2_Smooth)) = 0;            
            spikeRate1_Combined(isnan(spikeRate1_Combined)) = 0;
            spikeRate2_Combined(isnan(spikeRate2_Combined)) = 0;
            spikeRate1_Combined_Smooth(isnan(spikeRate1_Combined_Smooth)) = 0;
            spikeRate2_Combined_Smooth(isnan(spikeRate2_Combined_Smooth)) = 0;
            
            DelayFire_Assembly_onoff.off.(sessDirs{j}).spikeRate1{k} = spikeRate1;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).spikeRate1_Smooth{k} = spikeRate1_Smooth;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).spikeRateNorm1{k} = spikeRateNorm1;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).spikeRateNorm1_Smooth{k} = spikeRateNorm1_Smooth;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).spikeRate2{k} = spikeRate2;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).spikeRate2_Smooth{k} = spikeRate2_Smooth;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).spikeRateNorm2{k} = spikeRateNorm2;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).spikeRateNorm2_Smooth{k} = spikeRateNorm2_Smooth;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).spikeRate1_Combined{k} = spikeRate1_Combined;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).spikeRate2_Combined{k} = spikeRate2_Combined;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).spikeRate1_Combined_Smooth{k} = spikeRate1_Combined_Smooth;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).spikeRate2_Combined_Smooth{k} = spikeRate2_Combined_Smooth;
             
            DelayFire_Assembly_onoff.off.(sessDirs{j}).pre_spikeRate1{k} = pre_spikeRate1;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).pre_spikeRate1_Smooth{k} = pre_spikeRate1_Smooth;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).pre_spikeRateNorm1{k} = pre_spikeRateNorm1;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).pre_spikeRateNorm1_Smooth{k} = pre_spikeRateNorm1_Smooth;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).pre_spikeRate2{k} = pre_spikeRate2;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).pre_spikeRate2_Smooth{k} = pre_spikeRate2_Smooth;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).pre_spikeRateNorm2{k} = pre_spikeRateNorm2;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).pre_spikeRateNorm2_Smooth{k} = pre_spikeRateNorm2_Smooth;

            DelayFire_Assembly_onoff.off.(sessDirs{j}).post_spikeRate1{k} = post_spikeRate1;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).post_spikeRate1_Smooth{k} = post_spikeRate1_Smooth;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).post_spikeRateNorm1{k} = post_spikeRateNorm1;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).post_spikeRateNorm1_Smooth{k} = post_spikeRateNorm1_Smooth;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).post_spikeRate2{k} = post_spikeRate2;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).post_spikeRate2_Smooth{k} = post_spikeRate2_Smooth;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).post_spikeRateNorm2{k} = post_spikeRateNorm2;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).post_spikeRateNorm2_Smooth{k} = post_spikeRateNorm2_Smooth;
            
            DelayFire_Assembly_onoff.off.(sessDirs{j}).pre_delay_post_spikeRate1_Combined{k} = pre_delay_post_spikeRate1_Combined;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).pre_delay_post_spikeRate2_Combined{k} = pre_delay_post_spikeRate2_Combined;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).pre_delay_post_spikeRate1_Combined_Smooth{k} = pre_delay_post_spikeRate1_Combined_Smooth;
            DelayFire_Assembly_onoff.off.(sessDirs{j}).pre_delay_post_spikeRate2_Combined_Smooth{k} = pre_delay_post_spikeRate2_Combined_Smooth;
            
            figure(1)            
            subplot('Position',pos1)
            imagesc([pre_spikeRate1_Smooth,spikeRate1_Smooth,post_spikeRate1_Smooth])
            colormap(jet)
            TITLE1 = sprintf('%s%d%s%d%s%d%s%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-off_Assembly-',k,'-',sessDirs{j},'-eventTime');
            title({TITLE1},'Interpreter','None')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
            ylabel('Trials')
%             xlim([0 delayBinCount])
                       
            subplot('Position',pos2)
            imagesc(pre_delay_post_spikeRate1_Combined_Smooth)
            colormap(jet)     
            title('Combined time map')
            xlabel('Time (Sec)')
            ylabel('Trials')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
                                   
            subplot('Position',pos3)
            imagesc([pre_spikeRate2_Smooth,spikeRate2_Smooth,post_spikeRate2_Smooth])
            colormap(jet)
            title('Event strength','Interpreter','None')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
%             xlabel('Time (Sec)')
            ylabel('Trials')
%             xlim([0 delayBinCount])
                       
            subplot('Position',pos4)
            imagesc(pre_delay_post_spikeRate2_Combined_Smooth)
            colormap(jet)     
            title('Combined time map')
            xlabel('Time (Sec)')
            ylabel('Trials')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
            
            if p.savePlot == 1
                figure(1)
                figName = sprintf('%s%s%d%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-offDelayAssembly-',k,'-',sessDirs{j},'-DelayFiring');
                print(figName,'-dpng','-r300');            
            end          
            close all
        end
        
    end
    
    sessName2 = {'on10','off10','on30','off30'};
    for j = 1:length(sessName2)
        sess1 = strcat(sessName2{j},'_1');
        sess2 = strcat(sessName2{j},'_2');
        if contains(sessName2{j},'10')
            delayBinCount = ceil(10*10^3/p.timeBin);
        else
            delayBinCount = ceil(30*10^3/p.timeBin);
        end
        
        for k = 1:patNum_on
                       
            DelayFire_Assembly_onoff.on.(sessName2{j}).spikeRate1{k} = [DelayFire_Assembly_onoff.on.(sess1).spikeRate1{k};DelayFire_Assembly_onoff.on.(sess2).spikeRate1{k}];
            DelayFire_Assembly_onoff.on.(sessName2{j}).spikeRate1_Smooth{k} = [DelayFire_Assembly_onoff.on.(sess1).spikeRate1_Smooth{k};DelayFire_Assembly_onoff.on.(sess2).spikeRate1_Smooth{k}];
            DelayFire_Assembly_onoff.on.(sessName2{j}).spikeRate1_Combined_Smooth{k} = nanmean(DelayFire_Assembly_onoff.on.(sessName2{j}).spikeRate1_Smooth{k});
            DelayFire_Assembly_onoff.on.(sessName2{j}).spikeRate2{k} = [DelayFire_Assembly_onoff.on.(sess1).spikeRate2{k};DelayFire_Assembly_onoff.on.(sess2).spikeRate2{k}];
            DelayFire_Assembly_onoff.on.(sessName2{j}).spikeRate2_Smooth{k} = [DelayFire_Assembly_onoff.on.(sess1).spikeRate2_Smooth{k};DelayFire_Assembly_onoff.on.(sess2).spikeRate2_Smooth{k}];
            DelayFire_Assembly_onoff.on.(sessName2{j}).spikeRate2_Combined_Smooth{k} = nanmean(DelayFire_Assembly_onoff.on.(sessName2{j}).spikeRate2_Smooth{k});
            
            DelayFire_Assembly_onoff.on.(sessName2{j}).pre_spikeRate1{k} = [DelayFire_Assembly_onoff.on.(sess1).pre_spikeRate1{k};DelayFire_Assembly_onoff.on.(sess2).pre_spikeRate1{k}];
            DelayFire_Assembly_onoff.on.(sessName2{j}).pre_spikeRate1_Smooth{k} = [DelayFire_Assembly_onoff.on.(sess1).pre_spikeRate1_Smooth{k};DelayFire_Assembly_onoff.on.(sess2).pre_spikeRate1_Smooth{k}];
            DelayFire_Assembly_onoff.on.(sessName2{j}).pre_spikeRate1_Combined_Smooth{k} = nanmean(DelayFire_Assembly_onoff.on.(sessName2{j}).pre_spikeRate1_Smooth{k});
            DelayFire_Assembly_onoff.on.(sessName2{j}).pre_spikeRate2{k} = [DelayFire_Assembly_onoff.on.(sess1).pre_spikeRate2{k};DelayFire_Assembly_onoff.on.(sess2).pre_spikeRate2{k}];
            DelayFire_Assembly_onoff.on.(sessName2{j}).pre_spikeRate2_Smooth{k} = [DelayFire_Assembly_onoff.on.(sess1).pre_spikeRate2_Smooth{k};DelayFire_Assembly_onoff.on.(sess2).pre_spikeRate2_Smooth{k}];
            DelayFire_Assembly_onoff.on.(sessName2{j}).pre_spikeRate2_Combined_Smooth{k} = nanmean(DelayFire_Assembly_onoff.on.(sessName2{j}).pre_spikeRate2_Smooth{k});
            
            DelayFire_Assembly_onoff.on.(sessName2{j}).post_spikeRate1{k} = [DelayFire_Assembly_onoff.on.(sess1).post_spikeRate1{k};DelayFire_Assembly_onoff.on.(sess2).post_spikeRate1{k}];
            DelayFire_Assembly_onoff.on.(sessName2{j}).post_spikeRate1_Smooth{k} = [DelayFire_Assembly_onoff.on.(sess1).post_spikeRate1_Smooth{k};DelayFire_Assembly_onoff.on.(sess2).post_spikeRate1_Smooth{k}];
            DelayFire_Assembly_onoff.on.(sessName2{j}).post_spikeRate1_Combined_Smooth{k} = nanmean(DelayFire_Assembly_onoff.on.(sessName2{j}).post_spikeRate1_Smooth{k});
            DelayFire_Assembly_onoff.on.(sessName2{j}).post_spikeRate2{k} = [DelayFire_Assembly_onoff.on.(sess1).post_spikeRate2{k};DelayFire_Assembly_onoff.on.(sess2).post_spikeRate2{k}];
            DelayFire_Assembly_onoff.on.(sessName2{j}).post_spikeRate2_Smooth{k} = [DelayFire_Assembly_onoff.on.(sess1).post_spikeRate2_Smooth{k};DelayFire_Assembly_onoff.on.(sess2).post_spikeRate2_Smooth{k}];
            DelayFire_Assembly_onoff.on.(sessName2{j}).post_spikeRate2_Combined_Smooth{k} = nanmean(DelayFire_Assembly_onoff.on.(sessName2{j}).post_spikeRate1_Smooth{k});       
        
            h1 = figure(1);
            h1.Position = [100 100 900 600];
            % define where each subplot should be
            pos1 = [0.1 0.35 0.3 0.6];
            pos2 = [0.1 0.15 0.3 0.1];
            pos3 = [0.6 0.35 0.3 0.6];
            pos4 = [0.6 0.15 0.3 0.1];
            subplot('Position',pos1)
            imagesc([DelayFire_Assembly_onoff.on.(sessName2{j}).pre_spikeRate1_Smooth{k},DelayFire_Assembly_onoff.on.(sessName2{j}).spikeRate1_Smooth{k},DelayFire_Assembly_onoff.on.(sessName2{j}).post_spikeRate1_Smooth{k}])
            colormap(jet)
            TITLE1 = sprintf('%s%d%s%d%s%d%s%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-on_Assembly-',k,'-',sessName2{j},'-eventTime');
            title({TITLE1},'Interpreter','None')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
            ylabel('Trials')
%             xlim([0 delayBinCount])
                       
            subplot('Position',pos2)
            imagesc([DelayFire_Assembly_onoff.on.(sessName2{j}).pre_spikeRate1_Combined_Smooth{k},DelayFire_Assembly_onoff.on.(sessName2{j}).spikeRate1_Combined_Smooth{k},DelayFire_Assembly_onoff.on.(sessName2{j}).post_spikeRate1_Combined_Smooth{k}])
            colormap(jet)     
            title('Combined time map')
            xlabel('Time (Sec)')
            ylabel('Trials')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
                                   
            subplot('Position',pos3)
            imagesc([DelayFire_Assembly_onoff.on.(sessName2{j}).pre_spikeRate2_Smooth{k},DelayFire_Assembly_onoff.on.(sessName2{j}).spikeRate2_Smooth{k},DelayFire_Assembly_onoff.on.(sessName2{j}).post_spikeRate2_Smooth{k}])
            colormap(jet)
            title('Event strength','Interpreter','None')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
%             xlabel('Time (Sec)')
            ylabel('Trials')
%             xlim([0 delayBinCount])
                       
            subplot('Position',pos4)
            imagesc([DelayFire_Assembly_onoff.on.(sessName2{j}).pre_spikeRate2_Combined_Smooth{k},DelayFire_Assembly_onoff.on.(sessName2{j}).spikeRate2_Combined_Smooth{k},DelayFire_Assembly_onoff.on.(sessName2{j}).post_spikeRate2_Combined_Smooth{k}])
            colormap(jet)     
            title('Combined time map')
            xlabel('Time (Sec)')
            ylabel('Trials')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
            
            if p.savePlot == 1
                figure(1)
                figName = sprintf('%s%s%d%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-onDelayAssembly-',k,'-',sessName2{j},'-DelayFiring');
                print(figName,'-dpng','-r300');            
            end          
            close all
            
        end
        
        for k = 1:patNum_off
                       
            DelayFire_Assembly_onoff.off.(sessName2{j}).spikeRate1{k} = [DelayFire_Assembly_onoff.off.(sess1).spikeRate1{k};DelayFire_Assembly_onoff.off.(sess2).spikeRate1{k}];
            DelayFire_Assembly_onoff.off.(sessName2{j}).spikeRate1_Smooth{k} = [DelayFire_Assembly_onoff.off.(sess1).spikeRate1_Smooth{k};DelayFire_Assembly_onoff.off.(sess2).spikeRate1_Smooth{k}];
            DelayFire_Assembly_onoff.off.(sessName2{j}).spikeRate1_Combined_Smooth{k} = nanmean(DelayFire_Assembly_onoff.off.(sessName2{j}).spikeRate1_Smooth{k});
            DelayFire_Assembly_onoff.off.(sessName2{j}).spikeRate2{k} = [DelayFire_Assembly_onoff.off.(sess1).spikeRate2{k};DelayFire_Assembly_onoff.off.(sess2).spikeRate2{k}];
            DelayFire_Assembly_onoff.off.(sessName2{j}).spikeRate2_Smooth{k} = [DelayFire_Assembly_onoff.off.(sess1).spikeRate2_Smooth{k};DelayFire_Assembly_onoff.off.(sess2).spikeRate2_Smooth{k}];
            DelayFire_Assembly_onoff.off.(sessName2{j}).spikeRate2_Combined_Smooth{k} = nanmean(DelayFire_Assembly_onoff.off.(sessName2{j}).spikeRate2_Smooth{k});
            
            DelayFire_Assembly_onoff.off.(sessName2{j}).pre_spikeRate1{k} = [DelayFire_Assembly_onoff.off.(sess1).pre_spikeRate1{k};DelayFire_Assembly_onoff.off.(sess2).pre_spikeRate1{k}];
            DelayFire_Assembly_onoff.off.(sessName2{j}).pre_spikeRate1_Smooth{k} = [DelayFire_Assembly_onoff.off.(sess1).pre_spikeRate1_Smooth{k};DelayFire_Assembly_onoff.off.(sess2).pre_spikeRate1_Smooth{k}];
            DelayFire_Assembly_onoff.off.(sessName2{j}).pre_spikeRate1_Combined_Smooth{k} = nanmean(DelayFire_Assembly_onoff.off.(sessName2{j}).pre_spikeRate1_Smooth{k});
            DelayFire_Assembly_onoff.off.(sessName2{j}).pre_spikeRate2{k} = [DelayFire_Assembly_onoff.off.(sess1).pre_spikeRate2{k};DelayFire_Assembly_onoff.off.(sess2).pre_spikeRate2{k}];
            DelayFire_Assembly_onoff.off.(sessName2{j}).pre_spikeRate2_Smooth{k} = [DelayFire_Assembly_onoff.off.(sess1).pre_spikeRate2_Smooth{k};DelayFire_Assembly_onoff.off.(sess2).pre_spikeRate2_Smooth{k}];
            DelayFire_Assembly_onoff.off.(sessName2{j}).pre_spikeRate2_Combined_Smooth{k} = nanmean(DelayFire_Assembly_onoff.off.(sessName2{j}).pre_spikeRate2_Smooth{k});
            
            DelayFire_Assembly_onoff.off.(sessName2{j}).post_spikeRate1{k} = [DelayFire_Assembly_onoff.off.(sess1).post_spikeRate1{k};DelayFire_Assembly_onoff.off.(sess2).post_spikeRate1{k}];
            DelayFire_Assembly_onoff.off.(sessName2{j}).post_spikeRate1_Smooth{k} = [DelayFire_Assembly_onoff.off.(sess1).post_spikeRate1_Smooth{k};DelayFire_Assembly_onoff.off.(sess2).post_spikeRate1_Smooth{k}];
            DelayFire_Assembly_onoff.off.(sessName2{j}).post_spikeRate1_Combined_Smooth{k} = nanmean(DelayFire_Assembly_onoff.off.(sessName2{j}).post_spikeRate1_Smooth{k});
            DelayFire_Assembly_onoff.off.(sessName2{j}).post_spikeRate2{k} = [DelayFire_Assembly_onoff.off.(sess1).post_spikeRate2{k};DelayFire_Assembly_onoff.off.(sess2).post_spikeRate2{k}];
            DelayFire_Assembly_onoff.off.(sessName2{j}).post_spikeRate2_Smooth{k} = [DelayFire_Assembly_onoff.off.(sess1).post_spikeRate2_Smooth{k};DelayFire_Assembly_onoff.off.(sess2).post_spikeRate2_Smooth{k}];
            DelayFire_Assembly_onoff.off.(sessName2{j}).post_spikeRate2_Combined_Smooth{k} = nanmean(DelayFire_Assembly_onoff.off.(sessName2{j}).post_spikeRate1_Smooth{k});       
        
            h1 = figure(1);
            h1.Position = [100 100 900 600];
            % define where each subplot should be
            pos1 = [0.1 0.35 0.3 0.6];
            pos2 = [0.1 0.15 0.3 0.1];
            pos3 = [0.6 0.35 0.3 0.6];
            pos4 = [0.6 0.15 0.3 0.1];
            subplot('Position',pos1)
            imagesc([DelayFire_Assembly_onoff.off.(sessName2{j}).pre_spikeRate1_Smooth{k},DelayFire_Assembly_onoff.off.(sessName2{j}).spikeRate1_Smooth{k},DelayFire_Assembly_onoff.off.(sessName2{j}).post_spikeRate1_Smooth{k}])
            colormap(jet)
            TITLE1 = sprintf('%s%d%s%d%s%d%s%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-off_Assembly-',k,'-',sessName2{j},'-eventTime');
            title({TITLE1},'Interpreter','None')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
            ylabel('Trials')
%             xlim([0 delayBinCount])
                       
            subplot('Position',pos2)
            imagesc([DelayFire_Assembly_onoff.off.(sessName2{j}).pre_spikeRate1_Combined_Smooth{k},DelayFire_Assembly_onoff.off.(sessName2{j}).spikeRate1_Combined_Smooth{k},DelayFire_Assembly_onoff.off.(sessName2{j}).post_spikeRate1_Combined_Smooth{k}])
            colormap(jet)     
            title('Combined time map')
            xlabel('Time (Sec)')
            ylabel('Trials')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
                                   
            subplot('Position',pos3)
            imagesc([DelayFire_Assembly_onoff.off.(sessName2{j}).pre_spikeRate2_Smooth{k},DelayFire_Assembly_onoff.off.(sessName2{j}).spikeRate2_Smooth{k},DelayFire_Assembly_onoff.off.(sessName2{j}).post_spikeRate2_Smooth{k}])
            colormap(jet)
            title('Event strength','Interpreter','None')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
%             xlabel('Time (Sec)')
            ylabel('Trials')
%             xlim([0 delayBinCount])
                       
            subplot('Position',pos4)
            imagesc([DelayFire_Assembly_onoff.off.(sessName2{j}).pre_spikeRate2_Combined_Smooth{k},DelayFire_Assembly_onoff.off.(sessName2{j}).spikeRate2_Combined_Smooth{k},DelayFire_Assembly_onoff.off.(sessName2{j}).post_spikeRate2_Combined_Smooth{k}])
            colormap(jet)     
            title('Combined time map')
            xlabel('Time (Sec)')
            ylabel('Trials')
            xTick= [0 preBinCount,preBinCount+delayBinCount,preBinCount+delayBinCount+postBinCount];
            xTickLabel = [-1*p.preDelayTime 0 maxT maxT++p.postDelayTime]/1000;    
            set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
            
            if p.savePlot == 1
                figure(1)
                figName = sprintf('%s%s%d%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-offDelayAssembly-',k,'-',sessName2{j},'-DelayFiring');
                print(figName,'-dpng','-r300');            
            end          
            close all
        
        end
    end
    
    if p.writeToFile == 1
        save(fullfile(savedir2,'DelayFire_Assembly_onoff.mat'), 'DelayFire_Assembly_onoff');
    end
    clear DelayFire_Assembly_onoff
    fprintf('Finished position analysis for session %d\n',i);
end

end