% Population vector reactivation
% 
% Oct-29-2021, Li Yuan, UCSD
function Fig8TreadmillReactivation_WholeSession(inFile,AnalyzeSes)
close all

p.savePlot = 0;
p.writeToFile = 1;

% set cofire time window, unit: sec
p.timeWindow = 150./10^3; % time window for calculating population event
% by changing silentWindow to 0 means no silent period need to exist for a
% cell to be assigned in this population event
p.silentWindow = 100./10^3; % silent time window before a event start for a single cell 
% p.silentWindow = 100./10^3; % silent time window before a event start for a single cell
p.spkActiveThres = 2; % spike number allowed in a slient time window for a cell to be included in that pop event
p.timeWindowIncrement = 20./10^3; % time window increment
p.popEventThres = 1/5; % threshold for a significant population event


% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'

% Read in input information
sessInfo = SessInfoImport(inFile);

locOrder = {'Return','Delay','Stem','Choice','Reward'};

for i = AnalyzeSes(1:end)
    
%     if p.savePlot
%         % directory for plot figures
%         % generate a folder for each rat eah day under the current folder
%         savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures',sessInfo(i).animal,'-day',sessInfo(i).day,'\Cell pair cofire');
%         if ~exist(savedir, 'dir')
%             mkdir(savedir);
%         end
%     end
    
    % load time field file
    timeFieldFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayTimeField_Trial_Shuffle.mat');
    load(timeFieldFile);
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);

    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);

    % get valid cell ind
    % time cell / non-time cell
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;  
    
    DelayPopFire_WholeSes.rat = sessInfo(i).animal;
    DelayPopFire_WholeSes.day = sessInfo(i).day;
    DelayPopFire_WholeSes.PyrCellNum = sum(rateLabel);
    DelayPopFire_WholeSes.timeWindow = p.timeWindow;
    DelayPopFire_WholeSes.timeWindowIncrement = p.timeWindowIncrement;
    DelayPopFire_WholeSes.popEventThres = p.popEventThres;
    DelayPopFire_WholeSes.method = 'Cofire cell count in timebin';
    
    for k = 1:clusterNum
        DelayPopFire_WholeSes.tList{k} = TList{k}(1:end-2);
    end
    tList = DelayPopFire_WholeSes.tList;
    
    if length(unique(TList))~=clusterNum
        error('TTList file has repeated clusters')
    end    
    
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sesDirs = sessInfo(i).sessDirs;
    
    % Fig 8 maze
    for j = 1:length(sesDirs)       
        % load analyzed positions
        delayFile = fullfile(mainDir,sesDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        % load path zone time
        pathZoneFile = fullfile(mainDir,sesDirs{j}, 'PathZone.mat');
        load(pathZoneFile);
        % load session x,y,t
        pathDataFile = fullfile(mainDir,sesDirs{j},'pathData.mat');
        pathData = load(pathDataFile);
        % load turning directions and correctness
        tInfoFile = fullfile(mainDir,sesDirs{j}, 'trialInfo.mat');
        tInfo = load(tInfoFile);
             
        
        rateBin1 = pathData.t(1):p.timeWindowIncrement:pathData.t(end)-p.timeWindow-p.timeWindowIncrement;
        rateBin2 = pathData.t(1)+p.timeWindow:p.timeWindowIncrement:pathData.t(end)-p.timeWindowIncrement;
            
        if any(round((rateBin2-rateBin1),2)~=p.timeWindow)
            error('Error in generating time bins')
        end
        % -----------------------------------------------------------------
        trialNum = length(PathZone.posStartT.Return);
        % region start and end matrix
        % [region, trial start, trial end]
        % 1 return, 2 delay, 3 stem, 4 choice, 5 reward
        regTimeBorder = nan(5,trialNum,2);
        regTimeBorder(1,:,:) = [PathZone.posStartT.Return,PathZone.posEndT.Base]; % return+base
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        regTimeBorder(2,:,:) = [Fig8DelayZonePos.delayPos1.startT',Fig8DelayZonePos.delayPos1.endT']; % delay
        regTimeBorder(3,:,:) = [Fig8DelayZonePos.delayPos1.endT'+0.02,PathZone.posEndT.Center]; % stem
        regTimeBorder(4,:,:) = [PathZone.posStartT.Choice,PathZone.posEndT.Choice]; % choice
        regTimeBorder(5,:,:) = [PathZone.posStartT.Reward,PathZone.posEndT.Reward]; % reward
        regionTimeSum = sum(regTimeBorder(:,:,2)-regTimeBorder(:,:,1),2);       

        binCount = length(rateBin1);
        ts.spikeCount = zeros(clusterNum,binCount);
        ts.cellLabel = zeros(clusterNum,binCount);
        ts.spikeTime = cell(clusterNum,1);
        
        for k = 1:clusterNum
            % get each spike time, change unit to msec from sec
            % ts unit: ms
            tSp = Spike_Session.(sesDirs{j}){k};
            if rateLabel(k) == 1
                if ~isempty(tSp)
                    
                    % time bin activate label
                    % def 1, delay start at afloor barrier
                    fireTemp = zeros(1,binCount);
                    for n = 1:binCount
                        fireTemp(n) = sum(tSp>=rateBin1(n) & tSp<rateBin2(n));
                    end
                    ts.spikeCount(k,:) = fireTemp;
                    ts.cellLabel(k,:) = double(fireTemp>0);
                    ts.spikeTime{k} = tSp;
                    
                else
                    ts.spikeTime{k} = double.empty(0,1);
                end
            else
                ts.spikeTime{k} = double.empty(0,1);
            end
%             if sum(ts.spikeCount(k,:))~=10*length(tSp)
%                 error('Error in assign spikes to rate bins')
%             end
        end
        
%         ratePopMat_Smooth = zeros(trialNum(j),binCount);
%         cellPopMat_Smooth = zeros(trialNum(j),binCount);
       
        cellTemp = ts.cellLabel;
        tsTemp = ts.spikeCount;
        cellPopMat = sum(cellTemp,1);
        ratePopMat = sum(tsTemp,1)/p.timeWindow;
        
        sigCofire = cellPopMat >= floor(p.popEventThres*sum(rateLabel));
            
        % define event by consecutive sigCofire labels
        sigTemp = diff([0,sigCofire]);
        eventStart = find(sigTemp == 1);
        eventEnd = find(sigTemp == -1);
        % consider boundary conditions
        if sigCofire(end) == 1 && sigCofire(end-1) == 0
            eventStart = eventStart(1:end-1);
        end
        if sigCofire(end) == 1 && sigCofire(end-1) == 1
            eventEnd = [eventEnd,length(sigTemp)+1];
        end
        eventEnd = eventEnd-1;

        
        eventStartTs = rateBin1(eventStart);
        eventEndTs = rateBin1(eventEnd)+p.timeWindow;
        
        % remove time before first return
        indValid = eventStartTs>=regTimeBorder(1,1,1);
        eventStartTs = eventStartTs(indValid);
        eventEndTs = eventEndTs(indValid);
        
        if any((eventEndTs-eventStartTs)<0) || any((eventEndTs-eventStartTs)>=0.5) 
            error('TIME DETECTION IS WRONG')
        end
        eventNum = length(eventStartTs);       
        
        eventTsp = cell(length(eventStartTs),clusterNum);
        eventTsp2Start = nan(length(eventStartTs),clusterNum);
        eventTsp2StartNorm = nan(length(eventStartTs),clusterNum);
        eventTsAllspks = nan(length(eventStartTs),clusterNum);
        
        eventTrial = nan(length(eventStartTs),1);
        eventLoc = nan(length(eventStartTs),1);
        eventLocName = cell(length(eventStartTs),1);
        eventDirection = cell(length(eventStartTs),1);
        eventDirectionLabel = nan(length(eventStartTs),1);

        for kk = 1:eventNum
            startTs = eventStartTs(kk);
            endTs = eventEndTs(kk);
            
            % assign location tag to the event
            locTemp = startTs>=regTimeBorder(:,:,1)&startTs<=regTimeBorder(:,:,2);
            trialNum = find(sum(locTemp,1)==1);
            locTemp = sum(locTemp,2);
            % double check if event is only assigned to one location
            if sum(locTemp) ~= 1
                % easy fix
                % use next small difference to assign positions
                % 1 return, 2 delay, 3 stem, 4 choice, 5 reward
                posDiff1 = regTimeBorder(:,:,1)-startTs;
                posDiff2 = regTimeBorder(:,:,2)-startTs;
                posInd = abs(posDiff1) < abs(posDiff2);
                posDiff = zeros(size(posDiff1,1),size(posDiff1,2));
                posDiff(posInd) = posDiff1(posInd);
                posDiff(~posInd) = posDiff2(~posInd);
                [minDiff,minDiffInd] = min(min(abs(posDiff)));
                trialNum = minDiffInd;
                if minDiff < 15
                    [~,minDiffInd2] = min(abs(posDiff(:,minDiffInd)));
                    locTemp(minDiffInd2) = 1;
                end
                if sum(locTemp) ~= 1
                    error('Error in event position detection')
                end                
            end
            eventTrial(kk) = trialNum;
            eventLoc(kk) = find(locTemp==1);
            eventLocName(kk) = locOrder(eventLoc(kk));    
            eventDirection(kk) = tInfo.direction(trialNum);
            if strcmp(eventDirection(kk),'L')
                eventDirectionLabel(kk) = 0;
            elseif strcmp(eventDirection(kk),'R')
                eventDirectionLabel(kk) = 1;
            else
                error('Turn direction label is wrong')
            end

                            
            % get spikes within this event
            for k = 1:clusterNum
                tSp = Spike_Session.(sesDirs{j}){k};
                eventStartTemp = eventStartTs(kk);
                eventEndTemp = eventEndTs(kk);

                
                if rateLabel(k) == 1
                    spkTsTemp = tSp(tSp>=eventStartTemp & tSp<eventEndTemp);
                    % before event starts, apply a spike number
                    % threshold to remove random active cells
                    if ~isempty(spkTsTemp)
                        if sum(tSp>=spkTsTemp(1)-p.silentWindow & tSp<eventStartTemp)<=p.spkActiveThres
                            eventTsp{kk,k} = spkTsTemp-eventStartTemp;
                        end
                    end
                        
                    % first spike timing to event start
                    if ~isempty(eventTsp{kk,k})
                        % first spike distance from start
                        eventTsp2Start(kk,k) = eventTsp{kk,k}(1);
                        % normalize event length as 1
                        % first spike distance from start
                        eventTsp2StartNorm(kk,k) = eventTsp{kk,k}(1)/(endTs-startTs);
                        % all spikes average time from start
                        eventTsAllspks(kk,k) = mean(eventTsp{kk,k});
                    end
                end
            end
        end
        validInd = sum(~isnan(eventTsp2Start),2)>= floor(p.popEventThres*sum(rateLabel));
        
        DelayPopFire_WholeSes.(sesDirs{j}).eventNum = sum(validInd);
        DelayPopFire_WholeSes.(sesDirs{j}).eventStartTs = eventStartTs(validInd);
        DelayPopFire_WholeSes.(sesDirs{j}).eventEndTs = eventEndTs(validInd);
        DelayPopFire_WholeSes.(sesDirs{j}).locOrder = locOrder;
        DelayPopFire_WholeSes.(sesDirs{j}).locTime = regionTimeSum;
        
        DelayPopFire_WholeSes.(sesDirs{j}).eventTrial = eventTrial(validInd);
        DelayPopFire_WholeSes.(sesDirs{j}).eventDirectionLabel = eventDirectionLabel(validInd);
        DelayPopFire_WholeSes.(sesDirs{j}).eventDirection = eventDirection(validInd);
        DelayPopFire_WholeSes.(sesDirs{j}).eventTsp = eventTsp(validInd,:);
        DelayPopFire_WholeSes.(sesDirs{j}).eventTsp2Start = eventTsp2Start(validInd,:);
        DelayPopFire_WholeSes.(sesDirs{j}).eventTsp2StartNorm = eventTsp2StartNorm(validInd,:);
        DelayPopFire_WholeSes.(sesDirs{j}).eventTsAllspks = eventTsAllspks(validInd,:);
        DelayPopFire_WholeSes.(sesDirs{j}).eventLoc = eventLoc(validInd);
        DelayPopFire_WholeSes.(sesDirs{j}).eventLocName = eventLocName(validInd);       
    end

    % sleep box
    sesDirs = sessInfo(i).sleepDirs;
    if prod(~isempty(sesDirs))
        for j = 1:length(sesDirs)
            % load analyzed positions
            sleepFile = fullfile(strcat(sessInfo(i).mainDir,'\processedData'),'indataS');
            sleepPos = load(sleepFile);
            sleepTs = sleepPos.indata(j).t;
            
            rateBin1 = sleepTs(1):p.timeWindowIncrement:sleepTs(end)-p.timeWindow-p.timeWindowIncrement;
            rateBin2 = sleepTs(1)+p.timeWindow:p.timeWindowIncrement:sleepTs(end)-p.timeWindowIncrement;
            
            if any(round((rateBin2-rateBin1),2)~=p.timeWindow)
                error('Error in generating time bins')
            end
            % -----------------------------------------------------------------
            regionTimeSum = sleepTs(end)-sleepTs(1);
            
            binCount = length(rateBin1);
            ts.spikeCount = zeros(clusterNum,binCount);
            ts.cellLabel = zeros(clusterNum,binCount);
            ts.spikeTime = cell(clusterNum,1);
            
            for k = 1:clusterNum
                % get each spike time, change unit to msec from sec
                % ts unit: ms
                tSp = Spike_Session.(sesDirs{j}){k};
                if rateLabel(k) == 1
                    if ~isempty(tSp)
                        
                        % time bin activate label
                        % def 1, delay start at afloor barrier
                        fireTemp = zeros(1,binCount);
                        for n = 1:binCount
                            fireTemp(n) = sum(tSp>rateBin1(n) & tSp<rateBin2(n));
                        end
                        ts.spikeCount(k,:) = fireTemp;
                        ts.cellLabel(k,:) = double(fireTemp>0);
                        ts.spikeTime{k} = tSp;
                        
                    else
                        ts.spikeTime{k} = double.empty(0,1);
                    end
                else
                    ts.spikeTime{k} = double.empty(0,1);
                end
            end
            
            %         ratePopMat_Smooth = zeros(trialNum(j),binCount);
            %         cellPopMat_Smooth = zeros(trialNum(j),binCount);
            
            cellTemp = ts.cellLabel;
            tsTemp = ts.spikeCount;
            cellPopMat = sum(cellTemp,1);
            ratePopMat = sum(tsTemp,1)/p.timeWindow;
            sigCofire = cellPopMat >= floor(p.popEventThres*sum(rateLabel));
            
            % define event by consecutive sigCofire labels
            sigTemp = diff([0,sigCofire]);
            eventStart = find(sigTemp == 1);
            eventEnd = find(sigTemp == -1);
            % consider boundary conditions
            if sigCofire(end) == 1 && sigCofire(end-1) == 0
                eventStart = eventStart(1:end-1);
            end
            if sigCofire(end) == 1 && sigCofire(end-1) == 1
                eventEnd = [eventEnd,length(sigTemp)+1];
            end
            eventEnd = eventEnd-1;
            
            if any((eventEnd-eventStart)<0)
                error('Population event detection is wrong')
            end
            
            eventStartTs = rateBin1(eventStart);
            eventEndTs = rateBin1(eventEnd)+p.timeWindow;
            eventNum = length(eventStart);
            
            
            eventTsp = cell(length(eventStart),clusterNum);
            eventTsp2Start = nan(length(eventStart),clusterNum);
            eventTsp2StartNorm = nan(length(eventStart),clusterNum);
            eventTsAllspks = nan(length(eventStart),clusterNum);
            eventLoc = zeros(length(eventStart),1);
            eventLocName = cell(length(eventStart),1);
        
            for kk = 1:eventNum
                startTs = eventStartTs(kk);
                endTs = eventEndTs(kk);
                
                % assign location tag to the event
                eventLocName{kk} = 'sleep';
                
                % get spikes within this event
                for k = 1:clusterNum
                    tSp = Spike_Session.(sesDirs{j}){k};
                    eventStartTemp = eventStartTs(kk);
                    eventEndTemp = eventEndTs(kk);
                    
                    if rateLabel(k) == 1
                        spkTsTemp = tSp(tSp>=eventStartTemp & tSp<eventEndTemp);
                        
                        if ~isempty(spkTsTemp)
                            if sum(tSp>=spkTsTemp(1)-p.silentWindow & tSp<eventStartTemp)<=p.spkActiveThres
                                eventTsp{kk,k} = spkTsTemp-eventStartTemp;
                            end
                        end
                        
                        % first spike timing to event start
                        if ~isempty(eventTsp{kk,k})
                            % first spike distance from start
                            eventTsp2Start(kk,k) = eventTsp{kk,k}(1);
                            % normalize event length as 1
                            % first spike distance from start
                            eventTsp2StartNorm(kk,k) = eventTsp{kk,k}(1)/(endTs-startTs);
                            % all spikes average time from start
                            eventTsAllspks(kk,k) = mean(eventTsp{kk,k});
                            
                        end                        
                    end
                end
            end
            validInd = sum(~isnan(eventTsp2Start),2)>= floor(p.popEventThres*sum(rateLabel));
             
            DelayPopFire_WholeSes.(sesDirs{j}).eventNum = sum(validInd);
            DelayPopFire_WholeSes.(sesDirs{j}).eventStartTs = eventStartTs(validInd);
            DelayPopFire_WholeSes.(sesDirs{j}).eventEndTs = eventEndTs(validInd);
            DelayPopFire_WholeSes.(sesDirs{j}).locOrder = 'sleep';
            DelayPopFire_WholeSes.(sesDirs{j}).locTime = regionTimeSum;

            DelayPopFire_WholeSes.(sesDirs{j}).eventTsp = eventTsp(validInd,:);
            DelayPopFire_WholeSes.(sesDirs{j}).eventTsp2Start = eventTsp2Start(validInd,:);
            DelayPopFire_WholeSes.(sesDirs{j}).eventTsp2StartNorm = eventTsp2StartNorm(validInd,:);
            DelayPopFire_WholeSes.(sesDirs{j}).eventTsAllspks = eventTsAllspks(validInd,:);
            DelayPopFire_WholeSes.(sesDirs{j}).eventLoc = eventLoc(validInd);
            DelayPopFire_WholeSes.(sesDirs{j}).eventLocName = eventLocName(validInd);
        end
    end
    
    if p.writeToFile == 1
        save(fullfile(savedir2,'DelayPopFire_WholeSes.mat'), 'DelayPopFire_WholeSes');
    end
    
    clear DelayPopFire_WholeSes
    fprintf('Finished analysis for session %d\n',i);
end

end