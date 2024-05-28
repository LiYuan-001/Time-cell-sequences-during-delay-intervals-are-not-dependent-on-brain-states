% Visualize population vector reactivation
% 
% Oct-10-2021, Li Yuan, UCSD
function Fig8TreadmillDelay_Reactivation(inFile,AnalyzeSes)
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


% plot spikes for visualization
p.plot = 1;

% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'

% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    close all
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures',sessInfo(i).animal,'-day',sessInfo(i).day,'\Cell pair cofire');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end
    
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
    
%     if sum(rateLabel)<= 20
%         p.popEventThres = 1/4;
%     elseif sum(rateLabel)> 20
%         p.popEventThres = 1/5;
%     end
    
    
    for k = 1:clusterNum
        DelayPopFire_Delay.tList{k} = TList{k}(1:end-2);
    end
    tList = DelayPopFire_Delay.tList;
    
    % initiate the data
    DelayPopFire_Delay.rat = sessInfo(i).animal;
    DelayPopFire_Delay.day = sessInfo(i).day;
    DelayPopFire_Delay.timeBin = p.timeWindow;
    DelayPopFire_Delay.method = 'Cofire cell count in timebin';
    
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    

    for j = 1:length(sessDirs)
        
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        % load turning directions and correctness
        tInfoFile = fullfile(mainDir,sessDirs{j}, 'trialInfo.mat');
        tInfo = load(tInfoFile);
        % -----------------------------------------------------------------
        % use the mode of delay is always 10 sec or 30 sec
        % rather than identify real time in delay zone
        % to make plot cleaner
        % I can change to delay time in delay zone later
        % Li Yuan, 19-Aug-2020, UCSD
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
        delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
        delayTend2 = Fig8DelayZonePos.delayPos2.endT;
        
        trialNum = size(delayTstart1,2);
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
            delayTend2_2 = delayTstart2+maxT;
            rateBin1 = 0:p.timeWindowIncrement:maxT-p.timeWindow-p.timeWindowIncrement; 
            rateBin2 = p.timeWindow:p.timeWindowIncrement:maxT-p.timeWindowIncrement;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
            delayTend2_2 = delayTstart2+maxT;
            rateBin1 = 0:p.timeWindowIncrement:maxT-p.timeWindow; 
            rateBin2 = p.timeWindow:p.timeWindowIncrement:maxT;
        else
            error('Delay time is wrong')
        end       
  
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        binCount = length(rateBin1);
        ts_Def1.spikeCount = zeros(clusterNum,trialNum,binCount);
        ts_Def2.spikeCount = zeros(clusterNum,trialNum,binCount);  
        ts_Def1.cellLabel = zeros(clusterNum,trialNum,binCount);
        ts_Def2.cellLabel = zeros(clusterNum,trialNum,binCount); 
        
        if p.plot == 1
            % initiate the figure to plot spikes
            h = figure(j);
            h.Position = [100,100,1200,800];
            title(strcat('Cofire visualiztion ',sessDirs{j}));
        end
        
        for k = 1:clusterNum
            % get each spike time, change unit to msec from sec
            % ts unit: ms

            tSp = Spike_Session.(sessDirs{j}){k};
            if ~isempty(tSp)
                % only consider pupative pyramidal cells
                if rateLabel(k) == 1
                    for m = 1:trialNum
                        
                        % delay definition 1 start from barrier
                        ts_Delay1 = tSp(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                        ts_Delay1Temp= ts_Delay1-delayTstart1(m);
                        % delay def 2: start from entrance
                        ts_Delay2 = tSp(tSp>delayTstart2(m) & tSp<delayTend2_2(m));
                        ts_Delay2Temp = ts_Delay2-delayTstart2(m);
                        
                        % time bin activate label
                        % def 1, delay start at afloor barrier
                        fireTemp = zeros(1,binCount);
                        for n = 1:binCount
                            fireTemp(n) = sum(ts_Delay1Temp>rateBin1(n) & ts_Delay1Temp<rateBin2(n));
                        end
                        % number of spikes in each time bin
                        ts_Def1.spikeCount(k,m,:) = fireTemp;
                        % whether this cell fires in each time bin
                        ts_Def1.cellLabel(k,m,:) = double(fireTemp>0);
                        
                        % def 2, delay start at entrance                        
                        % delay definition 2
                        fireTemp = zeros(1,binCount);
                        for n = 1:binCount
                            fireTemp(n) = sum(ts_Delay2Temp>rateBin1(n) & ts_Delay2Temp<rateBin2(n))>0;
                        end
                        % number of spikes in each time bin
                        ts_Def2.spikeCount(k,m,:) = fireTemp;
                        % whether this cell fires in each time bin
                        ts_Def2.cellLabel(k,m,:) = double(fireTemp>0);
                        
                        if p.plot == 1
                            % plot all spikes for a cluster in delay
                            if ~isempty(ts_Delay1Temp)
                                figure(j)
                                if rateLabel(k) == 1
                                    plot(ts_Delay1Temp+(m-1)*maxT,-1*k,'k>');
                                    hold on
                                elseif rateLabel(k) == 0
                                    plot(ts_Delay1Temp+(m-1)*maxT,-1*k,'b>');
                                    hold on
                                end
                            end
                        end
                    
                    end
                end
            end
        end
        
        eventNum = 0;
        eventStartTs = [];
        eventEndTs = [];
        eventTsp = [];
        eventTsp2Start = [];
        eventTsp2StartNorm = [];
        eventTsAllspks = [];
        eventTurnDir = [];

        for m = 1:trialNum  
            
            eventDirection = tInfo.direction(m);
            if strcmp(eventDirection,'L')
                eventDirectionLabel= 0;
            elseif strcmp(eventDirection,'R')
                eventDirectionLabel = 1;
            else
                error('Turn direction label is wrong')
            end
            
            turnDir = eventDirectionLabel;
            
            % take the clusterNum*binNum cell firing label matrix out
            cellTemp = squeeze(ts_Def1.cellLabel(:,m,:));            
            tsTemp = squeeze(ts_Def1.spikeCount(:,m,:));    
            % get the number of cells firing in each time bin
            cellPopMat = sum(cellTemp,1);
            % get the spike rate from all cells in each time bin
            ratePopMat = sum(tsTemp,1)/p.timeWindow;
            
            if p.plot == 1
                figure(j)
                %             plot(rateBin1+p.timeWindow/2+(m-1)*maxT,2*ratePopMat(m,:)/max(ratePopMat(m,:))-1,'c')
                plot(rateBin1+p.timeWindow/2+(m-1)*maxT,5*cellPopMat/sum(rateLabel)+2,'m')
            end
            
            % step 2: 
            % get the bins which has significant number of cells firing
            sigCofire = cellPopMat >= floor(p.popEventThres*sum(rateLabel));
            
            % define event by consecutive sigCofire labels
            sigTemp = diff([0,sigCofire]);
            eventStartTemp = find(sigTemp == 1);
            eventEndTemp = find(sigTemp == -1);
            % consider boundary conditions
            if sigCofire(end) == 1 && sigCofire(end-1) == 0
                eventStartTemp = eventStartTemp(1:end-1);
            end
            if sigCofire(end) == 1 && sigCofire(end-1) == 1
                eventEndTemp = [eventEndTemp,length(sigTemp)+1];
            end
            eventEndTemp = eventEndTemp-1;
            
            if any((eventEndTemp-eventStartTemp)<0)
                error('Population event detection is wrong')
            end
            
            eventStartTsTemp = rateBin1(eventStartTemp);
            eventEndTsTemp = rateBin1(eventEndTemp)+p.timeWindow;
            
            eventStartTsTemp2 = eventStartTsTemp+delayTstart1(m);
            eventEndTsTemp2 = eventEndTsTemp+delayTstart1(m);
            
            eventTspTemp = cell(length(eventStartTemp),clusterNum);
            eventTsp2StartTemp = nan(length(eventStartTemp),clusterNum);
            eventTsp2StartNormTemp = nan(length(eventStartTemp),clusterNum);
            eventTsAllspksTemp = nan(length(eventStartTemp),clusterNum);
            eventTurnTemp = ones(length(eventStartTemp),1).*turnDir;
            
            if p.plot == 1
                startPlotTemp = eventStartTsTemp2;
                endPlotTemp = eventEndTsTemp2;
                if ~isempty(startPlotTemp)
                    for kk = 1:length(startPlotTemp)
                        startTs = startPlotTemp(kk)-delayTstart1(m)+(m-1)*maxT;
                        endTs = endPlotTemp(kk)-delayTstart1(m)+(m-1)*maxT;
                        y1 = -1*clusterNum-2;
                        y2 = 5;
                        h = patch([startTs endTs endTs startTs],[y1 y1 y2 y2],'g','FaceAlpha',0.1,'LineStyle','none');
                    end
                end
            end
            
            % step 3
            % remove cells which have certain amount of spikes preceding
            % the event and reaccess the event
            for kk = 1:length(eventStartTemp)         
                % get spikes within this event
                for k = 1:clusterNum
                    
                    tSp = Spike_Session.(sessDirs{j}){k};
                    eventStartTemp3 = eventStartTsTemp2(kk);
                    eventEndTemp3 = eventEndTsTemp2(kk);
                    
                    if rateLabel(k) == 1
                        % get the spikes belongs to the cell this event
                        spkTsTemp = tSp(tSp>=eventStartTemp3 & tSp<eventEndTemp3);
                        % before event starts, apply a spike number
                        % threshold to remove random active cells
                        if ~isempty(spkTsTemp)
%                             if sum(tSp>=spkTsTemp(1)-p.silentWindow & tSp<eventStartTemp3)<=p.spkActiveThres
                            if sum(tSp>=spkTsTemp(1)-p.silentWindow & tSp<eventStartTemp3)<=p.spkActiveThres
                                eventTspTemp{kk,k} = spkTsTemp-eventStartTemp3;
                            end
                        end                  
                        
                        % first spike timing to event start
                        if ~isempty(eventTspTemp{kk,k})
                            % first spike distance from start
                            eventTsp2StartTemp(kk,k) = eventTspTemp{kk,k}(1);
                            % normalize event length as 1
                            % first spike distance from start
                            eventTsp2StartNormTemp(kk,k) = eventTspTemp{kk,k}(1)/(eventEndTemp3-eventStartTemp3);
                            % all spikes average time from start
                            eventTsAllspksTemp(kk,k) = mean(eventTspTemp{kk,k});
                        end                        
                    end                    
                end                          
            end  
            
            % reaccess whether this event has enough cells after removing
            % some cells
            validInd = sum(~isnan(eventTsp2StartTemp),2)>= floor(p.popEventThres*sum(rateLabel));
            
            if p.plot == 1
                startPlotTemp = eventStartTsTemp2(validInd);
                endPlotTemp = eventEndTsTemp2(validInd);
                if ~isempty(startPlotTemp)
                    for kk = 1:length(startPlotTemp)
                        startTs = startPlotTemp(kk)-delayTstart1(m)+(m-1)*maxT;
                        endTs = endPlotTemp(kk)-delayTstart1(m)+(m-1)*maxT;
                        y1 = -1*clusterNum-2;
                        y2 = 5;
                        h = patch([startTs endTs endTs startTs],[y1 y1 y2 y2],'y','FaceAlpha',0.4,'LineStyle','none');
                    end
                end
            end
            
            % keep the valid event only
            eventNum = eventNum + sum(validInd);
            eventStartTs = [eventStartTs,eventStartTsTemp2(validInd)];
            eventEndTs = [eventEndTs,eventEndTsTemp2(validInd)];
            eventTsp = [eventTsp;eventTspTemp(validInd,:)];
            eventTsp2Start = [eventTsp2Start;eventTsp2StartTemp(validInd,:)];
            eventTsp2StartNorm = [eventTsp2StartNorm;eventTsp2StartNormTemp(validInd,:)];
            eventTsAllspks = [eventTsAllspks;eventTsAllspksTemp(validInd,:)];
            eventTurnDir = [eventTurnDir;eventTurnTemp(validInd,:)];
            
            
        end
        
        DelayPopFire_Delay.(sessDirs{j}).eventNum = eventNum;
        DelayPopFire_Delay.(sessDirs{j}).eventStartTs = eventStartTs;
        DelayPopFire_Delay.(sessDirs{j}).eventEndTs = eventEndTs;
        DelayPopFire_Delay.(sessDirs{j}).eventTsp = eventTsp;
        DelayPopFire_Delay.(sessDirs{j}).eventTsp2Start = eventTsp2Start;
        DelayPopFire_Delay.(sessDirs{j}).eventTsp2StartNorm = eventTsp2StartNorm;
        DelayPopFire_Delay.(sessDirs{j}).eventTsAllspks = eventTsAllspks;
        DelayPopFire_Delay.(sessDirs{j}).eventTurnDir = eventTurnDir;    
    end
    
    if p.writeToFile == 1
        save(fullfile(savedir2,'DelayPopFire_Delay.mat'), 'DelayPopFire_Delay');
    end
    
    clear DelayPopFire_Delay
    fprintf('Finished analysis for session %d\n',i);
end

end