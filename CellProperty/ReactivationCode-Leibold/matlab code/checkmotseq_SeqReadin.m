% Read in sequences and form input for checkmotseq.m
% run after code function Fig8TreadmillReactivation_WholeSession(inFile,AnalyzeSes)
% Nov-19-2021, Li Yuan, UCSD
function checkmotseq_SeqReadin(inFile,AnalyzeSes)

p.writeToFile = 1;
p.cellNum = 16;

% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    % load reactivation file
    reactFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_WholeSes.mat');
    load(reactFile);
    % load reactivation delay specific file
    reactDelayFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_Delay.mat');
    load(reactDelayFile);
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);

    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
%     [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
%     clusterNum = length(TList);

    % get valid cell ind
    % rate < 5 Hz in whole session
    % time cell / non-time cell
    clusterNum = length(SpikeProp.max_AvgRate);
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5;
    
    % get each phase names (no delay etc)
    sesDirs = sessInfo(i).sessDirs;
    
    
%     if sum(rateLabel) >= p.cellNum
        DelayPopFire_RegionSeq.timeWindow = DelayPopFire_WholeSes.timeWindow;
        DelayPopFire_RegionSeq.timeWindowIncrement = DelayPopFire_WholeSes.timeWindowIncrement;
        DelayPopFire_RegionSeq.popEventThres = DelayPopFire_WholeSes.popEventThres;
        
        DelayPopFire_RegionSeq.rat = sessInfo(i).animal;
        DelayPopFire_RegionSeq.day = sessInfo(i).day;
        DelayPopFire_RegionSeq.clusterNum = clusterNum;
        DelayPopFire_RegionSeq.PyrCellNum = sum(rateLabel);
    
        % Fig 8 maze & sleep box
        
        for j = 1:length(sesDirs)
      
            DelayPopFire_RegionSeq.(sesDirs{j}).eventNum = DelayPopFire_WholeSes.(sesDirs{j}).eventNum;
            DelayPopFire_RegionSeq.(sesDirs{j}).eventStartTs = DelayPopFire_WholeSes.(sesDirs{j}).eventStartTs;
            DelayPopFire_RegionSeq.(sesDirs{j}).eventEndTs = DelayPopFire_WholeSes.(sesDirs{j}).eventEndTs;
            DelayPopFire_RegionSeq.(sesDirs{j}).locOrder = DelayPopFire_WholeSes.(sesDirs{j}).locOrder;
            DelayPopFire_RegionSeq.(sesDirs{j}).locTime = DelayPopFire_WholeSes.(sesDirs{j}).locTime;
            

            DelayPopFire_RegionSeq.(sesDirs{j}).eventTsp = DelayPopFire_WholeSes.(sesDirs{j}).eventTsp;
            DelayPopFire_RegionSeq.(sesDirs{j}).firstSpkTs = DelayPopFire_WholeSes.(sesDirs{j}).eventTsp2Start;
            DelayPopFire_RegionSeq.(sesDirs{j}).firstSpkTsNorm = DelayPopFire_WholeSes.(sesDirs{j}).eventTsp2StartNorm;
            DelayPopFire_RegionSeq.(sesDirs{j}).eventspikeMean = DelayPopFire_WholeSes.(sesDirs{j}).eventTsAllspks;
            DelayPopFire_RegionSeq.(sesDirs{j}).eventLoc = DelayPopFire_WholeSes.(sesDirs{j}).eventLoc;
            DelayPopFire_RegionSeq.(sesDirs{j}).eventLocName = DelayPopFire_WholeSes.(sesDirs{j}).eventLocName;

            % 1 return, 2 delay, 3 stem, 4 choice, 5 reward
            for kk = 1:length(DelayPopFire_WholeSes.(sesDirs{j}).locOrder)
                validInd = DelayPopFire_WholeSes.(sesDirs{j}).eventLoc==kk;                
                DelayPopFire_RegionSeq.(sesDirs{j}).eventTsp2Start{kk} = DelayPopFire_WholeSes.(sesDirs{j}).eventTsp2Start(validInd,:);
                DelayPopFire_RegionSeq.(sesDirs{j}).eventTsAllspks{kk} = DelayPopFire_WholeSes.(sesDirs{j}).eventTsAllspks(validInd,:);
                DelayPopFire_RegionSeq.(sesDirs{j}).DirectionLabel{kk} = DelayPopFire_WholeSes.(sesDirs{j}).eventDirectionLabel(validInd);
                DelayPopFire_RegionSeq.(sesDirs{j}).eventInd{kk} = validInd;
            end             
        end
        
        sesDirs = sessInfo(i).sleepDirs;
        for j = 1:length(sesDirs)
      
            DelayPopFire_RegionSeq.(sesDirs{j}).eventNum = DelayPopFire_WholeSes.(sesDirs{j}).eventNum;
            DelayPopFire_RegionSeq.(sesDirs{j}).eventStartTs = DelayPopFire_WholeSes.(sesDirs{j}).eventStartTs;
            DelayPopFire_RegionSeq.(sesDirs{j}).eventEndTs = DelayPopFire_WholeSes.(sesDirs{j}).eventEndTs;
            DelayPopFire_RegionSeq.(sesDirs{j}).locOrder = DelayPopFire_WholeSes.(sesDirs{j}).locOrder;
            DelayPopFire_RegionSeq.(sesDirs{j}).locTime = DelayPopFire_WholeSes.(sesDirs{j}).locTime;

            DelayPopFire_RegionSeq.(sesDirs{j}).eventTsp = DelayPopFire_WholeSes.(sesDirs{j}).eventTsp;
            DelayPopFire_RegionSeq.(sesDirs{j}).firstSpkTs = DelayPopFire_WholeSes.(sesDirs{j}).eventTsp2Start;
            DelayPopFire_RegionSeq.(sesDirs{j}).firstSpkTsNorm = DelayPopFire_WholeSes.(sesDirs{j}).eventTsp2StartNorm;
            DelayPopFire_RegionSeq.(sesDirs{j}).eventTsAllspks = DelayPopFire_WholeSes.(sesDirs{j}).eventTsAllspks;
            DelayPopFire_RegionSeq.(sesDirs{j}).eventLoc = DelayPopFire_WholeSes.(sesDirs{j}).eventLoc;
            DelayPopFire_RegionSeq.(sesDirs{j}).eventLocName = DelayPopFire_WholeSes.(sesDirs{j}).eventLocName;
        end 
        
        % format used in checkmotseq.m 
        % here we define seq as eventTsp2Start together and in L/R turn
 
        % calculate pca for each day and region       
        rewardInd.on10_1 = DelayPopFire_RegionSeq.on10_1.eventInd{5};
        rewardInd.on10_2 = DelayPopFire_RegionSeq.on10_2.eventInd{5};
        rewardInd.on30_1 = DelayPopFire_RegionSeq.on30_1.eventInd{5};
        rewardInd.on30_2 = DelayPopFire_RegionSeq.on30_2.eventInd{5};
        rewardInd.off10_1 = DelayPopFire_RegionSeq.off10_1.eventInd{5};
        rewardInd.off10_2 = DelayPopFire_RegionSeq.off10_2.eventInd{5};
        rewardInd.off30_1 = DelayPopFire_RegionSeq.off30_1.eventInd{5};
        rewardInd.off30_2 = DelayPopFire_RegionSeq.off30_2.eventInd{5};
           
        rewardLeft.on10_1 = ~DelayPopFire_RegionSeq.on10_1.DirectionLabel{5};
        rewardLeft.on10_2 = ~DelayPopFire_RegionSeq.on10_2.DirectionLabel{5};
        rewardLeft.on30_1 = ~DelayPopFire_RegionSeq.on30_1.DirectionLabel{5};
        rewardLeft.on30_2 = ~DelayPopFire_RegionSeq.on30_2.DirectionLabel{5};
        rewardLeft.off10_1 = ~DelayPopFire_RegionSeq.off10_1.DirectionLabel{5};
        rewardLeft.off10_2 = ~DelayPopFire_RegionSeq.off10_2.DirectionLabel{5};
        rewardLeft.off30_1 = ~DelayPopFire_RegionSeq.off30_1.DirectionLabel{5};
        rewardLeft.off30_2 = ~DelayPopFire_RegionSeq.off30_2.DirectionLabel{5};
        

        % on
        % direction label
        on_eventTsp2Start.delayLeftInd = ~[DelayPopFire_Delay.on10_1.eventTurnDir;DelayPopFire_Delay.on30_1.eventTurnDir;...
            DelayPopFire_Delay.on10_2.eventTurnDir;DelayPopFire_Delay.on30_2.eventTurnDir];
        
        on_eventTsp2Start.rewardLeftInd = [rewardLeft.on10_1;rewardLeft.on30_1;rewardLeft.on10_2;rewardLeft.on30_2];
        
        
        % event start ts 
        on_eventTsp2Start.delay = [DelayPopFire_Delay.on10_1.eventTsp2Start;DelayPopFire_Delay.on30_1.eventTsp2Start;...
            DelayPopFire_Delay.on10_2.eventTsp2Start;DelayPopFire_Delay.on30_2.eventTsp2Start];
        on_eventTsp2Start.reward = [DelayPopFire_RegionSeq.on10_1.eventTsp2Start{5};DelayPopFire_RegionSeq.on30_1.eventTsp2Start{5};...
            DelayPopFire_RegionSeq.on10_2.eventTsp2Start{5};DelayPopFire_RegionSeq.on30_2.eventTsp2Start{5}];
        

        % event average ts
        on_eventTsAllspks.delay = [DelayPopFire_Delay.on10_1.eventTsAllspks;DelayPopFire_Delay.on30_1.eventTsAllspks;...
            DelayPopFire_Delay.on10_2.eventTsAllspks;DelayPopFire_Delay.on30_2.eventTsAllspks];
        on_eventTsAllspks.reward = [DelayPopFire_RegionSeq.on10_1.eventTsAllspks{5};DelayPopFire_RegionSeq.on30_1.eventTsAllspks{5};...
            DelayPopFire_RegionSeq.on10_2.eventTsAllspks{5};DelayPopFire_RegionSeq.on30_2.eventTsAllspks{5}];


        % off
        % direction label
        off_eventTsp2Start.delayLeftInd = [DelayPopFire_Delay.off10_1.eventTurnDir;DelayPopFire_Delay.off30_1.eventTurnDir;...
            DelayPopFire_Delay.off10_2.eventTurnDir;DelayPopFire_Delay.off30_2.eventTurnDir];
        off_eventTsp2Start.rewardLeftInd = [rewardLeft.off10_1;rewardLeft.off30_1;rewardLeft.off10_2;rewardLeft.off30_2];
        
        % event start ts
        off_eventTsp2Start.delay = [DelayPopFire_Delay.off10_1.eventTsp2Start;DelayPopFire_Delay.off30_1.eventTsp2Start;...
            DelayPopFire_Delay.off10_2.eventTsp2Start;DelayPopFire_Delay.off30_2.eventTsp2Start];
        off_eventTsp2Start.reward = [DelayPopFire_RegionSeq.off10_1.eventTsp2Start{5};DelayPopFire_RegionSeq.off30_1.eventTsp2Start{5};...
            DelayPopFire_RegionSeq.off10_2.eventTsp2Start{5};DelayPopFire_RegionSeq.off30_2.eventTsp2Start{5}];
        
        % event average ts
        off_eventTsAllspks.delay = [DelayPopFire_Delay.off10_1.eventTsAllspks;DelayPopFire_Delay.off30_1.eventTsAllspks;...
            DelayPopFire_Delay.off10_2.eventTsAllspks;DelayPopFire_Delay.off30_2.eventTsAllspks];
        off_eventTsAllspks.reward = [DelayPopFire_RegionSeq.off10_1.eventTsAllspks{5};DelayPopFire_RegionSeq.off30_1.eventTsAllspks{5};...
            DelayPopFire_RegionSeq.off10_2.eventTsAllspks{5};DelayPopFire_RegionSeq.off30_2.eventTsAllspks{5}];
        
        
        sleep_eventTsp2Start = [DelayPopFire_WholeSes.sleep1.eventTsp2Start;DelayPopFire_WholeSes.sleep2.eventTsp2Start];
            
        
        DelayPopFire_RegionSeq.on_eventTsp2Start = on_eventTsp2Start;
        DelayPopFire_RegionSeq.off_eventTsp2Start = off_eventTsp2Start;
        DelayPopFire_RegionSeq.on_eventTsAllspks = on_eventTsAllspks;
        DelayPopFire_RegionSeq.off_eventTsAllspks = off_eventTsAllspks;
        DelayPopFire_RegionSeq.sleep_eventTsp2Start = sleep_eventTsp2Start;
        DelayPopFire_RegionSeq.sleep1_eventTsp2Start = DelayPopFire_WholeSes.sleep1.eventTsp2Start;
        DelayPopFire_RegionSeq.sleep2_eventTsp2Start = DelayPopFire_WholeSes.sleep2.eventTsp2Start;
        
        if p.writeToFile
            save(fullfile(savedir2,'DelayPopFire_RegionSeq.mat'), 'DelayPopFire_RegionSeq');
        end       
%     end
    
    fprintf('Finished analysis for session %d\n',i);
    clear DelayPopFire_RegionSeq on_eventTsp2Start off_eventTsp2Start sleep_eventTsp2Start
end

    

end