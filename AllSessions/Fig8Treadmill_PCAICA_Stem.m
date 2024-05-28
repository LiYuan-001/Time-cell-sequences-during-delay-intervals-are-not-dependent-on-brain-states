% calculate cell assemblies throught the sleep-Fig8-sleep session
% treat it as a unit and plot rate map
% Li Yuan, UCSD, May-12-2023
% 
function Fig8Treadmill_PCAICA_Stem(inFile,AnalyzeSes)

close all
p.savePlot = 0;
p.writeToFile = 1;

% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'

p.savePlot = 1;
p.writeToFile = 1;

p.timeBinWidth = 25./10^3; % unit sec
% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'
p.spikeFs = 30000;

p.eventThres = 2; % p.eventThres * sd + mean

opts.Patterns.method = 'ICA';
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.number_of_iterations = 500;
opts.spikeBinLimit = 20; % this ICA method sometimes put rarely firing cells into an assembly. 
% So I added a restrction that cells have to be in at least 50 bins to remove that chance that 2 cells happen to have firing in the same bin and have high correlation

% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
%     sessDirs = sessInfo(i).sessDirs;
    sleepDirs = sessInfo(i).sleepDirs;
    
    if p.savePlot
        % Li Yuan, UCSD, 15-Nov-2019
        
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\Assembly spatial map');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
%     % load population file
%     reactFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_WholeSes.mat');
%     load(reactFile);
    % load reactivation delay specific file
%     reactDelayFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_Delay.mat');
%     load(reactDelayFile);
%     
    % initiate the data
    CellAssembly_Stem.rat = sessInfo(i).animal;
    CellAssembly_Stem.day = sessInfo(i).day;
    CellAssembly_Stem.binWidth = p.timeBinWidth;
    
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    
    % load spikes from each main session
    % get event timestamp
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    
    %% treadmill on or off delay specific patterns
    spkTrain_Fig8 = [];
    rateBin_Fig8_Temp = []; 
    
    spkTrain_Delay_on = [];
    rateBin_Delay_on_Temp = [];
    sessDirs = {'on10_1','on10_2','on30_1','on30_2'};
    for j = 1:length(sessDirs) 
        % delay area spike and bin for assembly detection
        % load maps for each cluster
        pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'PathZone.mat');
        load(pathZoneFile);  
        delayFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        trialNum = length(PathZone.posEndT.Center);   
        
        for m = 1:trialNum
            startT = Fig8DelayZonePos.delayPos1.endT(m)+1/30;
            endT = PathZone.posEndT.Center(m);
            rateBin1 = startT:p.timeBinWidth:endT-p.timeBinWidth;
            rateBin2 = startT+p.timeBinWidth:p.timeBinWidth:endT;
            binCount = length(rateBin1);
            spkTrainTemp = zeros(clusterNum,binCount);
            
            for k = 1:clusterNum
                if rateLabel(k) == 1
                    tSpTemp = Spike_Session.(sessDirs{j}){k};
                    fireTemp = zeros(1,binCount);
                    for n = 1:binCount
                        fireTemp(n) = sum(tSpTemp>rateBin1(n) & tSpTemp<=rateBin2(n));
                    end
                    %             spkTrain(k,:) = zscore(fireTemp);
                    spkTrainTemp(k,:) = fireTemp;
                end
            end
            spkTrain_Delay_on = [spkTrain_Delay_on,spkTrainTemp];
            rateBin_Delay_on_Temp = [rateBin_Delay_on_Temp,rateBin1];
        end
        
        %% spike and bin in whole fig8 maze for assembly strength fitting
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        startT = pathData.t(1);
        endT = pathData.t(end);
        rateBin1 = startT:p.timeBinWidth:endT-p.timeBinWidth;
        rateBin2 = startT+p.timeBinWidth:p.timeBinWidth:endT;
        binCount = length(rateBin1);
        spkTrainTemp = zeros(clusterNum,binCount);
        for k = 1:clusterNum
            if rateLabel(k) == 1
                tSpTemp = Spike_Session.(sessDirs{j}){k};
                fireTemp = zeros(1,binCount);
                for n = 1:binCount
                    fireTemp(n) = sum(tSpTemp>rateBin1(n) & tSpTemp<=rateBin2(n));
                end
                %             spkTrain(k,:) = zscore(fireTemp);
                spkTrainTemp(k,:) = fireTemp;
            end
        end       
        spkTrain_Fig8 = [spkTrain_Fig8,spkTrainTemp];
        rateBin_Fig8_Temp = [rateBin_Fig8_Temp,rateBin1];
    end
    
    %% delay off
    spkTrain_Delay_off = [];
    rateBin_Delay_off_Temp = [];
    sessDirs = {'off10_1','off10_2','off30_1','off30_2'};
    for j = 1:length(sessDirs)
        % delay area spike and bin for assembly detection
        pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'PathZone.mat');
        load(pathZoneFile);  
        delayFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        trialNum = length(PathZone.posEndT.Center);               
        
        for m = 1:trialNum
            startT = Fig8DelayZonePos.delayPos1.endT(m)+1/30;
            endT = PathZone.posEndT.Center(m);
            rateBin1 = startT:p.timeBinWidth:endT-p.timeBinWidth;
            rateBin2 = startT+p.timeBinWidth:p.timeBinWidth:endT;
            binCount = length(rateBin1);
            spkTrainTemp = zeros(clusterNum,binCount);
            
            for k = 1:clusterNum
                if rateLabel(k) == 1
                    tSpTemp = Spike_Session.(sessDirs{j}){k};
                    fireTemp = zeros(1,binCount);
                    for n = 1:binCount
                        fireTemp(n) = sum(tSpTemp>rateBin1(n) & tSpTemp<=rateBin2(n));
                    end
                    %             spkTrain(k,:) = zscore(fireTemp);
                    spkTrainTemp(k,:) = fireTemp;
                end
            end
            spkTrain_Delay_off = [spkTrain_Delay_off,spkTrainTemp];
            rateBin_Delay_off_Temp = [rateBin_Delay_off_Temp,rateBin1];
        end
        
        % spike and bin in whole fig8 maze for assembly strength fitting
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        startT = pathData.t(1);
        endT = pathData.t(end);
        rateBin1 = startT:p.timeBinWidth:endT-p.timeBinWidth;
        rateBin2 = startT+p.timeBinWidth:p.timeBinWidth:endT;
        binCount = length(rateBin1);
        spkTrainTemp = zeros(clusterNum,binCount);
        for k = 1:clusterNum
            if rateLabel(k) == 1
                tSpTemp = Spike_Session.(sessDirs{j}){k};
                fireTemp = zeros(1,binCount);
                for n = 1:binCount
                    fireTemp(n) = sum(tSpTemp>rateBin1(n) & tSpTemp<=rateBin2(n));
                end
                %             spkTrain(k,:) = zscore(fireTemp);
                spkTrainTemp(k,:) = fireTemp;
            end
        end       
        spkTrain_Fig8 = [spkTrain_Fig8,spkTrainTemp];
        rateBin_Fig8_Temp = [rateBin_Fig8_Temp,rateBin1];
    end
    
    
    % sort the spike train by time
    [rateBin_Delay_on,rateBin_DelayInd_on] = sort(rateBin_Delay_on_Temp);
    spkTrain_Delay_on = spkTrain_Delay_on(rateLabel,rateBin_DelayInd_on);
    
    % sort the spike train by time
    [rateBin_Delay_off,rateBin_DelayInd_off] = sort(rateBin_Delay_off_Temp);
    spkTrain_Delay_off = spkTrain_Delay_off(rateLabel,rateBin_DelayInd_off);
    
    % sort the spike train by time
    [rateBin_Fig8,rateBin_Fig8Ind] = sort(rateBin_Fig8_Temp);
    spkTrain_Fig8 = spkTrain_Fig8(rateLabel,rateBin_Fig8Ind);
    
    
    %% sleep sessions      
    posFile = fullfile(sessInfo(i).mainDir,'processedData','indataS.mat');
    load(posFile);
    for j =  1:length(sleepDirs)       
        pos = indata(j);
        startT = pos.t(1);
        endT = pos.t(end);
        rateBin1 = startT:p.timeBinWidth:endT-p.timeBinWidth;
        rateBin2 = startT+p.timeBinWidth:p.timeBinWidth:endT;
        binCount = length(rateBin1);
        spkTrainTemp = zeros(clusterNum,binCount);
        for k = 1:clusterNum
            if rateLabel(k) == 1
                tSpTemp = Spike_Session.(sleepDirs{j}){k};
                fireTemp = zeros(1,binCount);
                for n = 1:binCount
                    fireTemp(n) = sum(tSpTemp>rateBin1(n) & tSpTemp<=rateBin2(n));
                end
                %             spkTrain(k,:) = zscore(fireTemp);
                spkTrainTemp(k,:) = fireTemp;
            end
        end
        spkTrain_sleep.(sleepDirs{j}) = spkTrainTemp(rateLabel,:);
        rateBin_sleep.(sleepDirs{j}) = rateBin1;
    end
    
    %% get assemblies for treadmill on delay
    % ICA Assembly
    AssemblyTemplates = assembly_patterns(spkTrain_Delay_on,opts);
    % define if it is a cell assembly
    cellID = 1:sum(rateLabel);
    nPatterns = size(AssemblyTemplates,2);
    AssmblPtrnCellIDs = cell(1,nPatterns);
    validInd = zeros(1,nPatterns);
    assmblPtrnWgtsTemp = [];
    invalidCellTemp = (sum(spkTrain_Delay_on>0,2))<opts.spikeBinLimit;
    AssemblyTemplates(invalidCellTemp,:) = 0;
    
    for pIdx = 1:nPatterns
        assmblPtrnWgtsTemp(:,pIdx) = AssemblyTemplates(:,pIdx)./norm(AssemblyTemplates(:,pIdx));
        if max(assmblPtrnWgtsTemp(:,pIdx)) < abs(min(assmblPtrnWgtsTemp(:,pIdx)))
            assmblPtrnWgtsTemp(:,pIdx) = -assmblPtrnWgtsTemp(:,pIdx);
        end
        memberInds = find(assmblPtrnWgtsTemp(:,pIdx) > ...
            (mean(assmblPtrnWgtsTemp(:,pIdx)) + 1*std(assmblPtrnWgtsTemp(:,pIdx))));
        if length(memberInds)>1
            validInd(pIdx) = 1;
        end
        AssmblPtrnCellIDs{pIdx} = cellID(memberInds);
    end
    AssmblPtrnCellIDs = AssmblPtrnCellIDs(validInd==1);
    AssmblWght = assmblPtrnWgtsTemp(:,validInd==1);    
    
    % detect strength for Fig8 and sleep sessions
    AssmblStrength = assembly_activity(AssmblWght,[spkTrain_sleep.sleep1,spkTrain_Fig8,spkTrain_sleep.sleep2]);
        
    % calculate rate maps for assemblies
    % 1. assembly strength >= 2sd+mean counts
    patNum = size(AssmblWght,2);
    cellNum = zeros(patNum,1);
    event_Time = cell(patNum,1);
    event_strength = cell(patNum,1);
    event_Num = zeros(patNum,1);
    eventRate = zeros(patNum,1);
    
    for kk = 1:size(AssmblWght,2)
        % detect cell number, event num and event strength of this
        % pattern
        strengthTemp = AssmblStrength(kk,:);
        cellNum(kk) = length(AssmblPtrnCellIDs{kk});
        
        % thres 1
        thresVal = mean(strengthTemp) + p.eventThres * std(strengthTemp);
        % thres 2
%         thresVal = 5;
        
        strengthLabel = strengthTemp > thresVal;
        strengthLabelDiff = diff([0,strengthLabel]);
        
        event_Num(kk) = sum(strengthLabel==1);
        timeAll = length(strengthTemp) * p.timeBinWidth;
        eventRate(kk) = sum(strengthLabel==1)./timeAll;
        event_strength{kk} = strengthTemp(strengthLabel==1);
        rateBin_All = [rateBin_sleep.sleep1,rateBin_Fig8,rateBin_sleep.sleep2];
        event_Time{kk} = rateBin_All(strengthLabel==1);
    end
    
    CellAssembly_Stem.DelayOn.ValidCellNum = sum(rateLabel);
    CellAssembly_Stem.DelayOn.ValidCellLabel = rateLabel;
    CellAssembly_Stem.DelayOn.ValidCell = TList{rateLabel};
    CellAssembly_Stem.DelayOn.patNum = patNum;
    CellAssembly_Stem.DelayOn.AssmblPtrnCellIDs = AssmblPtrnCellIDs;
    CellAssembly_Stem.DelayOn.AssmblWght = AssmblWght;
    CellAssembly_Stem.DelayOn.AssmblStrength = AssmblStrength;
    CellAssembly_Stem.DelayOn.event_Time = event_Time;
    CellAssembly_Stem.DelayOn.event_strength = event_strength;
    CellAssembly_Stem.DelayOn.event_Num = event_Num;
    CellAssembly_Stem.DelayOn.binTime = [rateBin_sleep.sleep1,rateBin_Fig8,rateBin_sleep.sleep2];
    
    %% get assemblies for treadmill off delay
    % ICA Assembly
    AssemblyTemplates = assembly_patterns(spkTrain_Delay_off,opts);
    % define if it is a cell assembly
    cellID = 1:sum(rateLabel);
    nPatterns = size(AssemblyTemplates,2);
    AssmblPtrnCellIDs = cell(1,nPatterns);
    validInd = zeros(1,nPatterns);
    assmblPtrnWgtsTemp = [];
    invalidCellTemp = (sum(spkTrain_Delay_off>0,2))<opts.spikeBinLimit;
    AssemblyTemplates(invalidCellTemp,:) = 0;
    
    for pIdx = 1:nPatterns
        assmblPtrnWgtsTemp(:,pIdx) = AssemblyTemplates(:,pIdx)./norm(AssemblyTemplates(:,pIdx));
        if max(assmblPtrnWgtsTemp(:,pIdx)) < abs(min(assmblPtrnWgtsTemp(:,pIdx)))
            assmblPtrnWgtsTemp(:,pIdx) = -assmblPtrnWgtsTemp(:,pIdx);
        end
        memberInds = find(assmblPtrnWgtsTemp(:,pIdx) > ...
            (mean(assmblPtrnWgtsTemp(:,pIdx)) + 1*std(assmblPtrnWgtsTemp(:,pIdx))));
        if length(memberInds)>1
            validInd(pIdx) = 1;
        end
        AssmblPtrnCellIDs{pIdx} = cellID(memberInds);
    end
    AssmblPtrnCellIDs = AssmblPtrnCellIDs(validInd==1);
    AssmblWght = assmblPtrnWgtsTemp(:,validInd==1);    
    
    % detect strength for Fig8 and sleep sessions
    AssmblStrength = assembly_activity(AssmblWght,[spkTrain_sleep.sleep1,spkTrain_Fig8,spkTrain_sleep.sleep2]);
        
    % calculate rate maps for assemblies
    % 1. assembly strength >= 2sd+mean counts
    patNum = size(AssmblWght,2);
    cellNum = zeros(patNum,1);
    event_Time = cell(patNum,1);
    event_strength = cell(patNum,1);
    event_Num = zeros(patNum,1);
    eventRate = zeros(patNum,1);
    
    for kk = 1:size(AssmblWght,2)
        % detect cell number, event num and event strength of this
        % pattern
        strengthTemp = AssmblStrength(kk,:);
        cellNum(kk) = length(AssmblPtrnCellIDs{kk});
        
        % thres 1
        thresVal = mean(strengthTemp) + p.eventThres * std(strengthTemp);
        % thres 2
%         thresVal = 5;
        
        strengthLabel = strengthTemp > thresVal;
        strengthLabelDiff = diff([0,strengthLabel]);
        
        event_Num(kk) = sum(strengthLabel==1);
        timeAll = length(strengthTemp) * p.timeBinWidth;
        eventRate(kk) = sum(strengthLabel==1)./timeAll;
        event_strength{kk} = strengthTemp(strengthLabel==1);
        rateBin_All = [rateBin_sleep.sleep1,rateBin_Fig8,rateBin_sleep.sleep2];
        event_Time{kk} = rateBin_All(strengthLabel==1);
    end
    
    CellAssembly_Stem.DelayOff.ValidCellNum = sum(rateLabel);
    CellAssembly_Stem.DelayOff.ValidCellLabel = rateLabel;
    CellAssembly_Stem.DelayOff.ValidCell = TList{rateLabel};
    CellAssembly_Stem.DelayOff.patNum = patNum;
    CellAssembly_Stem.DelayOff.AssmblPtrnCellIDs = AssmblPtrnCellIDs;
    CellAssembly_Stem.DelayOff.AssmblWght = AssmblWght;
    CellAssembly_Stem.DelayOff.AssmblStrength = AssmblStrength;
    CellAssembly_Stem.DelayOff.event_Time = event_Time;
    CellAssembly_Stem.DelayOff.event_strength = event_strength;
    CellAssembly_Stem.DelayOff.event_Num = event_Num;
    CellAssembly_Stem.DelayOff.binTime = [rateBin_sleep.sleep1,rateBin_Fig8,rateBin_sleep.sleep2];
    
    if p.writeToFile
        fileName = sprintf('%s%d%s','CellAssembly_Stem-',p.timeBinWidth*10^3,'ms.mat');
        save(fullfile(savedir2,fileName), 'CellAssembly_Stem');
    end
    clear CellAssembly_Stem
    
    fprintf('Finished analysis for session %d\n',i)
   
end
