% calculate cell assemblies throught the sleep-Fig8-sleep session
% treat it as a unit and plot rate map
% Li Yuan, UCSD, May-02-2022
% 
function Fig8Treadmill_PCAICA(inFile,AnalyzeSes)

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

p.eventThres = 1; % p.eventThres * sd + mean

opts.Patterns.method = 'ICA';
opts.threshold.method = 'MarcenkoPastur';
opts.Patterns.number_of_iterations = 500;

% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    sessDirs = sessInfo(i).sessDirs;
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
    CellAssembly_WholeSes.rat = sessInfo(i).animal;
    CellAssembly_WholeSes.day = sessInfo(i).day;
    CellAssembly_WholeSes.binWidth = p.timeBinWidth;
    
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    
    % load spikes from each main session
    % get event timestamp
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    
    spkTrain_Fig8 = [];
    rateBin_Fig8_Temp = [];   
    for j = 1:length(sessDirs) 
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
    [rateBin_Fig8,rateBin_Fig8Ind] = sort(rateBin_Fig8_Temp);
    spkTrain_Fig8 = spkTrain_Fig8(rateLabel,rateBin_Fig8Ind);
    
    % sleep sessions      
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
    
    % get assemblies
    %% ICA Assembly
    AssemblyTemplates = assembly_patterns(spkTrain_Fig8,opts);
    % define if it is a cell assembly
    cellID = 1:sum(rateLabel);
    nPatterns = size(AssemblyTemplates,2);
    AssmblPtrnCellIDs = cell(1,nPatterns);
    validInd = zeros(1,nPatterns);
    assmblPtrnWgtsTemp = [];
    for pIdx = 1:nPatterns
        assmblPtrnWgtsTemp(:,pIdx) = AssemblyTemplates(:,pIdx)./norm(AssemblyTemplates(:,pIdx));
        if max(assmblPtrnWgtsTemp(:,pIdx)) < abs(min(assmblPtrnWgtsTemp(:,pIdx)))
            assmblPtrnWgtsTemp(:,pIdx) = -assmblPtrnWgtsTemp(:,pIdx);
        end
%         memberInds = find(assmblPtrnWgtsTemp(:,pIdx) > ...
%             (mean(assmblPtrnWgtsTemp(:,pIdx)) + 1*std(assmblPtrnWgtsTemp(:,pIdx))));
        memberInds = find(assmblPtrnWgtsTemp(:,pIdx) > (1/sqrt(size(assmblPtrnWgtsTemp,1))));
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
    
    CellAssembly_WholeSes.ValidCellNum = sum(rateLabel);
    CellAssembly_WholeSes.ValidCellLabel = rateLabel;
    CellAssembly_WholeSes.ValidCell = TList{rateLabel};
    CellAssembly_WholeSes.patNum = patNum;
    CellAssembly_WholeSes.AssmblPtrnCellIDs = AssmblPtrnCellIDs;
    CellAssembly_WholeSes.AssmblWght = AssmblWght;
    CellAssembly_WholeSes.AssmblStrength = AssmblStrength;
    CellAssembly_WholeSes.event_Time = event_Time;
    CellAssembly_WholeSes.event_strength = event_strength;
    CellAssembly_WholeSes.event_Num = event_Num;
    CellAssembly_WholeSes.binTime = [rateBin_sleep.sleep1,rateBin_Fig8,rateBin_sleep.sleep2];
    
    if p.writeToFile
        fileName = sprintf('%s%d%s','CellAssembly_WholeSes-',p.timeBinWidth*10^3,'ms.mat');
        save(fullfile(savedir2,fileName), 'CellAssembly_WholeSes');
    end
    clear CellAssembly_WholeSes
    
    fprintf('Finished analysis for session %d\n',i)
   
end