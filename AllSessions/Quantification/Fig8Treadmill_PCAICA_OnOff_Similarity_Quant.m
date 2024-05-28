
function Fig8Treadmill_PCAICA_OnOff_Similarity_Quant(inFile,AnalyzeSes)
close all
% % Read in input information
sessInfo = SessInfoImport(inFile);

AssemblySimilarityOn = [];
AssemblySimilarityOff = [];

AssemblySimilarityOn_Alloff = [];
AssemblySimilarityOff_Allon = [];

AssemblyCellNum_On = [];
AssemblyCellNum_Off = [];

AssemblyDensity_On = [];
AssemblyDensity_Off = [];

% get each phase names (no delay etc)
sessDirs = {'on10','off10','on30','off30'};
for n = 1:length(sessDirs)
    timeMap_oncells.(sessDirs{n}) = [];
    timeMap_offcells.(sessDirs{n}) = [];
end

for i = AnalyzeSes(1:end)
    
    
    Delay_OnOff_File = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_DelayOnOff_LR-25ms.mat');
    load(Delay_OnOff_File);
    % load delayFire map
    delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap_2Session.mat');
    load(delayFile);
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;   
    validCell = find(rateLabel==1);
    
    onWeight = CellAssembly_DelayLR.DelayOn.AssmblWght;
    offWeight = CellAssembly_DelayLR.DelayOff.AssmblWght;
    cellNum = size(onWeight,1);
    
    AssemblyDensity_On = [AssemblyDensity_On,cellNum/CellAssembly_DelayLR.DelayOn.patNum];
    AssemblyDensity_Off = [AssemblyDensity_Off,cellNum/CellAssembly_DelayLR.DelayOff.patNum];
    
    onCellInd = [];
    offCellInd = [];
    for m = 1:CellAssembly_DelayLR.DelayOn.patNum
        AssemblyCellNum_On = [AssemblyCellNum_On,length(CellAssembly_DelayLR.DelayOn.AssmblPtrnCellIDs{m})];
        for k = 1:length(CellAssembly_DelayLR.DelayOn.AssmblPtrnCellIDs{m})
            cellInd = validCell(CellAssembly_DelayLR.DelayOn.AssmblPtrnCellIDs{m}(k));
            onCellInd = [onCellInd,cellInd];
        end
    end
    onCellInd = unique(onCellInd);
    for n = 1:length(sessDirs)
        for k = 1:length(onCellInd)
            timeMap_oncells.(sessDirs{n}) = [timeMap_oncells.(sessDirs{n});Fig8DelayTimeMap_2Session.(sessDirs{n}).spikeRate1_Combined_Smooth{onCellInd(k)}];
        end
    end
            
    for m = 1:CellAssembly_DelayLR.DelayOff.patNum
        AssemblyCellNum_Off = [AssemblyCellNum_Off,length(CellAssembly_DelayLR.DelayOff.AssmblPtrnCellIDs{m})];
        for k = 1:length(CellAssembly_DelayLR.DelayOff.AssmblPtrnCellIDs{m})
            cellInd = validCell(CellAssembly_DelayLR.DelayOff.AssmblPtrnCellIDs{m}(k));
            offCellInd = [offCellInd,cellInd];
        end
    end
    offCellInd = unique(offCellInd);
    for n = 1:length(sessDirs)
        for k = 1:length(offCellInd)
            timeMap_offcells.(sessDirs{n}) = [timeMap_offcells.(sessDirs{n});Fig8DelayTimeMap_2Session.(sessDirs{n}).spikeRate1_Combined_Smooth{offCellInd(k)}];
        end
    end
    
    for m = 1:CellAssembly_DelayLR.DelayOn.patNum
        sim_Temp = [];
        pattern_On = CellAssembly_DelayLR.DelayOn.AssmblWght(:,m);
        for j = 1:CellAssembly_DelayLR.DelayOff.patNum
            pattern_Off = CellAssembly_DelayLR.DelayOff.AssmblWght(:,j);
            sim_Temp(j) = dot(pattern_On,pattern_Off)/(norm(pattern_On)*norm(pattern_Off));  
            
%             figure
%             subplot(1,2,1)
%             stem(onWeight(:,m),'k');
%             hold on
%             stem(CellAssembly_DelayLR.DelayOn.AssmblPtrnCellIDs{m},...
%                 onWeight(CellAssembly_DelayLR.DelayOn.AssmblPtrnCellIDs{m},m),'r')
%             xlim([0.5 cellNum+0.5])
%             view([90 -90])
%             title('Delay on assembly')
%             
%             subplot(1,2,2)
%             stem(offWeight(:,j),'k');
%             hold on
%             stem(CellAssembly_DelayLR.DelayOff.AssmblPtrnCellIDs{j},...
%                 offWeight(CellAssembly_DelayLR.DelayOff.AssmblPtrnCellIDs{j},j),'r')
%             xlim([0.5 cellNum+0.5])
%             view([90 -90])
%             TITLE1 = 'Delay off assembly';
%             TITLE2 = sprintf('%s%1.2f','Cos similarity = ',sim_Temp(j))
%             title({TITLE1;TITLE2})
            
        end
        AssemblySimilarityOn = [AssemblySimilarityOn,max(sim_Temp)];
        AssemblySimilarityOn_Alloff = [AssemblySimilarityOn_Alloff,sim_Temp];
        [max_sim,maxInd] = max(sim_Temp);
        
    end
    
    for m = 1:CellAssembly_DelayLR.DelayOff.patNum
        sim_Temp = [];
        pattern_Off = CellAssembly_DelayLR.DelayOff.AssmblWght(:,m);
        for j = 1:CellAssembly_DelayLR.DelayOn.patNum
            pattern_On = CellAssembly_DelayLR.DelayOn.AssmblWght(:,j);
            sim_Temp(j) = dot(pattern_On,pattern_Off)/(norm(pattern_On)*norm(pattern_Off));  
        end
        AssemblySimilarityOff = [AssemblySimilarityOff,max(sim_Temp)];
        AssemblySimilarityOff_Allon = [AssemblySimilarityOff_Allon,sim_Temp];
    end
end

figure
boxplot(AssemblyDensity_On,'Position',1)
hold on
boxplot(AssemblyDensity_Off,'Position',2)
xlim([0 3])
ylim([0 7])

figure
for n = 1:length(sessDirs)
    subplot(2,length(sessDirs),n)
    map = timeMap_oncells.(sessDirs{n});
    [~,peakIdx] = max(map,[],2);
    [~,peak_Sort] = sort(peakIdx);
    cellMapTemp_Sort = map(peak_Sort,:);
    cellMapTemp_Sort_Norm = cellMapTemp_Sort./max(cellMapTemp_Sort,[],2);
    
    imagesc(cellMapTemp_Sort_Norm);hold on
    colormap(jet)
%     vertplot(preBinLength-0.5, 0, size(cellMapTemp_Sort_Norm,1), 'k--');
    TITLE1 = sessDirs{n};
    TITLE2 = 'on assembly cells';
    title({TITLE1;TITLE2},'Interpreter','None')
    axis on
    set(gca, 'xtick', [1 size(cellMapTemp_Sort_Norm,2)]);
    set(gca, 'xticklabels', [0 10]);
    caxis(gca,[0 1])
    
    subplot(2,length(sessDirs),n+length(sessDirs))
    map = timeMap_offcells.(sessDirs{n});
    [~,peakIdx] = max(map,[],2);
    [~,peak_Sort] = sort(peakIdx);
    cellMapTemp_Sort = map(peak_Sort,:);
    cellMapTemp_Sort_Norm = cellMapTemp_Sort./max(cellMapTemp_Sort,[],2);
    
    imagesc(cellMapTemp_Sort_Norm);hold on
    colormap(jet)
%     vertplot(preBinLength-0.5, 0, size(cellMapTemp_Sort_Norm,1), 'k--');
    TITLE1 = sessDirs{n};
    TITLE2 = 'off assembly cells';
    title({TITLE1;TITLE2},'Interpreter','None')
    axis on
    set(gca, 'xtick', [1 size(cellMapTemp_Sort_Norm,2)]);
    set(gca, 'xticklabels', [0 10]);
    caxis(gca,[0 1])
end


figure
cellNumAvg_on = mean(AssemblyCellNum_On);
cellNumSte_on = std(AssemblyCellNum_On)./sqrt(length(AssemblyCellNum_On)-1);

cellNumAvg_off = mean(AssemblyCellNum_Off);
cellNumSte_off = std(AssemblyCellNum_Off)./sqrt(length(AssemblyCellNum_Off)-1);

errorbar(1,cellNumAvg_on,cellNumSte_on,'ro')
hold on
errorbar(2,cellNumAvg_off,cellNumSte_off,'ko')
xlim([0 3])

figure
boxplot(AssemblyCellNum_On,'Position',1)
hold on
boxplot(AssemblyCellNum_Off,'Position',2)
xlim([0 3])
ylim([0 7])

figure
Violin(AssemblySimilarityOn_Alloff,1)
Violin(AssemblySimilarityOff_Allon,2)

% plot wenn of the delay on cell distribution
figure
setListData = {find(delayActive_on10 & SigArmChoice.on10.delayRateLRDiffAll); find(delayActive_off10 & SigArmChoice.off10.delayRateLRDiffAll); ...
    find(delayActive_on30 & SigArmChoice.on30.delayRateLRDiffAll); find(delayActive_off30 & SigArmChoice.off30.delayRateLRDiffAll)};
setLabels = ["on 10"; "off 10"; "on 30"; "off 30"];
h = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);
title('Delay selective assemblies')

figure
setListData = {find(stemActive_on10 & SigArmChoice.on10.stemRateLRDiffAll); find(stemActive_off10 & SigArmChoice.off10.stemRateLRDiffAll); ...
    find(stemActive_on30 & SigArmChoice.on30.stemRateLRDiffAll); find(stemActive_off30 & SigArmChoice.off30.stemRateLRDiffAll)};
setLabels = ["on 10"; "off 10"; "on 30"; "off 30"];
h = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);
title('Stem selective assemblies')

end