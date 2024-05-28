
function Fig8Treadmill_PCAICA_OnOff_Quant(inFile,AnalyzeSes)
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

AssemblyNum_On = [];
AssemblyNum_Off = [];

pyrNumAll = [];

% get each phase names (no delay etc)
sessDirs = {'on10','off10','on30','off30'};
for n = 1:length(sessDirs)
    timeMap_oncells.(sessDirs{n}) = [];
    timeMap_offcells.(sessDirs{n}) = [];
    
    timeMap_onAssembly.(sessDirs{n}) = [];
    timeMap_offAssembly.(sessDirs{n}) = [];
end

for i = AnalyzeSes(1:end)
    
    
    Delay_OnOff_File = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_DelayOnOff_LR-25ms.mat');
    load(Delay_OnOff_File);
    % load delayFire map
    delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap_2Session.mat');
    load(delayFile);
    % load assembly delay firing map
    delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayFire_Assembly_onoff.mat');
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
    pyrNumAll = [pyrNumAll,cellNum];
    
    AssemblyDensity_On = [AssemblyDensity_On,cellNum/CellAssembly_DelayLR.DelayOn.patNum];
    AssemblyDensity_Off = [AssemblyDensity_Off,cellNum/CellAssembly_DelayLR.DelayOff.patNum];
    
    AssemblyNum_On = [AssemblyNum_On,CellAssembly_DelayLR.DelayOn.patNum];        
    AssemblyNum_Off = [AssemblyNum_Off,CellAssembly_DelayLR.DelayOff.patNum];
    
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
        end
        AssemblySimilarityOn = [AssemblySimilarityOn,max(sim_Temp)];
        AssemblySimilarityOn_Alloff = [AssemblySimilarityOn_Alloff,sim_Temp];
        [max_sim,maxInd] = max(sim_Temp);    
        
        for n = 1:length(sessDirs)
            timeMap_onAssembly.(sessDirs{n}) = [timeMap_onAssembly.(sessDirs{n});DelayFire_Assembly_onoff.on.(sessDirs{n}).spikeRate1_Combined_Smooth{m}];
        end
    end
    
%     stem(AssemblyWeight,'k')
%         hold on
%         stem(AssmblPtrnCellIDs,AssemblyWeight(AssmblPtrnCellIDs),'r')
%         set(gca,'XDir','reverse')
%         xlim([0.5 cellNum+0.5])
%         view([90 -90])
        
    for m = 1:CellAssembly_DelayLR.DelayOff.patNum
        sim_Temp = [];
        pattern_Off = CellAssembly_DelayLR.DelayOff.AssmblWght(:,m);
        for j = 1:CellAssembly_DelayLR.DelayOn.patNum
            pattern_On = CellAssembly_DelayLR.DelayOn.AssmblWght(:,j);
            sim_Temp(j) = dot(pattern_On,pattern_Off)/(norm(pattern_On)*norm(pattern_Off));  
        end
        AssemblySimilarityOff = [AssemblySimilarityOff,max(sim_Temp)];
        AssemblySimilarityOff_Allon = [AssemblySimilarityOff_Allon,sim_Temp];
        
        for n = 1:length(sessDirs)
            timeMap_offAssembly.(sessDirs{n}) = [timeMap_offAssembly.(sessDirs{n});DelayFire_Assembly_onoff.off.(sessDirs{n}).spikeRate1_Combined_Smooth{m}];
        end
    end
end

figure
Violin(AssemblySimilarityOn,1)
hold on
Violin(AssemblySimilarityOff,2)


% get 2d hist for pat num in pyr cell num
xbin = [1:2:100];
ybin = [0:1:20];
[count_2d.imec0.on,xedges,yedges] = histcounts2(pyrNumAll,AssemblyNum_On,xbin,ybin);
[count_2d.imec0.off,xedges,yedges] = histcounts2(pyrNumAll,AssemblyNum_Off,xbin,ybin);

figure
imagesc(xedges(2:end),yedges(2:end),count_2d.imec0.on')
axis xy
caxis([0 4])
figure
imagesc(xedges(2:end),yedges(2:end),count_2d.imec0.off')
axis xy
caxis([0 4])


figure
boxplot(AssemblyCellNum_On,'Position',1);
hold on
boxplot(AssemblyCellNum_Off,'Position',2);

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
for n = 1:length(sessDirs)
    subplot(2,length(sessDirs),n)
    map = timeMap_onAssembly.(sessDirs{n});
    [~,peakIdx] = max(map,[],2);
    [~,peak_Sort] = sort(peakIdx);
    cellMapTemp_Sort = map(peak_Sort,:);
    cellMapTemp_Sort_Norm = cellMapTemp_Sort./max(cellMapTemp_Sort,[],2);
    
    imagesc(cellMapTemp_Sort_Norm);hold on
    colormap(jet)
%     vertplot(preBinLength-0.5, 0, size(cellMapTemp_Sort_Norm,1), 'k--');
    TITLE1 = sessDirs{n};
    TITLE2 = 'on assembly map';
    title({TITLE1;TITLE2},'Interpreter','None')
    axis on
    set(gca, 'xtick', [1 size(cellMapTemp_Sort_Norm,2)]);
    set(gca, 'xticklabels', [0 10]);
    caxis(gca,[0 1])
    
    subplot(2,length(sessDirs),n+length(sessDirs))
    map = timeMap_offAssembly.(sessDirs{n});
    [~,peakIdx] = max(map,[],2);
    [~,peak_Sort] = sort(peakIdx);
    cellMapTemp_Sort = map(peak_Sort,:);
    cellMapTemp_Sort_Norm = cellMapTemp_Sort./max(cellMapTemp_Sort,[],2);
    
    imagesc(cellMapTemp_Sort_Norm);hold on
    colormap(jet)
%     vertplot(preBinLength-0.5, 0, size(cellMapTemp_Sort_Norm,1), 'k--');
    TITLE1 = sessDirs{n};
    TITLE2 = 'off assembly map';
    title({TITLE1;TITLE2},'Interpreter','None')
    axis on
    set(gca, 'xtick', [1 size(cellMapTemp_Sort_Norm,2)]);
    set(gca, 'xticklabels', [0 10]);
    caxis(gca,[0 1])
end

end