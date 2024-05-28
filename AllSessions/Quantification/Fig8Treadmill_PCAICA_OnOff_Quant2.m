
function Fig8Treadmill_PCAICA_OnOff_Quant2(inFile,AnalyzeSes)
close all

p.avgRateThres = 0.5;
% % Read in input information
sessInfo = SessInfoImport(inFile);

assembly_ReturnSelective.on = [];
assembly_ChoiceSelective.on = [];
assembly_ReturnSelective.off = [];
assembly_ChoiceSelective.off = [];
assembly_ReturnSelectiveChance = [];
assembly_ChoiceSelectiveChance = [];

for i = AnalyzeSes(1:end)  
    Delay_OnOff_File = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_DelayOnOff_LR-25ms.mat');
    load(Delay_OnOff_File);
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    Fig8ArmChoice = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8TreadmillArmChoice_SameDelay.mat');
    armChoice = load(Fig8ArmChoice);  
    armRateFile = fullfile(sessInfo(i).mainDir,'Cell Property', 'Fig8TreadmillArmRate.mat');
    load(armRateFile); 
    
    % get each phase names (no delay etc)
    sessDirs = {'on10','off10','on30','off30'};
    for j = 1:length(sessDirs)
        returnSelection.(sessDirs{j}) = [];
        choiceSelection.(sessDirs{j}) = [];
    end
    rate_Return = [];
    return_onvalidCell = [];
    rate_Choice = [];
    choice_onvalidCell = [];
    

    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;   
    validCell = find(rateLabel==1);
    
    %% determine whether this cell is L / R specific cell on return and
    % choice arm
    for j = 1:length(sessDirs)
        for k = 1:length(rateLabel)
            if rateLabel(k) == 1
                rateComp = armChoice.Fig8TreadmillArmChoice_SameDelay.(sessDirs{j}).returnRateLRDiffAll(k);
                shfValue = armChoice.Fig8TreadmillArmChoice_SameDelay.(sessDirs{j}).shf_returnRateLRDiffAll95(k);
                if rateComp>0 && shfValue > 0 && ((abs(rateComp)-abs(shfValue))>0)
                    selectTemp = 1;
                elseif rateComp<0 && shfValue < 0 && ((abs(rateComp)-abs(shfValue))>0)
                    selectTemp = -1;
                else
                    selectTemp = 0;
                end
                returnSelection.(sessDirs{j}) = [returnSelection.(sessDirs{j});selectTemp];
                
                rateComp = armChoice.Fig8TreadmillArmChoice_SameDelay.(sessDirs{j}).ChoRateLRDiffAll(k);
                shfValue = armChoice.Fig8TreadmillArmChoice_SameDelay.(sessDirs{j}).shf_ChoRateLRDiffAll95(k);
                if rateComp>0 && shfValue > 0 && ((abs(rateComp)-abs(shfValue))>0)
                    selectTemp = 1;
                elseif rateComp<0 && shfValue < 0 && ((abs(rateComp)-abs(shfValue))>0)
                    selectTemp = -1;
                else
                    selectTemp = 0;
                end
                choiceSelection.(sessDirs{j}) = [choiceSelection.(sessDirs{j});selectTemp];
                
            end
        end
    end
    
    rate_Temp = [Fig8TreadmillArmRate.on10.rateReturn(validCell,1)';Fig8TreadmillArmRate.off10.rateReturn(validCell,1)';...
        Fig8TreadmillArmRate.on30.rateReturn(validCell,1)';Fig8TreadmillArmRate.off30.rateReturn(validCell,1)'];
    rate_Return = [rate_Return,rate_Temp];
    return_onvalidCellTemp = (sum(rate_Temp>p.avgRateThres,1))>0;
    return_onvalidCell = [return_onvalidCell,return_onvalidCellTemp];
    rate_Temp = [Fig8TreadmillArmRate.on10.rateChoice(validCell,1)';Fig8TreadmillArmRate.off10.rateChoice(validCell,1)';...
        Fig8TreadmillArmRate.on30.rateChoice(validCell,1)';Fig8TreadmillArmRate.off30.rateChoice(validCell,1)'];
    rate_Choice = [rate_Choice,rate_Temp];
    choice_onvalidCellTemp = (sum(rate_Temp>p.avgRateThres,1))>0;
    choice_onvalidCell = [choice_onvalidCell,choice_onvalidCellTemp];
    
    returnTemp = [returnSelection.on10,returnSelection.off10,returnSelection.on30,returnSelection.off30];
    return_L = (sum(returnTemp == 1,2)) > 0;
    return_R = (sum(returnTemp == -1,2)) > 0;
    choiceTemp = [choiceSelection.on10,choiceSelection.off10,choiceSelection.on30,choiceSelection.off30];
    choice_L = (sum(choiceTemp == 1,2)) > 0;
    choice_R = (sum(choiceTemp == -1,2)) > 0;
    
    %%
    assembly_return_L = 0;
    assembly_return_R = 0;
    assembly_choice_L = 0;
    assembly_choice_R = 0;
    assembly_cellNum = 0;
    for m = 1:CellAssembly_DelayLR.DelayOff.patNum
       assemblyID = CellAssembly_DelayLR.DelayOff.AssmblPtrnCellIDs{m};
       assembly_return_L = assembly_return_L + sum(return_L(assemblyID));
       assembly_return_R = assembly_return_R + sum(return_R(assemblyID));
       assembly_choice_L = assembly_choice_L + sum(choice_L(assemblyID));
       assembly_choice_R = assembly_choice_R + sum(choice_R(assemblyID));
       assembly_cellNum = assembly_cellNum + length(assemblyID);
    end
    returnSelective_off = (assembly_return_L+assembly_return_R)/assembly_cellNum;
    returnSelectiveBase = sum(return_L+return_R)/length(return_L);
    choiceSelective_off = (assembly_choice_L+assembly_choice_R)/assembly_cellNum;
    choiceSelectiveBase = sum(choice_L+choice_R)/length(return_L);
    
    
    assembly_return_L = 0;
    assembly_return_R = 0;
    assembly_choice_L = 0;
    assembly_choice_R = 0;
    assembly_cellNum = 0;
    for m = 1:CellAssembly_DelayLR.DelayOn.patNum
       assemblyID = CellAssembly_DelayLR.DelayOn.AssmblPtrnCellIDs{m};
       assembly_return_L = assembly_return_L + sum(return_L(assemblyID));
       assembly_return_R = assembly_return_R + sum(return_R(assemblyID));
       assembly_choice_L = assembly_choice_L + sum(choice_L(assemblyID));
       assembly_choice_R = assembly_choice_R + sum(choice_R(assemblyID));
       assembly_cellNum = assembly_cellNum + length(assemblyID);
    end
    returnSelective_on = (assembly_return_L+assembly_return_R)/assembly_cellNum;
    choiceSelective_on = (assembly_choice_L+assembly_choice_R)/assembly_cellNum;
    
    assembly_ReturnSelective.on = [assembly_ReturnSelective.on,returnSelective_on];
    assembly_ChoiceSelective.on = [assembly_ChoiceSelective.on,choiceSelective_on];
    assembly_ReturnSelective.off = [assembly_ReturnSelective.off,returnSelective_off];
    assembly_ChoiceSelective.off = [assembly_ChoiceSelective.off,choiceSelective_off];
    assembly_ReturnSelectiveChance = [assembly_ReturnSelectiveChance,returnSelectiveBase];
    assembly_ChoiceSelectiveChance = [assembly_ChoiceSelectiveChance,choiceSelectiveBase];

%     
%     ReturnSelective_off_shuffle = zeros(100,1);
%     for n = 1:100
%         assembly_return_L = zeros(CellAssembly_DelayLR.DelayOff.patNum,1);
%         assembly_return_R = zeros(CellAssembly_DelayLR.DelayOff.patNum,1);
%         for m = 1:CellAssembly_DelayLR.DelayOff.patNum
%             cellNum = length(CellAssembly_DelayLR.DelayOff.AssmblPtrnCellIDs{m});
%             assemblyID = ceil(rand(1,cellNum)*sum(rateLabel));
%             if sum(return_L(assemblyID)==1) >= 1
%                 assembly_return_L(m) = 1;
%             elseif sum(return_R(assemblyID)==1) >= 1
%                 assembly_return_R(m) = 1;
%             end
%         end
%         ReturnSelective_temp = sum((assembly_return_L + assembly_return_R)>=1)/CellAssembly_DelayLR.DelayOff.patNum;
%         ReturnSelective_off_shuffle(n) = ReturnSelective_temp;
%     end
    
end
figure
subplot(1,2,1)
plot([assembly_ReturnSelective.on;assembly_ReturnSelectiveChance]./assembly_ReturnSelectiveChance,'ko--')
xlim([0 3])
ylim([0 4])
subplot(1,2,2)
plot([assembly_ReturnSelective.off;assembly_ReturnSelectiveChance]./assembly_ReturnSelectiveChance,'ko--')
xlim([0 3])
ylim([0 4])

figure
subplot(1,2,1)
plot([assembly_ChoiceSelective.on;assembly_ChoiceSelectiveChance]./assembly_ChoiceSelectiveChance,'ko--')
xlim([0 3])
ylim([0 4])
subplot(1,2,2)
plot([assembly_ChoiceSelective.off;assembly_ChoiceSelectiveChance]./assembly_ChoiceSelectiveChance,'ko--')
xlim([0 3])
ylim([0 4])

figure
Violin(assembly_ReturnSelective.on./assembly_ReturnSelectiveChance,1)
Violin(assembly_ReturnSelective.off./assembly_ReturnSelectiveChance,2)
Violin(assembly_ChoiceSelective.on./assembly_ChoiceSelectiveChance,3)
Violin(assembly_ChoiceSelective.off./assembly_ChoiceSelectiveChance,4)


% figure
% Violin(AssemblySimilarityOn,1)
% hold on
% Violin(AssemblySimilarityOff,2)


figure
for n = 1:length(sessDirs)
    subplot(2,length(sessDirs),n)
    map = timeMap_oncells.(sessDirs{n});
    [~,peakvalidCell] = max(map,[],2);
    [~,peak_Sort] = sort(peakvalidCell);
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
    [~,peakvalidCell] = max(map,[],2);
    [~,peak_Sort] = sort(peakvalidCell);
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
    [~,peakvalidCell] = max(map,[],2);
    [~,peak_Sort] = sort(peakvalidCell);
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
    [~,peakvalidCell] = max(map,[],2);
    [~,peak_Sort] = sort(peakvalidCell);
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