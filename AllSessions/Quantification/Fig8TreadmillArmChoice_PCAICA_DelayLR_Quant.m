
function Fig8TreadmillArmChoice_PCAICA_DelayLR_Quant(inFile,AnalyzeSes)
close all
figureCount = 0;
ArmPCAChoice.rate = [];

sessDirs = {'on10','off10','on30','off30'};
assemblyCat = {'on','off'};
armCompareNames = {'delayRateLRDiffAll','returnRateLRDiffAll','stemRateLRDiffAll','ChoRateLRDiffAll'};
shuffle_armCompareNames = {'shf_delayRateLRDiffAll95','shf_returnRateLRDiffAll95','shf_stemRateLRDiffAll95','shf_ChoRateLRDiffAll95'};
assemblyNum_All = 0;

for j = 1:length(sessDirs)
    for n = 1:length(assemblyCat)
        for k = 1:length(armCompareNames)
            % arm compare value
            ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).(armCompareNames{k}) = [];
            ArmPCA_Strength_Choice.(sessDirs{j}).(assemblyCat{n}).(armCompareNames{k}) = [];
            
            ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).(shuffle_armCompareNames{k}) = [];
            ArmPCA_Strength_Choice.(sessDirs{j}).(assemblyCat{n}).(shuffle_armCompareNames{k}) = [];
        end
    end
%     Assembly_Delay_Rate.(sessDirs{j}) = [];
%     Assembly_Stem_Rate.(sessDirs{j}) = [];
%     Assembly_Reward_Rate.(sessDirs{j}) = [];
%     Assembly_Choice_Rate.(sessDirs{j}) = [];
%     AssemblyCellNum = []; 
end

% delayActive = [];
% returnActive = [];
% stemActive = [];
% choiceActive = [];

% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    
    
    Fig8ArmChoice = fullfile(sessInfo(i).mainDir,'Cell Property','Assembly_ArmChoice_DelayLR_25ms.mat');
    load(Fig8ArmChoice);
    
    if isfield(Assembly_ArmChoice_DelayLR,'on10')
        
%         % load assembly file
%         assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_WholeSes-25ms.mat');
%         load(assemblyFile);
%     
%         assemblyRateFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8TreadmillArmRate_SameDelay_Assembly-25ms.mat');
%         load(assemblyRateFile);
%           
%         assemblyNum = length(Assembly_ArmChoice_SameDelay.on10.delayRateLRDiffAll);
%         
%         for m = 1:CellAssembly_WholeSes.patNum
%             AssemblyCellNum = [AssemblyCellNum,length(CellAssembly_WholeSes.AssmblPtrnCellIDs{m})];
%             assemblyNum_All = assemblyNum_All +1;
%             nameTemp = sprintf('%s%d%s%d%s%d','Rat-',CellAssembly_WholeSes.rat,'Day-',CellAssembly_WholeSes.day,'Assembly-',m);
%             AssemblyName{assemblyNum_All} = nameTemp;
%         end
%   
        for j = 1:length(sessDirs)
            for n = 1:length(assemblyCat)
                %             % delay Rate
                %             Assembly_Delay_Rate.(sessDirs{j}) = [Assembly_Delay_Rate.(sessDirs{j}),...
                %                 Fig8TreadmillArmRate_SameDelay_Assembly.(sessDirs{j}).rateDelay];
                %             Assembly_Stem_Rate.(sessDirs{j}) = [Assembly_Stem_Rate.(sessDirs{j}),...
                %                 Fig8TreadmillArmRate_SameDelay_Assembly.(sessDirs{j}).rateStem];
                %             Assembly_Reward_Rate.(sessDirs{j}) = [Assembly_Reward_Rate.(sessDirs{j}),...
                %                 Fig8TreadmillArmRate_SameDelay_Assembly.(sessDirs{j}).rateReward];
                %             Assembly_Choice_Rate.(sessDirs{j}) = [Assembly_Choice_Rate.(sessDirs{j}),...
                %                 Fig8TreadmillArmRate_SameDelay_Assembly.(sessDirs{j}).rateChoice];
                
                % rate comparison
                ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).delayRateLRDiffAll = [ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).delayRateLRDiffAll,...
                    Assembly_ArmChoice_DelayLR.(sessDirs{j}).(assemblyCat{n}).delayRateLRDiffAll];
                
                ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).stemRateLRDiffAll = [ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).stemRateLRDiffAll,...
                    Assembly_ArmChoice_DelayLR.(sessDirs{j}).(assemblyCat{n}).stemRateLRDiffAll];
                
                ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).ChoRateLRDiffAll = [ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).ChoRateLRDiffAll,...
                    Assembly_ArmChoice_DelayLR.(sessDirs{j}).(assemblyCat{n}).stemChoRateLRDiffAll];
                
                ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).returnRateLRDiffAll = [ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).returnRateLRDiffAll,...
                    Assembly_ArmChoice_DelayLR.(sessDirs{j}).(assemblyCat{n}).returnRateLRDiffAll];
                
                % 95% value
                
                ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).shf_delayRateLRDiffAll95 = [ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).shf_delayRateLRDiffAll95,...
                    Assembly_ArmChoice_DelayLR.(sessDirs{j}).(assemblyCat{n}).shf_delayRateLRDiffAll95];
                
                ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).shf_stemRateLRDiffAll95 = [ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).shf_stemRateLRDiffAll95,...
                    Assembly_ArmChoice_DelayLR.(sessDirs{j}).(assemblyCat{n}).shf_stemRateLRDiffAll95];
                
                ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).shf_ChoRateLRDiffAll95 = [ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).shf_ChoRateLRDiffAll95,...
                    Assembly_ArmChoice_DelayLR.(sessDirs{j}).(assemblyCat{n}).shf_ChoRateLRDiffAll95];
                
                ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).shf_returnRateLRDiffAll95 = [ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).shf_returnRateLRDiffAll95,...
                    Assembly_ArmChoice_DelayLR.(sessDirs{j}).(assemblyCat{n}).shf_returnRateLRDiffAll95];
            end
        end
    end
        
end


for n = 1:length(assemblyCat)
    for k = 1:length(armCompareNames)
        for j = 1:length(sessDirs)            
            % All trials
            rateComp = ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).(armCompareNames{k});
            shfValue = ArmPCAChoice.(sessDirs{j}).(assemblyCat{n}).(shuffle_armCompareNames{k});
            SigArmChoice.(sessDirs{j}).(assemblyCat{n}).(armCompareNames{k}) = (abs(rateComp)-abs(shfValue))>0;
        end
    end
end
    
    figure
    bar(1,100*sum(SigArmChoice.on10.on.returnRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on10.on.shf_returnRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    hold on
    bar(2,100*sum(SigArmChoice.off10.on.returnRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off10.on.shf_returnRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    bar(3,100*sum(SigArmChoice.on30.on.returnRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on30.on.shf_returnRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    bar(4,100*sum(SigArmChoice.off30.on.returnRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off30.on.shf_returnRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    
    bar(6,100*sum(SigArmChoice.on10.on.delayRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on10.on.shf_delayRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    hold on
    bar(7,100*sum(SigArmChoice.off10.on.delayRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off10.on.shf_delayRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    bar(8,100*sum(SigArmChoice.on30.on.delayRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on30.on.shf_delayRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    bar(9,100*sum(SigArmChoice.off30.on.delayRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off30.on.shf_delayRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    
    bar(11,100*sum(SigArmChoice.on10.on.stemRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on10.on.shf_stemRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    hold on
    bar(12,100*sum(SigArmChoice.off10.on.stemRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off10.on.shf_stemRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    bar(13,100*sum(SigArmChoice.on30.on.stemRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on30.on.shf_stemRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    bar(14,100*sum(SigArmChoice.off30.on.stemRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off30.on.shf_stemRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    
    bar(16,100*sum(SigArmChoice.on10.on.ChoRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on10.on.shf_ChoRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    hold on
    bar(17,100*sum(SigArmChoice.off10.on.ChoRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off10.on.shf_ChoRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    bar(18,100*sum(SigArmChoice.on30.on.ChoRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on30.on.shf_ChoRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    bar(19,100*sum(SigArmChoice.off30.on.ChoRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off30.on.shf_ChoRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    
    set(gca,'XTick',[1,2,3,4,6,7,8,9,11,12,13,14,16:19],'XTickLabel',{'on10','off10','on30','off30'});

    
    figure
    bar(1,100*sum(SigArmChoice.on10.off.returnRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on10.off.shf_returnRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    hold on
    bar(2,100*sum(SigArmChoice.off10.off.returnRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off10.off.shf_returnRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    bar(3,100*sum(SigArmChoice.on30.off.returnRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on30.off.shf_returnRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    bar(4,100*sum(SigArmChoice.off30.off.returnRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off30.off.shf_returnRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    
    bar(6,100*sum(SigArmChoice.on10.off.delayRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on10.off.shf_delayRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    hold on
    bar(7,100*sum(SigArmChoice.off10.off.delayRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off10.off.shf_delayRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    bar(8,100*sum(SigArmChoice.on30.off.delayRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on30.off.shf_delayRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    bar(9,100*sum(SigArmChoice.off30.off.delayRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off30.off.shf_delayRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    
    bar(11,100*sum(SigArmChoice.on10.off.stemRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on10.off.shf_stemRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    hold on
    bar(12,100*sum(SigArmChoice.off10.off.stemRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off10.off.shf_stemRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    bar(13,100*sum(SigArmChoice.on30.off.stemRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on30.off.shf_stemRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    bar(14,100*sum(SigArmChoice.off30.off.stemRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off30.off.shf_stemRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    
    bar(16,100*sum(SigArmChoice.on10.off.ChoRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on10.off.shf_ChoRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    hold on
    bar(17,100*sum(SigArmChoice.off10.off.ChoRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off10.off.shf_ChoRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    bar(18,100*sum(SigArmChoice.on30.off.ChoRateLRDiffAll)/sum(~isnan(ArmPCAChoice.on30.off.shf_ChoRateLRDiffAll95)),'FaceColor',[1,0,0.2]);
    bar(19,100*sum(SigArmChoice.off30.off.ChoRateLRDiffAll)/sum(~isnan(ArmPCAChoice.off30.off.shf_ChoRateLRDiffAll95)),'FaceColor',[0.6,0.6,0.6]);
    
    set(gca,'XTick',[1,2,3,4,6,7,8,9,11,12,13,14,16:19],'XTickLabel',{'on10','off10','on30','off30'});
end