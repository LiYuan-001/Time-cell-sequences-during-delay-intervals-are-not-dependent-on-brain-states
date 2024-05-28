close all
sessInfo = SessInfoImport('Fig8Treadmill_OnOff.xlsx');
% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'
p.spikeFs = 32000;
for i = 3 % use 1043 day 6 or 1044 day 3
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    sessDirs = sessInfo(i).sessDirs;
    % load assembly file
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_WholeSes-25ms.mat');
    load(assemblyFile);
    % load spikes from each main session
    % get event timestamp
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
        
    cellNum = CellAssembly_WholeSes.ValidCellNum;
    cellInd = find(CellAssembly_WholeSes.ValidCellLabel == 1);
    binTime = CellAssembly_WholeSes.binTime;
    
    for j = 1
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        delayTstart = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend = delayTstart + 10;
        
        % plot assembly weight
        figure
        AssemblyWeight = CellAssembly_WholeSes.AssmblWght(:,1);
        AssmblPtrnCellIDs = CellAssembly_WholeSes.AssmblPtrnCellIDs{1};
        h1 = subplot(3,6,1,'align');
        stem(AssemblyWeight,'k')
        hold on
        stem(AssmblPtrnCellIDs,AssemblyWeight(AssmblPtrnCellIDs),'r')
        set(gca,'XDir','reverse')
        xlim([0.5 cellNum+0.5])
        view([90 -90])
        ylim([-1 1])

        % plot another assembly
        AssemblyWeight = CellAssembly_WholeSes.AssmblWght(:,4);
        AssmblPtrnCellIDs = CellAssembly_WholeSes.AssmblPtrnCellIDs{4};
        h1 = subplot(3,6,2,'align');
        stem(AssemblyWeight,'k')
        hold on
        stem(AssmblPtrnCellIDs,AssemblyWeight(AssmblPtrnCellIDs),'b')
        view([90 -90])
        ylim([-1 1])

        
        % trial 7
        startTime = delayTstart(1);
        endTime = delayTend(1);
        binInd = binTime>=startTime & binTime<=endTime;
        
        % plot strength
        subplot(3,3,3)
        AssmblStrength = CellAssembly_WholeSes.AssmblStrength(1,binInd);
        plot(AssmblStrength,'r');
        hold on
        AssmblStrength = CellAssembly_WholeSes.AssmblStrength(4,binInd);
        plot(AssmblStrength,'b');
        
        % plot spikes
        subplot(3,3,2)
        AssmblPtrnCellIDs = CellAssembly_WholeSes.AssmblPtrnCellIDs{1};
        for k = 1:length(cellInd)
            cellInd2 = cellInd(k);
            tsp = Spike_Session.(sessDirs{j}){cellInd2};
            tsp2 = tsp(tsp>=startTime & tsp<=endTime)-startTime;
            
            xPoints = [tsp2';tsp2'];
            yPoints = [k+zeros(size(tsp2'))-0.3;k+zeros(size(tsp2'))+0.3];               
            if any(AssmblPtrnCellIDs == k)               
                if ~isempty(tsp2)
                    plot(xPoints,yPoints,'r')
                    hold on
                end
            else
                if ~isempty(tsp2)
                    plot(xPoints,yPoints,'Color',[0.5,0.5,0.5])
                    hold on
                end
            end
        end
        xlim([0 10])
        ylim([0.5 cellNum+0.5])
        set(gca,'YDir','Normal')
        hold on
        
        AssmblPtrnCellIDs = CellAssembly_WholeSes.AssmblPtrnCellIDs{4};
        for k = 1:length(cellInd)
            cellInd2 = cellInd(k);
            tsp = Spike_Session.(sessDirs{j}){cellInd2};
            tsp2 = tsp(tsp>=startTime & tsp<=endTime)-startTime;
            
            xPoints = [tsp2';tsp2'];
            yPoints = [k+zeros(size(tsp2'))-0.3;k+zeros(size(tsp2'))+0.3];               
            if any(AssmblPtrnCellIDs == k)               
                if ~isempty(tsp2)
                    plot(xPoints,yPoints,'b')
                    hold on
                end
%             else
%                 if ~isempty(tsp2)
%                     plot(xPoints,yPoints,'Color',[0.5,0.5,0.5])
%                     hold on
%                 end
            end
        end       
                
        
%         h2 = subplot(6,3,2,'align');
%         imagesc(delay.seqTotal_L(:,800:1200))
%         title('Left delay strength')
%         h3 = subplot(6,3,5,'align');
%         plot(delay.AssmblStrength_LL(1,800:1200))
%         h4 = subplot(6,3,3,'align');
%         imagesc(delay.seqTotal_R(:,800:1200))
%         title('Right delay strength')
%         h5 = subplot(6,3,6,'align');
%         plot(delay.AssmblStrength_LR(1,800:1200))
%         linkaxes([h3,h5],'y');
%         ylim([-2 30])
%         linkaxes([h2,h3,h4,h5],'x');
    end
    
    for j = 6
        delay = CellAseembly.(sessDirs{j}).delay;
        reward = CellAseembly.(sessDirs{j}).reward;
        h1 = subplot(6,6,14,'align');
        stem(delay.AssmblWght_L(:,1))
        hold on
        stem(delay.AssmblPtrnCellIDs_L{1},delay.AssmblWght_L(delay.AssmblPtrnCellIDs_L{1},1),'r')
        set(gca,'XDir','reverse')
        xlim([0.5 cellNum+0.5])
        view([90 -90])
        %                     xlim([0.5 35.5])
        title('off10 Delay L assemble example')
        h2 = subplot(6,3,8,'align');
        imagesc(delay.seqTotal_L(:,800:1200))
        title('Left delay strength')
        h3 = subplot(6,3,11,'align');
        plot(delay.AssmblStrength_LL(1,800:1200))
        h4 = subplot(6,3,9,'align');
        imagesc(delay.seqTotal_R(:,800:1200))
        title('Right delay strength')
        h5 = subplot(6,3,12,'align');
        plot(delay.AssmblStrength_LR(1,800:1200))
        linkaxes([h3,h5],'y');
        ylim([-2 30])
        linkaxes([h2,h3,h4,h5],'x');
    end
    for j = 1
        delay = CellAseembly.(sessDirs{j}).delay;
        reward = CellAseembly.(sessDirs{j}).reward;
        h1 = subplot(6,6,26,'align');
        stem(reward.AssmblWght_L(:,3))
        hold on
        stem(reward.AssmblPtrnCellIDs_L{3},reward.AssmblWght_L(reward.AssmblPtrnCellIDs_L{3},3),'r')
        set(gca,'XDir','reverse')
        xlim([0.5 cellNum+0.5])
        view([90 -90])
        %                     xlim([0.5 35.5])
        title('on10 Reward L assemble example')
        h2 = subplot(6,3,14,'align');
        imagesc(reward.seqTotal_L(:,1100:1500))
        title('Left reward strength')
        h3 = subplot(6,3,17,'align');
        plot(reward.AssmblStrength_LL(3,1100:1500))
        h4 = subplot(6,3,15,'align');
        imagesc(reward.seqTotal_R(:,1100:1500))
        title('Right reward strength')
        h5 = subplot(6,3,18,'align');
        plot(reward.AssmblStrength_LR(3,1100:1500))
        linkaxes([h3,h5],'y');
        ylim([-2 30])
        linkaxes([h2,h3,h4,h5],'x');
    end
    
end