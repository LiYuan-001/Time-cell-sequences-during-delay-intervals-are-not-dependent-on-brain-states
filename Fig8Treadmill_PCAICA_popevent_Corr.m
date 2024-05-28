function Fig8Treadmill_PCAICA_popevent_Corr(inFile,AnalyzeSes)

% set parameters for analysis
p.savePlot = 0;
p.writeToFile = 0;
binWidth = 25/1000; % unit sec

% Read in input information
sessInfo = SessInfoImport(inFile);
sessDirs2 = {'off10_1','off30_1','off10_2','off30_2'};

AssemblySimilarityOn = [];
AssemblySimilarityOff = [];
assembly_ID_on = [];
assembly_ID_off = [];

on_corr = [];
off_corr = [];

for i = AnalyzeSes(1:end)
    close all
    % selectivity for the whole map
    Assembly_Pop_Corr.session = sessInfo(i).mainDir;
    Assembly_Pop_Corr.rat = sessInfo(i).animal;
    Assembly_Pop_Corr.day = sessInfo(i).day;
    
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\Assembly_ArmChoice');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end
    
    if p.writeToFile
        savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
        if ~exist(savedir2, 'dir')
            mkdir(savedir2);
        end
    end
    
    display(['About to do session ' sessInfo(i).mainDir]);
    % get each phase names (no delay etc)
%     sessDirs = sessInfo(i).sessDirs;
    
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load assembly file
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_DelayOnOff_LR-25ms.mat');
    load(assemblyFile);
    % load population events file
    popFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_Delay.mat');
    load(popFile);
    
    
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    binTime = CellAssembly_DelayLR.DelayOn.binTime;
    
    patNum_on = CellAssembly_DelayLR.DelayOn.patNum; 
    AssmblPtrnCellIDs_on = CellAssembly_DelayLR.DelayOn.AssmblPtrnCellIDs;
    AssmblWght_on = CellAssembly_DelayLR.DelayOn.AssmblWght;
    AssmblStrength_on = CellAssembly_DelayLR.DelayOn.AssmblStrength;
    event_Time_on = CellAssembly_DelayLR.DelayOn.event_Time;
    event_strength_on = CellAssembly_DelayLR.DelayOn.event_strength;
    event_Num_on = CellAssembly_DelayLR.DelayOn.event_Num;
    
    patNum_off = CellAssembly_DelayLR.DelayOff.patNum; 
    AssmblPtrnCellIDs_off = CellAssembly_DelayLR.DelayOff.AssmblPtrnCellIDs;
    AssmblWght_off = CellAssembly_DelayLR.DelayOff.AssmblWght;
    AssmblStrength_off = CellAssembly_DelayLR.DelayOff.AssmblStrength;
    event_Time_off = CellAssembly_DelayLR.DelayOff.event_Time;
    event_strength_off = CellAssembly_DelayLR.DelayOff.event_strength;
    event_Num_off = CellAssembly_DelayLR.DelayOff.event_Num;
    
    popLabel = [];
    assemblyLabel_off = cell(1,patNum_off);
    assemblyLabel_on = cell(1,patNum_on);
    
    if sum(rateLabel) >= 20
        for m = 1:patNum_on
            sim_Temp = [];
            pattern_On = AssmblWght_on(:,m);
            for k = 1:patNum_off
                pattern_Off = AssmblWght_off(:,k);
                sim_Temp(k) = dot(pattern_On,pattern_Off)/(norm(pattern_On)*norm(pattern_Off));
            end
            AssemblySimilarityOn = [AssemblySimilarityOn,max(sim_Temp)];
            assembly_ID_on = [assembly_ID_on,CellAssembly_DelayLR.rat*100+CellAssembly_DelayLR.day*10+m];
        end
        
        for m = 1:patNum_off
            sim_Temp = [];
            pattern_Off = AssmblWght_off(:,m);
            for k = 1:patNum_on
                pattern_On = AssmblWght_on(:,k);
                sim_Temp(k) = dot(pattern_Off,pattern_On)/(norm(pattern_Off)*norm(pattern_On));
            end
            AssemblySimilarityOff = [AssemblySimilarityOff,max(sim_Temp)];
            assembly_ID_off = [assembly_ID_off,CellAssembly_DelayLR.rat*100+CellAssembly_DelayLR.day*10+m];
        end
        
        for j = 1:length(sessDirs2)
            eventStartTs = DelayPopFire_Delay.(sessDirs2{j}).eventStartTs;
            eventEndTs = DelayPopFire_Delay.(sessDirs2{j}).eventEndTs;
            
            % load maps for each cluster
            pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs2{j}, 'PathZone.mat');
            load(pathZoneFile);
            delayFile = fullfile(sessInfo(i).mainDir,sessDirs2{j}, 'Fig8DelayZonePos.mat');
            load(delayFile);
            
            % def1: delay starts at barrier
            % def2: delay starts at entrance
            delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
            delayTend1 = Fig8DelayZonePos.delayPos1.endT;
            %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
            %         delayTend2 = Fig8DelayZonePos.delayPos2.endT;
            
            trialNum = size(delayTstart1,2);
            if strcmp(sessDirs2{j}(end-3:end-2),'10')
                maxT = 10;
                delayTend1_2 = delayTstart1+maxT;
            elseif strcmp(sessDirs2{j}(end-3:end-2),'30')
                maxT = 30;
                delayTend1_2 = delayTstart1+maxT;
            else
                error('Delay time is wrong')
            end
            
            for k = 1:trialNum
                binTemp = delayTstart1(k):binWidth:delayTend1_2(k);
                popLabelTemp = zeros(1,length(binTemp));
                
                for mm = 1:length(binTemp)-1
                    % population event label
                    if any(binTemp(mm) >= eventStartTs(:) & binTemp(mm)<= eventEndTs(:))
                        popLabelTemp(mm) = 1;
                    end
                end
                popLabel = [popLabel,popLabelTemp];
                
                % assemblyEvent label
                for nn = 1:patNum_off
                    assemblyLabelTemp = zeros(1,length(binTemp));
                    eventTs = event_Time_off{nn};
                    for mm = 1:length(binTemp)-1
                        if any(eventTs>=binTemp(mm) & eventTs<binTemp(mm+1))
                            assemblyLabelTemp(mm) = 1;
                        end
                    end
                    assemblyLabel_off{nn} = [assemblyLabel_off{nn},assemblyLabelTemp];
                end
                
                % assemblyEvent label
                for nn = 1:patNum_on
                    assemblyLabelTemp = zeros(1,length(binTemp));
                    eventTs = event_Time_on{nn};
                    for mm = 1:length(binTemp)-1
                        if any(eventTs>=binTemp(mm) & eventTs<binTemp(mm+1))
                            assemblyLabelTemp(mm) = 1;
                        end
                    end
                    assemblyLabel_on{nn} = [assemblyLabel_on{nn},assemblyLabelTemp];
                end
                
            end
        end
        for nn = 1:patNum_off
            corrVal = corr(assemblyLabel_off{nn}',popLabel');
            off_corr = [off_corr,corrVal];
        end
        for nn = 1:patNum_on
            corrVal = corr(assemblyLabel_on{nn}',popLabel');
            on_corr = [on_corr,corrVal];
        end
    end    
end

figure
Violin(on_corr,1,'ShowData',false);
Violin(off_corr,2,'ShowData',false);
end