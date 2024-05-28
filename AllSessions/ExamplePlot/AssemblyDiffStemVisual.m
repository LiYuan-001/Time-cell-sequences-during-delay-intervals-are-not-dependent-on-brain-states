%  Li Yuan, UCSD, 02-May-2023

close all
% Read in input information
sessInfo = SessInfoImport('V:\LiYuan\Codes\Fig8MazeTreadmill_V2\Fig8Treadmill_OnOff.xlsx');
sessDirs = sessInfo(1).sessDirs;
sessDirs2 = {'on10','off10','on30','off30'};

%%
i = 7;
% load spikes
% load average firing rate file
SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
load(SpikeShapeFile);
% load assembly file
assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_DelayOnOff_LR-25ms.mat');
load(assemblyFile);
Fig8ArmChoice = fullfile(sessInfo(i).mainDir,'Cell Property','Assembly_ArmChoice_DelayLR_25ms.mat');
load(Fig8ArmChoice);

clusterNum = length(SpikeProp.max_AvgRate);
rateCluster = SpikeProp.max_AvgRate;
rateLabel = rateCluster<5 & rateCluster>0.1;
tList = SpikeProp.tList(rateLabel);
% binTime = CellAssembly_DelayLR.binWidth;
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
AssmblWght_onff = CellAssembly_DelayLR.DelayOff.AssmblWght;
AssmblStrength_off = CellAssembly_DelayLR.DelayOff.AssmblStrength;
event_Time_off = CellAssembly_DelayLR.DelayOff.event_Time;
event_strength_off = CellAssembly_DelayLR.DelayOff.event_strength;
event_Num_off = CellAssembly_DelayLR.DelayOff.event_Num;
    
   for sesGroup = 1:4

            sesNumSum = 0;
            posLabel = [];
            turnLabel = [];

            behaveType = strsplit(sessDirs{sesGroup},'_');
            behaveType = behaveType{1};
            
            delayLR = Assembly_ArmChoice_DelayLR.(sessDirs2{sesGroup}).off.delayRateLRDiffAll;
            delayLR_95 = Assembly_ArmChoice_DelayLR.(sessDirs2{sesGroup}).off.shf_delayRateLRDiffAll95;
            stemLR = Assembly_ArmChoice_DelayLR.(sessDirs2{sesGroup}).off.stemRateLRDiffAll;
            stemLR_95 = Assembly_ArmChoice_DelayLR.(sessDirs2{sesGroup}).off.shf_stemRateLRDiffAll95;
            
            for j = [sesGroup,sesGroup+4]
                
                k = 1;
                % get each spike time
%                 binTime = event_Time_off{k};
                strength = AssmblStrength_off(k,:);
                 
                % load maps for each cluster
                pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'PathZone.mat');
                load(pathZoneFile);
                delayFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
                load(delayFile);
                
                % def1: delay starts at barrier
                % def2: delay starts at entrance
                delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
                delayTend1 = Fig8DelayZonePos.delayPos1.endT;
                %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
                %         delayTend2 = Fig8DelayZonePos.delayPos2.endT;
                
                trialNum = size(delayTstart1,2);
                if strcmp(sessDirs{j}(end-3:end-2),'10')
                    maxT = 10;
                    delayTend1_2 = delayTstart1+maxT;
                elseif strcmp(sessDirs{j}(end-3:end-2),'30')
                    maxT = 30;
                    delayTend1_2 = delayTstart1+maxT;
                else
                    error('Delay time is wrong')
                end
                
                        
                map_1D = load(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'));
                
                for m = 1:trialNum
                    if map_1D.ratesByECLR.ECLR(m) == 1
                        posLabel(sesNumSum+m) = 2;
                        turnLabel(sesNumSum+m) = 1;
                    elseif map_1D.ratesByECLR.ECLR(m) == 2
                        posLabel(sesNumSum+m) = 1;
                        turnLabel(sesNumSum+m) = 1;
                    elseif map_1D.ratesByECLR.ECLR(m) == 3
                        posLabel(sesNumSum+m) = 1;
                        turnLabel(sesNumSum+m) = 0;
                    elseif map_1D.ratesByECLR.ECLR(m) == 4
                        posLabel(sesNumSum+m) = 2;
                        turnLabel(sesNumSum+m) = 0;
                    else
                        error('Turning label ERROR');
                    end
                    
                    spikeInd = (binTime>=delayTstart1(m) & binTime<delayTend1_2(m));
                    
                    figure(1)
                    if turnLabel(sesNumSum+m) == 1
                        subplot(4,4,4*sesGroup-3)
                        plot((1:sum(spikeInd))*0.025,strength(spikeInd));
                        hold on
                        title('Delay L')
                        ylim([0 200])
                        %                             ylim([0 10])
                    else
                        subplot(4,4,4*sesGroup-2)
                        plot((1:sum(spikeInd))*0.025,strength(spikeInd));
                        hold on
                        title('Delay R')
                        ylim([0 200])
                        %                             ylim([0 10])
                    end
                    % stem
%                     figure(2)
                    spikeInd = (binTime>=(Fig8DelayZonePos.delayPos1.endT(m)+1/30) & binTime<PathZone.posEndT.Center(m));
                    if turnLabel(sesNumSum+m) == 1
                        subplot(4,4,4*sesGroup-1)
                        plot((1:sum(spikeInd))*0.025,strength(spikeInd));
                        hold on
                        title('Stem L')
                        ylim([0 100])
                        %                             ylim([0 10])
                    else
                        subplot(4,4,4*sesGroup)
                        plot((1:sum(spikeInd))*0.025,strength(spikeInd));
                        hold on
                        title('Stem R')
                        ylim([0 100])
                        %                             ylim([0 10])
                    end                   
                end
                sesNumSum = sesNumSum + sum(map_1D.ratesByECLR.valid);
            end
            
            subplot(4,4,4*sesGroup-2)
            TEXT1 = sprintf('%s%1.2f','L-R = ',delayLR(k));
            TEXT2 = sprintf('%s%1.2f','95% = ',delayLR_95(k));
            text(5,20,TEXT1);
            text(5,10,TEXT2);
            
            subplot(4,4,4*sesGroup)
            TEXT1 = sprintf('%s%1.2f','L-R = ',stemLR(k));
            TEXT2 = sprintf('%s%1.2f','95% = ',stemLR_95(k));
            text(1,20,TEXT1);
            text(1,10,TEXT2);
   end         
