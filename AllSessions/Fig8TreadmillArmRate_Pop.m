function Fig8TreadmillArmRate_Pop(inFile,AnalyzeSes)

%  Li Yuan, UCSD, 11-May-2022
%  Calculate firing rate of each arm
%  for each assembly
% -------------------------------------------------------------------------
% set parameters for analysis

p.savePlot = 0;
p.writeToFile = 1;
% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'
% Read in input information
sessInfo = SessInfoImport(inFile);

sessDirs2 = {'on10','on30','off10','off30'};

for i = AnalyzeSes(1:end)
    close all
    % selectivity for the whole map
    Fig8TreadmillArmRate.session = sessInfo(i).mainDir;
    Fig8TreadmillArmRate.rat = sessInfo(i).animal;
    Fig8TreadmillArmRate.day = sessInfo(i).day;
    
    
    if p.writeToFile
        savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
        if ~exist(savedir2, 'dir')
            mkdir(savedir2);
        end
    end
    
    display(['About to do session ' sessInfo(i).mainDir]);
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;    
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);

    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);
    if length(unique(TList))~=clusterNum
        error('TTList file has repeated clusters')
    end
    
    % get time stamps for the assembly
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.max_AvgRate;
%     rateLabel = rateCluster<5 & rateCluster>0.1;    
%     tList = SpikeProp.tList(rateLabel);
    
    if prod(~isempty(char(sessDirs{1})))
        for sesGroup = 1:4
            sesNumSum = 0;
            posLabel = [];
            turnLabel = [];
            time_All = [];
            
            behaveType = strsplit(sessDirs{sesGroup},'_');
            behaveType = behaveType{1};
            
            for j = [sesGroup,sesGroup+4]
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
                
                if exist(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'), 'file')==2
                    map_1D = load(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'));
                    posLabel = [];
                    turnLabel = [];
                    for m = 1:sum(map_1D.ratesByECLR.valid)
                        if map_1D.ratesByECLR.ECLR(m) == 1
                            posLabel(m,1) = 2;
                            turnLabel(m) = 1;
                        elseif map_1D.ratesByECLR.ECLR(m) == 2
                            posLabel(m,1) = 1;
                            turnLabel(m) = 1;
                        elseif map_1D.ratesByECLR.ECLR(m) == 3
                            posLabel(m,1) = 1;
                            turnLabel(m) = 0;
                        elseif map_1D.ratesByECLR.ECLR(m) == 4
                            posLabel(m,1) = 2;
                            turnLabel(m) = 0;
                        else
                            error('Turning label ERROR');
                        end
                        % return
                        timeTemp = PathZone.posEndT.Return(m)-PathZone.posStartT.Return(m);
                        time_All(m,1) = timeTemp;
                        % base
                        timeTemp = PathZone.posEndT.Base(m)-PathZone.posStartT.Base(m);
                        time_All(m,2) = timeTemp;
                        % delay
                        timeTemp = delayTend1_2(m)-delayTstart1(m);
                        time_All(m,3) = timeTemp;
                        % stem
                        timeTemp = PathZone.posEndT.Center(m)-(Fig8DelayZonePos.delayPos1.endT(m)+1/30);
                        time_All(m,4) = timeTemp;
                        % choice
                        timeTemp = PathZone.posEndT.Choice(m)-PathZone.posStartT.Choice(m);
                        time_All(m,5) = timeTemp;
                        % reward
                        timeTemp = PathZone.posEndT.Reward(m)-PathZone.posStartT.Reward(m);
                        time_All(m,6) = timeTemp;
                        
                    end
                end
                
                % find Left / Right turn trials
                leftIdx = find(posLabel(:,1)==1); % idx on all trials
                rightIdx = find(posLabel(:,1)==2); % idx on all trials
                
                LcorrectIdx = find(turnLabel' == 1 & posLabel(:,1)==2); % idx on all trials
                LerrorIdx = find(turnLabel' ==0 & posLabel(:,1)==1); % idx on all trials
                RcorrectIdx = find(turnLabel' == 1 & posLabel(:,1)==1); % idx on all trials
                RerrorIdx = find(turnLabel' == 0 & posLabel(:,1)==2); % idx on all trials
                
                
                for k = 1:clusterNum
                    % initiate variables
                    % return, delay, stem
                    spike_All = nan(sum(map_1D.ratesByECLR.valid),6);
                    % load path
                    pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'PathZone.mat');
                    load(pathZoneFile);
                    delayFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
                    load(delayFile);
                    
                    % get each spike time
                    tSp = Spike_Session.(sessDirs{j}){k};
                    
                    for m = 1:sum(map_1D.ratesByECLR.valid)
                        % get spike count
                        % return
                        spikeInd = (tSp>=PathZone.posStartT.Return(m) & tSp<PathZone.posEndT.Return(m));
                        spike_All(m,1) = sum(spikeInd);
                        % base
                        spikeInd = (tSp>=PathZone.posStartT.Base(m) & tSp<PathZone.posEndT.Base(m));
                        spike_All(m,2) = sum(spikeInd);
                        % delay
                        spikeInd = (tSp>=delayTstart1(m) & tSp<delayTend1_2(m));
                        spike_All(m,3) = sum(spikeInd);
                        % stem
                        spikeInd = (tSp>=(Fig8DelayZonePos.delayPos1.endT(m)+1/30) & tSp<PathZone.posEndT.Center(m));
                        spike_All(m,4) = sum(spikeInd);
                        % choice
                        spikeInd = (tSp>=PathZone.posStartT.Choice(m) & tSp<PathZone.posEndT.Choice(m));
                        spike_All(m,5) = sum(spikeInd);
                        % reward
                        spikeInd = (tSp>=PathZone.posStartT.Reward(m) & tSp<PathZone.posEndT.Reward(m));
                        spike_All(m,6) = sum(spikeInd);
                    end
                    
                    
                    % firing rate in different arms
                    rateReturn_L = nansum(nansum(spike_All(leftIdx,1)))./sum(sum(time_All(leftIdx,1)));
                    rateReturn_R= nansum(nansum(spike_All(rightIdx,1)))./sum(sum(time_All(rightIdx,1)));
                    
                    rateReturn = nansum(nansum(spike_All(:,1)))./sum(sum(time_All(:,1)));
                    % base
                    rateBase = nansum(nansum(spike_All(:,2)))./sum(sum(time_All(:,2)));
                    % delay zone
                    rateDelay = nansum(nansum(spike_All(:,3)))./sum(sum(time_All(:,3)));
                    % stem zone
                    rateStem = nansum(nansum(spike_All(:,4)))./sum(sum(time_All(:,4)));
                    % choice zone
                    rateChoice = nansum(nansum(spike_All(:,5)))./sum(sum(time_All(:,5)));
                    % reward zone
                    rateReward = nansum(nansum(spike_All(:,6)))./sum(sum(time_All(:,6)));
                    
                    Fig8TreadmillArmRate.(sessDirs{j}).rateReturn(k,:) = [rateReturn,sum(sum(time_All(:,1))),nansum(nansum(spike_All(:,1)))];
                    Fig8TreadmillArmRate.(sessDirs{j}).rateReturn_L(k,:) = [rateReturn_L,sum(sum(time_All(leftIdx,1))),nansum(nansum(spike_All(leftIdx,1)))];
                    Fig8TreadmillArmRate.(sessDirs{j}).rateReturn_R(k,:) = [rateReturn_R,sum(sum(time_All(rightIdx,1))),nansum(nansum(spike_All(rightIdx,1)))];
                    Fig8TreadmillArmRate.(sessDirs{j}).rateBase(k,:) = [rateBase,sum(sum(time_All(:,2))),nansum(nansum(spike_All(:,2)))];
                    Fig8TreadmillArmRate.(sessDirs{j}).rateDelay(k,:) = [rateDelay,sum(sum(time_All(:,3))),nansum(nansum(spike_All(:,3)))];
                    Fig8TreadmillArmRate.(sessDirs{j}).rateStem(k,:) = [rateStem,sum(sum(time_All(:,4))),nansum(nansum(spike_All(:,4)))];
                    Fig8TreadmillArmRate.(sessDirs{j}).rateChoice(k,:) = [rateChoice,sum(sum(time_All(:,5))),nansum(nansum(spike_All(:,5)))];
                    Fig8TreadmillArmRate.(sessDirs{j}).rateReward(k,:) = [rateReward,sum(sum(time_All(:,6))),nansum(nansum(spike_All(:,6)))];

                end
            end
                
        end
        
        if length(sessDirs) == 8
            sesGroupTemp = 1:4;
        else
            sesGroupTemp = 1:2;
        end
        
        for sesGroup = sesGroupTemp
            behaveType = strsplit(sessDirs{sesGroup},'_');
            behaveType = behaveType{1};
            
            sess_1 = strcat(behaveType,'_1');
            sess_2 = strcat(behaveType,'_2');
            
            % firing rate in different arms
            % RETURN
            spkNumTemp = Fig8TreadmillArmRate.(sess_1).rateReturn(:,3) + Fig8TreadmillArmRate.(sess_2).rateReturn(:,3);
            timeTemp = Fig8TreadmillArmRate.(sess_1).rateReturn(:,2) + Fig8TreadmillArmRate.(sess_2).rateReturn(:,2);
            rateTemp = spkNumTemp./timeTemp;
            Fig8TreadmillArmRate.(behaveType).rateReturn = [rateTemp,timeTemp,spkNumTemp];
            
            % RETURN L
            spkNumTemp = Fig8TreadmillArmRate.(sess_1).rateReturn_L(:,3) + Fig8TreadmillArmRate.(sess_2).rateReturn_L(:,3);
            timeTemp = Fig8TreadmillArmRate.(sess_1).rateReturn_L(:,2) + Fig8TreadmillArmRate.(sess_2).rateReturn_L(:,2);
            rateTemp = spkNumTemp./timeTemp;
            Fig8TreadmillArmRate.(behaveType).rateReturn_L = [rateTemp,timeTemp,spkNumTemp];
            
            % RETURN R
            spkNumTemp = Fig8TreadmillArmRate.(sess_1).rateReturn_R(:,3) + Fig8TreadmillArmRate.(sess_2).rateReturn_R(:,3);
            timeTemp = Fig8TreadmillArmRate.(sess_1).rateReturn_R(:,2) + Fig8TreadmillArmRate.(sess_2).rateReturn_R(:,2);
            rateTemp = spkNumTemp./timeTemp;
            Fig8TreadmillArmRate.(behaveType).rateReturn_R = [rateTemp,timeTemp,spkNumTemp];
            
            % base
            spkNumTemp = Fig8TreadmillArmRate.(sess_1).rateBase(:,3) + Fig8TreadmillArmRate.(sess_2).rateBase(:,3);
            timeTemp = Fig8TreadmillArmRate.(sess_1).rateBase(:,2) + Fig8TreadmillArmRate.(sess_2).rateBase(:,2);
            rateTemp = spkNumTemp./timeTemp;
            Fig8TreadmillArmRate.(behaveType).rateBase = [rateTemp,timeTemp,spkNumTemp];
            
            % delay zone
            spkNumTemp = Fig8TreadmillArmRate.(sess_1).rateDelay(:,3) + Fig8TreadmillArmRate.(sess_2).rateDelay(:,3);
            timeTemp = Fig8TreadmillArmRate.(sess_1).rateDelay(:,2) + Fig8TreadmillArmRate.(sess_2).rateDelay(:,2);
            rateTemp = spkNumTemp./timeTemp;
            Fig8TreadmillArmRate.(behaveType).rateDelay = [rateTemp,timeTemp,spkNumTemp];
            
            % stem zone
            spkNumTemp = Fig8TreadmillArmRate.(sess_1).rateStem(:,3) + Fig8TreadmillArmRate.(sess_2).rateStem(:,3);
            timeTemp = Fig8TreadmillArmRate.(sess_1).rateStem(:,2) + Fig8TreadmillArmRate.(sess_2).rateStem(:,2);
            rateTemp = spkNumTemp./timeTemp;
            Fig8TreadmillArmRate.(behaveType).rateStem = [rateTemp,timeTemp,spkNumTemp];
            
            % choice zone
            spkNumTemp = Fig8TreadmillArmRate.(sess_1).rateChoice(:,3) + Fig8TreadmillArmRate.(sess_2).rateChoice(:,3);
            timeTemp = Fig8TreadmillArmRate.(sess_1).rateChoice(:,2) + Fig8TreadmillArmRate.(sess_2).rateChoice(:,2);
            rateTemp = spkNumTemp./timeTemp;
            Fig8TreadmillArmRate.(behaveType).rateChoice = [rateTemp,timeTemp,spkNumTemp];
            
            % reward zone
            spkNumTemp = Fig8TreadmillArmRate.(sess_1).rateReward(:,3) + Fig8TreadmillArmRate.(sess_2).rateReward(:,3);
            timeTemp = Fig8TreadmillArmRate.(sess_1).rateReward(:,2) + Fig8TreadmillArmRate.(sess_2).rateReward(:,2);
            rateTemp = spkNumTemp./timeTemp;
            Fig8TreadmillArmRate.(behaveType).rateReward = [rateTemp,timeTemp,spkNumTemp];
            
        end
        
        % sleep sessions
        sleepDirs = sessInfo(i).sleepDirs;
        posFile = fullfile(sessInfo(i).mainDir,'processedData','indataS.mat');
        load(posFile);
        for j =  1:length(sleepDirs)
            pos = indata(j);
            startT = pos.t(1);
            endT = pos.t(end);
            
            for k = 1:clusterNum
                % get spikes
                tSp = Spike_Session.(sleepDirs{j}){k};
                spikeInd = (tSp>=startT & tSp<endT);
                rateSleep = sum(spikeInd)/(endT-startT);
                Fig8TreadmillArmRate.(sleepDirs{j}).rateSleep(k,:) = [rateSleep,endT-startT,sum(spikeInd)];
            end
        end
    
        
        if p.writeToFile
            save(fullfile(savedir2,'Fig8TreadmillArmRate.mat'), 'Fig8TreadmillArmRate');
        end
        clear Fig8TreadmillArmRate map_1D
        close all
        fprintf('Finished Arm analysis for session %d\n',i);
    end

end
