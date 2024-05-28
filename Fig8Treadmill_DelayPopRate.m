function Fig8Treadmill_DelayPopRate(inFile,AnalyzeSes)

%  Li Yuan, UCSD, 11-May-2022
%  Calculate firing rate of each arm
%  for each assembly
% -------------------------------------------------------------------------
% set parameters for analysis

p.savePlot = 0;
p.writeToFile = 1;
% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'
p.timeBinWidth = 500./10^3; % unit sec

% Read in input information
sessInfo = SessInfoImport(inFile);

sessDirs2 = {'on10','off10','on30','off30'};

for i = AnalyzeSes(1:end)
    close all
    % selectivity for the whole map
    Fig8TreadmillDelayRatePop.session = sessInfo(i).mainDir;
    Fig8TreadmillDelayRatePop.rat = sessInfo(i).animal;
    Fig8TreadmillDelayRatePop.day = sessInfo(i).day;
    
    
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
    rateLabel = rateCluster<5 & rateCluster>0.1;    
%     tList = SpikeProp.tList(rateLabel);
    Fig8TreadmillDelayRatePop.clusterNum = clusterNum;
    
    if prod(~isempty(char(sessDirs{1})))
        for sesGroup = 1:4            
            if strcmp(sessDirs2{sesGroup}(end-1:end),'10')
                maxT = 10;
                binCount = maxT/p.timeBinWidth;
            elseif strcmp(sessDirs2{sesGroup}(end-1:end),'30')
                maxT = 30;
                binCount = maxT/p.timeBinWidth;
            else
                error('Delay time is wrong')
            end
                
            spike_all = zeros(clusterNum,binCount,0);
            
            behaveType = strsplit(sessDirs{sesGroup},'_');
            behaveType = behaveType{1};
            
            for j = [sesGroup,sesGroup+4]
                if ~strcmp(sessDirs2{sesGroup},sessDirs{j}(1:end-2))
                    error('session name error')
                end
             % load maps for each cluster;
                delayFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
                load(delayFile);
                
                % def1: delay starts at barrier
                % def2: delay starts at entrance
                delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
                delayTend1 = Fig8DelayZonePos.delayPos1.endT;
                delayTend1_2 = delayTstart1+maxT;
                %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
                %         delayTend2 = Fig8DelayZonePos.delayPos2.endT;
                
                trialNum = size(delayTstart1,2);
                
                % initiate the spike count temp
                spike_Temp = zeros(clusterNum,binCount,trialNum);
                
                for m = 1:trialNum                                   
                    rateBin1 = delayTstart1(m):p.timeBinWidth:delayTend1_2(m)-p.timeBinWidth;
                    rateBin2 = delayTstart1(m)+p.timeBinWidth:p.timeBinWidth:delayTend1_2(m);
                    for k = 1:clusterNum
                        if rateLabel(k) == 1
                            tSpTemp = Spike_Session.(sessDirs{j}){k};
                            fireTemp = zeros(1,binCount);
                            for n = 1:binCount
                                fireTemp(n) = sum(tSpTemp>rateBin1(n) & tSpTemp<=rateBin2(n));
                            end
                            %             spkTrain(k,:) = zscore(fireTemp);
                            spike_Temp(k,:,m) = fireTemp;
                        end
                    end
                end
                spike_all = cat(3,spike_all,spike_Temp);
            end
            delayRate = (mean(spike_all,1))/p.timeBinWidth;
            delayRatio = (sum(spike_all>0,1))/sum(rateLabel);
            Fig8TreadmillDelayRatePop.(behaveType).delayRate  = squeeze(delayRate);
            Fig8TreadmillDelayRatePop.(behaveType).delayRatio  = squeeze(delayRatio);
        end

        
        if p.writeToFile
            save(fullfile(savedir2,'Fig8TreadmillDelayRatePop.mat'), 'Fig8TreadmillDelayRatePop');
        end
        clear Fig8TreadmillArmRate map_1D
        close all
        fprintf('Finished Arm analysis for session %d\n',i);
    end

end
