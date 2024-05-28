function Fig8TreadmillCrossCorr(inFile,AnalyzeSes)
close all

p.savePlot = 1;
p.writeToFile = 1;

% set crossCorr time bin and time length
p.binSize = 100./1000; % unit sec
p.xcorr_width_msec = 5000./1000; % unit: sec

% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'

% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\Cell pair crosscorr');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end
    
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);

    for k = 1:clusterNum
        DelayCrossCorr.tList{k} = TList{k}(1:end-2);
    end
    tList = DelayCrossCorr.tList;
    
    % initiate the data
    DelayCrossCorr.rat = sessInfo(i).animal;
    DelayCrossCorr.day = sessInfo(i).day;
    DelayCrossCorr.timeBin = p.binSize;
              
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    
    % get each phase names (no delay etc)
    SessDirs = sessInfo(i).sessDirs;
    

    for j = 1:length(SessDirs)
        
        % load analyzed positions
        delayFile = fullfile(mainDir,SessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
%         pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
%         pathData = load(pathFile);
%         
        % -----------------------------------------------------------------
        % use the mode of delay is always 10 sec or 30 sec
        % rather than identify real time in delay zone
        % to make plot cleaner
        % I can change to delay time in delay zone later
        % Li Yuan, 19-Aug-2020, UCSD
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
        delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
        delayTend2 = Fig8DelayZonePos.delayPos2.endT;
        
        trialNum = size(delayTstart1,2);
        if contains(SessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
            delayTend2_2 = delayTstart2+maxT;
        elseif contains(SessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
            delayTend2_2 = delayTstart2+maxT;
        else
            error('Delay time is wrong')
        end       
                
        % initiate timestamps
        ts_Delay1_Combine = cell(clusterNum,1);
        ts_Delay1_Combine_wholeDelay = cell(clusterNum,1);
%         ts_Delay1_Trial = cell(clusterNum,trialNum);
%         ts_Delay2_Combine = cell(clusterNum,1);
%         ts_Delay2_Trial = cell(clusterNum,trialNum);
        
        ts_NonDelay1_Combine = cell(clusterNum,1);
        ts_NonDelay2_Combine = cell(clusterNum,1);
        % -----------------------------------------------------------------
        for k = 1:clusterNum           
            % get each spike time
            tSp = Spike_Session.(SessDirs{j}){k};
            if ~isempty(tSp)
                for m = 1:trialNum
                    
                    % delay definition 1
                    ts_Delay1 = tSp(tSp>delayTstart1(m) & tSp<delayTend1_2(m));
                    %                 ts_Delay1 = ts_Delay1-delayTstart1(m);
                    %                 % each trial ts
                    %                 ts_Delay1_Trial(k,m) = ts_Delay1;
                    % add trials into combined ts
                    ts_Delay1_CombineTemp = ts_Delay1;
                    ts_Delay1_Combine{k} = [ts_Delay1_Combine{k};ts_Delay1_CombineTemp];
                    
                    ts_Delay1 = tSp(tSp>delayTstart1(m) & tSp<delayTend1(m));
                    ts_Delay1_Combine_wholeDelay{k} = [ts_Delay1_Combine_wholeDelay{k};ts_Delay1];
                    
%                     % delay definition 2
%                     ts_Delay2 = tSp(tSp>delayTstart2(m) & tSp<delayTend2_2(m));
%                     %                 ts_Delay2 = ts_Delay2-delayTstart2(m);
%                     %                 % each trial ts
%                     %                 ts_Delay2_Trial(k,m) = ts_Delay2;
%                     % add trials into combined ts
%                     ts_Delay2_CombineTemp = ts_Delay2;
%                     ts_Delay2_Combine{k} = [ts_Delay2_Combine{k};ts_Delay2_CombineTemp];
                    
                end
                % no delay spikes
                ts_nonDelayIdx = ~ismember(tSp,ts_Delay1_Combine_wholeDelay{k});
                ts_NonDelay1_Combine{k} = [ts_NonDelay1_Combine{k};tSp(ts_nonDelayIdx)];
%                 % no delay spikes
%                 ts_nonDelayIdx = ~ismember(tSp,ts_Delay2_Combine{k});
%                 ts_NonDelay2_Combine{k} = [ts_NonDelay2_Combine{k};tSp(ts_nonDelayIdx)];
            else
                % get all as 0x1 double
                ts_Delay1_Combine{k} =  double.empty(0,1);
%                 ts_Delay2_Combine{k} = double.empty(0,1);
                ts_NonDelay1_Combine{k} = double.empty(0,1);
%                 ts_NonDelay2_Combine{k} = double.empty(0,1);
            end
        end
        
        % calculate cofire and cross correlation for each pair
        if clusterNum > 1
            cellPairCount = 0;
            for k = 1:clusterNum-1
                for m = k+1:clusterNum
                    cellPairCount = cellPairCount + 1;
                    
                    % calculate Cross-Corr of delay area spikes
                    % use CrossCorr function from mclust 4.4
                    nbins = ceil(p.xcorr_width_msec/p.binSize)*2+1;
                    if sum(ts_Delay1_Combine{k})~=0 && sum(ts_Delay1_Combine{m})~=0
                        [Y_delay1,xdim] = ccf(ts_Delay1_Combine{k}, ts_Delay1_Combine{m}, p.binSize, p.xcorr_width_msec);
                    else
                        Y_delay1 = nan(nbins,1);
                        xdim = -1*p.xcorr_width_msec:p.binSize:p.xcorr_width_msec;
                    end
%                     if sum(ts_Delay2_Combine{k})~=0 && sum(ts_Delay2_Combine{m})~=0
%                         [Y_delay2,xdim] = ccf(ts_Delay2_Combine{k}, ts_Delay2_Combine{m}, p.binSize, p.xcorr_width_msec);
%                     else
%                         Y_delay2 = nan(nbins,1);
%                         xdim = -1*p.xcorr_width_msec:p.binSize:p.xcorr_width_msec;
%                     end
                    if sum(ts_NonDelay1_Combine{k})~=0 && sum(ts_NonDelay1_Combine{m})~=0
                        [Y_nodelay,xdim] = ccf(ts_NonDelay1_Combine{k}, ts_NonDelay1_Combine{m}, p.binSize, p.xcorr_width_msec);
                    else
                        Y_nodelay = nan(nbins,1);
                        xdim = -1*p.xcorr_width_msec:p.binSize:p.xcorr_width_msec;
                    end
                    
                    DelayCrossCorr.(SessDirs{j}).CellPairIDX1(cellPairCount) = k;
                    DelayCrossCorr.(SessDirs{j}).CellPairIDX2(cellPairCount) = m;
                    DelayCrossCorr.(SessDirs{j}).CellPair{cellPairCount} = {tList{k} tList{m}};
                    DelayCrossCorr.(SessDirs{j}).Delay1{cellPairCount} = Y_delay1;
%                     DelayCrossCorr.(SessDirs{j}).Delay2{cellPairCount} = Y_delay2;
                    DelayCrossCorr.(SessDirs{j}).NoDelay{cellPairCount} = Y_nodelay;
                    
                end
            end
        end     
    end
    
    pairNum = length(DelayCrossCorr.(SessDirs{j}).CellPairIDX1);
    for n = 1:pairNum
        for j=1:length(SessDirs)
            h=figure(n);
            if j == 1
                h.Position = [50 50 1800 900];
                TITLE1 =  sprintf('%s%d%s%d%s','Theta lock -',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',SessDirs{j});
            else
                TITLE1 = SessDirs{j};
            end
            TITLE2 = ['Delay area XCorr ' DelayCrossCorr.(SessDirs{j}).CellPair{n}{1} ' Vs ' DelayCrossCorr.(SessDirs{j}).CellPair{n}{2}];
            TITLE3 = ['Non-Delay area XCorr ' DelayCrossCorr.(SessDirs{j}).CellPair{n}{1} ' Vs ' DelayCrossCorr.(SessDirs{j}).CellPair{n}{2}];
            hh(j) = subplot(4,4,j);
            plot(xdim,DelayCrossCorr.(SessDirs{j}).Delay1{n},'Color',[1,0.2,0.2])
            xlabel('sec')
            ylabel('Count')
            title({TITLE1;TITLE2},'Interpreter','None');
            hhh(j) = subplot(4,4,j+8);
            plot(xdim,DelayCrossCorr.(SessDirs{j}).NoDelay{n},'Color',[0,0,1])
            xlabel('sec')
            ylabel('Count')
            title({TITLE1;TITLE3},'Interpreter','None');
        end
        linkaxes(hh,'y');
        linkaxes(hhh,'y');
        
        if p.savePlot == 1           
            figName = sprintf('%s%s%d%s%d%s%s%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',...
                DelayCrossCorr.(SessDirs{j}).CellPair{n}{1},'vs',DelayCrossCorr.(SessDirs{j}).CellPair{n}{2},'-PairwiseCrossCorr');
            print(figName,'-dpng','-r300');         
        end
        close(h)
    end
    
    if p.writeToFile == 1
        save(fullfile(savedir2,'DelayCrossCorr.mat'), 'DelayCrossCorr');
    end
    clear DelayCrossCorr
    fprintf('Finished analysis for session %d\n',i);
end

end

