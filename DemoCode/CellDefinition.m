% delay active cell, time cell, time selective cell, consistently on cell
% definition code
close all
clear
savedir_1 = fullfile(cd,'ExamplePlot\DelayFireStability');
mkdir(savedir_1);
savedir_2 = fullfile(cd,'ExamplePlot\TimeCell_Definition');
mkdir(savedir_2);
delayFile ='Fig8DelayTimeMap_2Session.mat';
load(delayFile);

shuffleTimes = 1000;
% sessName2 = {'on10_1','off10_1','on30_1','off30_1','on10_2','off10_2','on30_2','off30_2'};
sessName2 = {'on10','off10','on30','off30'};

clusterNum = length(Fig8DelayTimeMap_2Session.tList);


%% calculate delay fire stability and define time cells
for j = 1:length(sessName2)
    maxT = sessName2{j}(end-1:end);
    delayCorr_Def1_Trial = nan(1,clusterNum);
    delayCorr_Def1_Trial_Shuffle = nan(shuffleTimes,clusterNum);
    delayCorr_Def1_Trial_Shuffle95 = nan(1,clusterNum);
    delayCorr_Def1_Trial_Shuffle99 = nan(1,clusterNum);
    delayCorr_Def1_Trial_Shuffle_Label95 = zeros(1,clusterNum);
    delayCorr_Def1_Trial_Shuffle_Label99 = zeros(1,clusterNum);
    
    
    for k = 1:clusterNum
        % def1, delay starts at barrier
        % get map for 2 sessions combined
        spikeRate1_2session = Fig8DelayTimeMap_2Session.(sessName2{j}).spikeRate1_Smooth{k};
        trialNum = size(spikeRate1_2session,1);
        % calculate correlation in a session trial by trial (only for
        % trials with values)
        % in each session, at least 5 trials should be active
        trialSumRate_1 = sum(spikeRate1_2session~=0,2);
        validTrial_1 = find(trialSumRate_1~=0);
        % count valid trials in first session
        halfTrialNumber = floor(trialNum/2);
        validSum1 = sum(trialSumRate_1(1:halfTrialNumber)~=0)>=5;
        validSum2 = sum(trialSumRate_1(halfTrialNumber+1:end)~=0)>=5;
        validState_1 = validSum1+validSum2 >1;
        if ~validState_1
            delayCorr_Def1_Trial(k) = NaN;
        else
            corrTemp_1 = nan(nchoosek(length(validTrial_1),2),1);
            countTemp_1 = 0;
            for m = 1:length(validTrial_1)-1
                for mm = m+1:length(validTrial_1)
                    countTemp_1 = countTemp_1+1;
                    corrTemp_1(countTemp_1) = Xcorrelate(spikeRate1_2session(validTrial_1(m),:),spikeRate1_2session(validTrial_1(mm),:));
                end
            end
            delayCorr_Def1_Trial(k) = median(corrTemp_1,'omitnan');
        end
          
        %% shuffle each trial rate bins (move in one direction and loop) and then do delayCorr_Half_Rotate
        binCount = size(spikeRate1_2session,2);
        for n = 1:shuffleTimes
            % shf from 3rd bin to (max-3)th bin
            shf_1 = ceil(rand(trialNum,1)*(binCount-6))+3;
            spikeRate1_2session_Shuffle = nan(trialNum,binCount);
            for m = 1:trialNum
                if shf_1 < binCount
                    shfOrder_1 = [shf_1(m)+1:binCount,1:shf_1(m)];
                else
                    shfOrder_1 = [shf_1(m):binCount,1:shf_1(m)-1];
                end
                spikeRate1_2session_Shuffle(m,:) = spikeRate1_2session(m,shfOrder_1);
            end
            
            % trial by trial shuffled correlation
            % def 1
            if ~validState_1
                delayCorr_Def1_Trial_Shuffle(n,k) = NaN;
            else
                corrTemp_1 = nan(nchoosek(length(validTrial_1),2),1);
                countTemp_1 = 0;
                for m = 1:length(validTrial_1)-1
                    for mm = m+1:length(validTrial_1)
                        countTemp_1 = countTemp_1+1;
                        corrTemp_1(countTemp_1) = Xcorrelate(spikeRate1_2session_Shuffle(validTrial_1(m),:),spikeRate1_2session_Shuffle(validTrial_1(mm),:));
                    end
                end
                delayCorr_Def1_Trial_Shuffle(n,k) = median(corrTemp_1,'omitnan');
            end
        end
        
        % get 95% and 99%
        % def 1
        shuffleTemp_1 = delayCorr_Def1_Trial_Shuffle(:,k);
        tempSort_1 = sort(shuffleTemp_1(~isnan(shuffleTemp_1)));
        if length(tempSort_1) > shuffleTimes/10
            delayCorr_Def1_Trial_Shuffle95(k) = tempSort_1(ceil(length(tempSort_1)*0.95));
            delayCorr_Def1_Trial_Shuffle99(k) = tempSort_1(ceil(length(tempSort_1)*0.99));
            if delayCorr_Def1_Trial(k) > delayCorr_Def1_Trial_Shuffle95(k)
                delayCorr_Def1_Trial_Shuffle_Label95(k) = 1;
            end
            if delayCorr_Def1_Trial(k) > delayCorr_Def1_Trial_Shuffle99(k)
                delayCorr_Def1_Trial_Shuffle_Label99(k) = 1;
            end
        else
            delayCorr_Def1_Trial_Shuffle95(k) = NaN;
            delayCorr_Def1_Trial_Shuffle99(k) = NaN;
        end

        % plot heat map
        if j == 1
            h = figure(k);
            h.Position = [100 100 800 600];
        end
        
        % def 1
        figure(k)
        subplot(2,length(sessName2),j)
        imagesc([0,binCount],[1 size(spikeRate1_2session,1)],spikeRate1_2session/max(max(spikeRate1_2session)))
        %             hold on
        %             imagesc([0,binCount],size(spikeRate2_2session,1)+1.5,spikeRate2_Combined_Smooth/max(spikeRate2_Combined_Smooth))
        colormap(jet)
        xlabel('Time (Sec)')
        ylabel('Trials')
        xlim([0 binCount])
        xTick= [0 binCount];
        xTickLabel = [0 str2double(maxT)];
        %             axis off
        set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
        
        if j == 1
            TITLE1 = sprintf('%s%d%s','Cell-',k,'DelayStability');
            TITLE2 = 'RateMap_NormbyTrial_Onset at barrier';
        else
            TITLE1 = sessName2{j};
            TITLE2 = [];
        end
        title({TITLE1;TITLE2},'Interpreter','None')
        xlabel(sprintf('%s%1.2f','Stability: ',delayCorr_Def1_Trial(k)));
        % plot shuffle distribution
        
        subplot(2,length(sessName2),j+length(sessName2))
        if ~isnan(delayCorr_Def1_Trial_Shuffle95(k))
            Violin(delayCorr_Def1_Trial_Shuffle(:,k),1,'ShowData',false);
            hold on
            plot([0.7,1.3],ones(1,2).*delayCorr_Def1_Trial_Shuffle95(k),'r-')
            plot(1,delayCorr_Def1_Trial(k),'r*')
            hold off
        else
            text(0.1,0.5,'Comparison not exist','FontSize',8);
        end
        %             axis off
        ylabel('shuffle compare')
        xlabel('Time shuffle')
        
    end

    
    DelayFireStability_TrialbyTrial_2Session.(sessName2{j}).delayCorr_Def1_Trial = delayCorr_Def1_Trial;
    %         Fig8DelayTimeMap_2SessionStability_Shuffle.(sessName2{j}).delayCorr_Trial = delayCorr_Trial;
    DelayFireStability_TrialbyTrial_2Session.(sessName2{j}).delayCorr_Def1_Trial_Shuffle = delayCorr_Def1_Trial_Shuffle;
    DelayFireStability_TrialbyTrial_2Session.(sessName2{j}).delayCorr_Def1_Trial_Shuffle95 = delayCorr_Def1_Trial_Shuffle95;
    DelayFireStability_TrialbyTrial_2Session.(sessName2{j}).delayCorr_Def1_Trial_Shuffle99 = delayCorr_Def1_Trial_Shuffle99;
    DelayFireStability_TrialbyTrial_2Session.(sessName2{j}).delayCorr_Def1_Trial_Shuffle_Label95 = delayCorr_Def1_Trial_Shuffle_Label95;
    DelayFireStability_TrialbyTrial_2Session.(sessName2{j}).delayCorr_Def1_Trial_Shuffle_Label99 = delayCorr_Def1_Trial_Shuffle_Label99;
end
    
    for k = 1:clusterNum
        figure(k)
        figName = sprintf('%s%s%d%s',savedir_1,'\Cell-',k,'-DelayStability');
        print(figName,'-dpng','-r300');
    end
close all

%% define time cells
p.peak = 0; % unit Hz
p.sdThreshold = 2;
p.fieldTreshold = 0.2;
for j = 1:length(sessName2)
        
        shuffleSig_1 = DelayFireStability_TrialbyTrial_2Session.(sessName2{j}).delayCorr_Def1_Trial_Shuffle_Label95;
        timeField_1 = nan(1,clusterNum);
        fieldLabel_1 = cell(1,clusterNum);
        
        for k = 1:clusterNum
            
            % get map for 2 sessions combined
            % def 1
            spikeRate1_2session = Fig8DelayTimeMap_2Session.(sessName2{j}).spikeRate1_Smooth{k};
            mapMean_Def1_2Session = nanmean(spikeRate1_2session,1);
            avgRate1 = mean(mapMean_Def1_2Session);
            peakRate_1 = max(mapMean_Def1_2Session);
            std_Def1_2Session = std(mapMean_Def1_2Session);         

%             if shuffleSig_1(k) == 1 && (peakRate_1 >= (avgRate1 + p.sdThreshold*std_Def1_2Session))
            if shuffleSig_1(k) == 1
                % calculate information/spike, info/sec, coherence, sparseness
                [nFields_Def1,FieldBinX_Def1,fieldLabelTemp_Def1] = timeFieldSearch(mapMean_Def1_2Session,p);
            else
                nFields_Def1 = 0;
                fieldLabelTemp_Def1 = zeros(length(mapMean_Def1_2Session),1);
            end
            timeField_1(k) = nFields_Def1;
            fieldLabel_1{k} = fieldLabelTemp_Def1;
            
            
            if j == 1                
                h = figure(k);
                h.Position = [100 100 1600 600];
            end
            
            figure(k)
            subplot(1,length(sessName2),j)
            imagesc([0,size(spikeRate1_2session,2)],[20-size(spikeRate1_2session,1)+1 size(spikeRate1_2session,1)],spikeRate1_2session/max(max(spikeRate1_2session)))
            hold on
            imagesc([0,size(mapMean_Def1_2Session,2)],size(spikeRate1_2session,1)+2,mapMean_Def1_2Session/max(mapMean_Def1_2Session))
            colormap(jet)
            if nFields_Def1 ==1
                plot(FieldBinX_Def1-1.5,ones(1,length(FieldBinX_Def1))*23,'r','LineWidth',3)
            end
            
            if j == 1
                TITLE1 = sprintf('%s%d%s','Cell-',k,'TimeField');
            else
                TITLE1 = sessName2{j};
            end
            title({TITLE1},'Interpreter','None')
            axis off
            text2 = sprintf('%s%d','Shuffled significance: ',shuffleSig_1(k));
            text1 = sprintf('%s%2.2f%s','Peak rate : ',peakRate_1,' Hz');
            text3 = sprintf('%s%d','TimeField: ',nFields_Def1);
            text([1,1,1],[24,25,26],{text1,text2,text3},'FontSize',10,'Interpreter','None')
            ylim([-2 size(spikeRate1_2session,1)+10])  
            
        end 
        DelayField_TrialbyTrial_2Session.(sessName2{j}).timeField_1 = timeField_1;        
        DelayField_TrialbyTrial_2Session.(sessName2{j}).fieldLabel_1 = fieldLabel_1;             
        
end

    for k = 1:clusterNum
        figure(k)
        figName = sprintf('%s%s%d%s',savedir_2,'\Cell-',k,'-TimeField');
        print(figName,'-dpng','-r300');
    end
close all

%% plot all delay active cell, time cells etc
% load average firing rate file
p.avgRateThres = 0.5;
SpikeShapeFile = 'SpikeShape.mat';
load(SpikeShapeFile);
% load arm rate file
armRateFile = 'Fig8TreadmillArmRate.mat';
load(armRateFile);
idx = find(SpikeProp.AvgRate.Fig8Rate>0.1 & SpikeProp.AvgRate.Fig8Rate<5);
avgRate = SpikeProp.AvgRate.Fig8Rate;

    
for j = 1:length(sessName2)
    timeMap_Def1.(sessName2{j}) = [];
    pre_Map_Def1.(sessName2{j}) = [];
    post_Map_Def1.(sessName2{j}) = [];
    endField_1.(sessName2{j}) = [];
    
    fieldLabel_Def1.(sessName2{j}) = DelayField_TrialbyTrial_2Session.(sessName2{j}).timeField_1(idx);
    % average correlation value of the 2 blocks
    corr_Def1.(sessName2{j}) = DelayFireStability_TrialbyTrial_2Session.(sessName2{j}).delayCorr_Def1_Trial(idx);
    corrLabel_Def1.(sessName2{j}) = DelayFireStability_TrialbyTrial_2Session.(sessName2{j}).delayCorr_Def1_Trial_Shuffle_Label95(idx);
    
    %         blockFieldTemp = DelayField_TrialbyTrial_Shuffle.(sess1).timeField_1(idx) + DelayField_TrialbyTrial_Shuffle.(sess2).timeField_1(idx);
    %         block_field_Def1.(sessName2{j}) = [block_field_Def1.(sessName2{j}),blockFieldTemp];
    %         blockFieldTemp = DelayField_TrialbyTrial_Shuffle.(sess1).timeField_2(idx) + DelayField_TrialbyTrial_Shuffle.(sess2).timeField_2(idx);
    %         block_field_Def2.(sessName2{j}) = [block_field_Def2.(sessName2{j}),blockFieldTemp];
    
    for k = idx
        timeMap_Def1.(sessName2{j}) = [timeMap_Def1.(sessName2{j});Fig8DelayTimeMap_2Session.(sessName2{j}).spikeRate1_Combined_Smooth{k}];
        pre_Map_Def1.(sessName2{j}) = [pre_Map_Def1.(sessName2{j});Fig8DelayTimeMap_2Session.(sessName2{j}).pre_spikeRate1_Combined_Smooth{k}];
        post_Map_Def1.(sessName2{j}) = [post_Map_Def1.(sessName2{j});Fig8DelayTimeMap_2Session.(sessName2{j}).post_mapMean_Def1_2Session{k}(1:20)];
        
        if sum(DelayField_TrialbyTrial_2Session.(sessName2{j}).fieldLabel_1{k}(end-6:end))>0
            endField_1.(sessName2{j}) = [endField_1.(sessName2{j}),1];
        else
            endField_1.(sessName2{j}) = [endField_1.(sessName2{j}),0];
        end       
    end
end

rate_DelayTemp = [Fig8TreadmillArmRate.on10.rateDelay(idx,1)';Fig8TreadmillArmRate.off10.rateDelay(idx,1)';...
    Fig8TreadmillArmRate.on30.rateDelay(idx,1)';Fig8TreadmillArmRate.off30.rateDelay(idx,1)'];
delay_onIdx = (sum(rate_DelayTemp>p.avgRateThres,1))>0;

h2 = figure;
h2.Position = [100 100 1200 400];
% on10 off10
on10Idx_Time = (fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
on10_mean = mean(cellMapTemp1,2);

off10Idx_Time = (fieldLabel_Def1.off10>0 & corr_Def1.off10>0 & delay_onIdx);
cellMapTemp2 = timeMap_Def1.off10(off10Idx_Time,:);
off10_mean = mean(cellMapTemp2,2);

plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),...
    1,[1,4,1],[1,4,2],{'Delay start at barrier','Time cells on10'},{'-3','0','10','13'},{'Time cells off10'},{'-3','0','10','13'});

% on30 and off30
on30Idx_Time = (fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx);
cellMapTemp1 = timeMap_Def1.on30(on30Idx_Time,:);
on30_mean = mean(cellMapTemp1,2);

off30Idx_Time = (fieldLabel_Def1.off30>0 & corr_Def1.off30>0 & delay_onIdx);
cellMapTemp2 = timeMap_Def1.off30(off30Idx_Time,:);
off30_mean = mean(cellMapTemp2,2);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on30(on30Idx_Time,:),post_Map_Def1.on30(on30Idx_Time,:),pre_Map_Def1.off30(off30Idx_Time,:),post_Map_Def1.off30(off30Idx_Time,:),...
    1,[1,4,3],[1,4,4],{'Delay start at barrier','Time cells on30'},{'-3','0','30','33'},{'Time cells off30'},{'-3','0','30','33'});

% all time limited cells
h1 = figure;
h1.Position = [100 100 1200 400];
% on10 off10
on10Idx_Time = (fieldLabel_Def1.on10>0 & corr_Def1.on10>0 & delay_onIdx & ~endField_1.on10);
cellMapTemp1 = timeMap_Def1.on10(on10Idx_Time,:);
on10_mean = mean(cellMapTemp1,2);

off10Idx_Time = (fieldLabel_Def1.off10>0 & corr_Def1.off10>0  & delay_onIdx & ~endField_1.off10);
cellMapTemp2 = timeMap_Def1.off10(off10Idx_Time,:);
off10_mean = mean(cellMapTemp2,2);

plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx_Time,:),post_Map_Def1.on10(on10Idx_Time,:),pre_Map_Def1.off10(off10Idx_Time,:),post_Map_Def1.off10(off10Idx_Time,:),...
    1,[1,4,1],[1,4,2],{'Delay start at barrier','Time-limited cells on10'},{'-3','0','10','13'},{'Time-limited cells off10'},{'-3','0','10','13'});

% on30 and off30
on30Idx_Time = (fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & ~endField_1.on30);
cellMapTemp1 = timeMap_Def1.on30(on30Idx_Time,:);
on30_mean = mean(cellMapTemp1,2);

off30Idx_Time = (fieldLabel_Def1.off30>0 & corr_Def1.off30>0  & delay_onIdx & ~endField_1.off30);
cellMapTemp2 = timeMap_Def1.off30(off30Idx_Time,:);
off30_mean = mean(cellMapTemp2,2);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on30(on30Idx_Time,:),post_Map_Def1.on30(on30Idx_Time,:),pre_Map_Def1.off30(off30Idx_Time,:),post_Map_Def1.off30(off30Idx_Time,:),...
    1,[1,4,3],[1,4,4],{'Delay start at barrier','Time-limited cells on30'},{'-3','0','30','33'},{'Time-limited cells off30'},{'-3','0','30','33'});

% consistent on cells
h1 = figure;
h1.Position = [100 100 1200 400];
% on10 off10
on10Idx = (fieldLabel_Def1.on10>0 & corr_Def1.on10>0  & delay_onIdx & endField_1.on10);
cellMapTemp1 = timeMap_Def1.on10(on10Idx,:);

off10Idx = (fieldLabel_Def1.off10>0 & corr_Def1.off10>0  & delay_onIdx & endField_1.off10);
cellMapTemp2 = timeMap_Def1.off10(off10Idx,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx,:),post_Map_Def1.on10(on10Idx,:),pre_Map_Def1.off10(off10Idx,:),post_Map_Def1.off10(off10Idx,:),...
    1,[1,4,1],[1,4,2],{'Delay start at barrier','consistent on on10'},{'-3','0','10','13'},{'consistent on off10'},{'-3','0','10','13'});


% on30 and off30
on30Idx = (fieldLabel_Def1.on30>0 & corr_Def1.on30>0 & delay_onIdx & endField_1.on30);
cellMapTemp1 = timeMap_Def1.on30(on30Idx,:);

off30Idx = (fieldLabel_Def1.off30>0 & corr_Def1.off30>0  & delay_onIdx & endField_1.off30);
cellMapTemp2 = timeMap_Def1.off30(off30Idx,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on30(on30Idx,:),post_Map_Def1.on30(on30Idx,:),pre_Map_Def1.off30(off30Idx,:),post_Map_Def1.off30(off30Idx,:),...
    1,[1,4,3],[1,4,4],{'Delay start at barrier','consistent on on30'},{'-3','0','30','33'},{'consistent on off30'},{'-3','0','30','33'});


% delayActive cells
h1 = figure;
h1.Position = [100 100 1200 400];
% on10 off10
on10Idx = (delay_onIdx == 1);
cellMapTemp1 = timeMap_Def1.on10(on10Idx,:);

off10Idx = (delay_onIdx == 1);
cellMapTemp2 = timeMap_Def1.off10(off10Idx,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on10(on10Idx,:),post_Map_Def1.on10(on10Idx,:),pre_Map_Def1.off10(off10Idx,:),post_Map_Def1.off10(off10Idx,:),...
    1,[1,4,1],[1,4,2],{'Delay start at barrier','Delay on on10'},{'-3','0','10','13'},{'Delay on off10'},{'-3','0','10','13'});


% on30 and off30
on30Idx = (delay_onIdx == 1);
cellMapTemp1 = timeMap_Def1.on30(on30Idx,:);

off30Idx = (delay_onIdx == 1);
cellMapTemp2 = timeMap_Def1.off30(off30Idx,:);
plotTimeMap2_Pre_Post(cellMapTemp1,cellMapTemp2,pre_Map_Def1.on30(on30Idx,:),post_Map_Def1.on30(on30Idx,:),pre_Map_Def1.off30(off30Idx,:),post_Map_Def1.off30(off30Idx,:),...
    1,[1,4,3],[1,4,4],{'Delay start at barrier','Delay on on30'},{'-3','0','30','33'},{'Delay on off30'},{'-3','0','30','33'});



% Xcorrelate calculates the pearson correlation between the rate map from two
% cells. Pixels that are visited less than 150 ms in the maps are also omitted
% from both cells.
function [corrValue,pix] = Xcorrelate(map1,map2)
% Both the rate maps and the time maps will have equal size
[M,N] = size(map1);
[M2,N2] = size(map2);

if M~=M2 || N~=N2
    error('Correlation sizes do not match!');
else
% % Exclude bins that was visited less than 150 ms in either room
% for ii = 1:N
%     for jj = 1:M
%         if timeMap1(ii,jj) < 0.150 | timeMap2(ii,jj) < 0.150
%             map1(ii,jj) = NaN;
%             map2(ii,jj) = NaN;
%         end
%     end
% end

% Transform the 2-D maps to 1-D arrays by assembling the columns from the
% maps
    A = reshape(map1,M*N,1);
    B = reshape(map2,M*N,1);
    % Find the pixels containing NaN in A
    An = isnan(A);
    index = find(An==0);
    % Remove the bins with NaN in A from both maps
    A = A(index);
    B = B(index);
    % Find the pixels that still contain NaN in B
    Bn = isnan(B);
    index = find(Bn==0);
    % Remove the bins with NaN in B from both maps
    A = A(index);
    B = B(index);
    % Number of bins used in the correlation
    pix = length(A);

    % Calculate the correlation for the two rate maps.
    if length(A) ~= length(B)
        error('Correlation sizes do not match!');
    elseif isempty(A) || isempty(B)
            corrValue = NaN; % Cannot do calculation on empty arrays
        else
            corrC = corrcoef(A,B,'Rows','pairwise');
            corrValue = corrC(1,2);
    end
end
end

% this is to detect single field
function [nFields,FieldBinX,fieldLabel] = timeFieldSearch(map,p)

% Counter for the number of fields
FieldBinX = [];

% Allocate memory to the arrays
M = length(map);
% Array that contain the bins of the map this algorithm has visited
visited = zeros(M,1);
fieldLabel = visited;

nanInd = isnan(map);
visited(nanInd) = 1;

% Array that will contain the bin positions to the current placefield

% Find the current maximum
[peak,maxId] = max(map);
visited(map<p.fieldTreshold*peak) = 1;
% sd = std(map,'omitnan');

% if peak >= p.peak && peak>= nanmean(map) + p.sdThreshold*sd
if peak >= p.peak
    % look for time field
    
    % Find the bins that construct the peak field
    [binsX,visited] = recursiveBins_Lin(map,visited,[],maxId,M);
    binsX = sort(binsX);
    if length(binsX)>= 1 % Minimum size of a placefield
        nFields = 1;
        % Total rate
        visited(binsX) = 1; % modified by Li
        % Bins in Field
        FieldBinX = binsX;
        fieldLabel(binsX) = 1;
    end
else
    nFields = 0;
end
end


function [binsX,visited] = recursiveBins_Lin(map,visited,binsX,ii,M)
% If outside boundaries of map -> return.
if ii<1 || ii>M
    return;
end
% If all bins are visited -> return.
if prod(prod(visited))
    return;
end
if visited(ii) % This bin has been visited before
    return;
else
    binsX = [binsX;ii];
    visited(ii) = 1;
    % Call this function again in each of the 24 neighbour bins
    [binsX,visited] = recursiveBins_Lin(map,visited,binsX,ii+1,M);
    [binsX,visited] = recursiveBins_Lin(map,visited,binsX,ii-1,M);
end
end

% sortMethod = 1, sort by map its own
% sort method = 2, map2 sort by map1
function plotTimeMap2_Pre_Post(map1,map2,pre_map1,post_map1,pre_map2,post_map2,sortMethod,pos1,pos2,TITLEText1,TICKText1,TITLEText2,TICKText2)
if sortMethod == 1
    [~,peakIdx] = max(map1,[],2);
    [~,peak_Sort] = sort(peakIdx);
    cellMapTemp1_Sort = map1(peak_Sort,:);   
    pre_map1_sort = pre_map1(peak_Sort,:);
    post_map1_sort = post_map1(peak_Sort,:);
    
    [~,peakIdx] = max(map2,[],2);
    [~,peak_Sort] = sort(peakIdx);
    cellMapTemp2_Sort = map2(peak_Sort,:);
    pre_map2_sort = pre_map2(peak_Sort,:);
    post_map2_sort = post_map2(peak_Sort,:);

else
    [~,peakIdx] = max(map1,[],2);
    [~,peak_Sort] = sort(peakIdx);
    cellMapTemp1_Sort = map1(peak_Sort,:);
    cellMapTemp2_Sort = map2(peak_Sort,:);
    pre_map1_sort = pre_map1(peak_Sort,:);
    post_map1_sort = post_map1(peak_Sort,:);
    pre_map2_sort = pre_map2(peak_Sort,:);
    post_map2_sort = post_map2(peak_Sort,:);
end

% cellMapTemp1_Sort_Norm = cellMapTemp1_Sort;
% cellMapTemp2_Sort_Norm = cellMapTemp2_Sort;

cellMapTemp1_Sort_Norm = [pre_map1_sort,cellMapTemp1_Sort,post_map1_sort]./max(cellMapTemp1_Sort,[],2);
cellMapTemp2_Sort_Norm = [pre_map2_sort,cellMapTemp2_Sort,post_map2_sort]./max(cellMapTemp2_Sort,[],2);

subplot(pos1(1),pos1(2),pos1(3))
imagesc(cellMapTemp1_Sort_Norm)
colormap(jet)
title(TITLEText1,'Interpreter','None')
axis on
set(gca, 'xtick', [0,size(pre_map1_sort,2),size([pre_map1_sort,cellMapTemp1_Sort],2),size(cellMapTemp1_Sort_Norm,2)]);
set(gca, 'xticklabels', TICKText1);
clim([0 1])

subplot(pos2(1),pos2(2),pos2(3))
imagesc(cellMapTemp2_Sort_Norm)
colormap(jet)
title(TITLEText2,'Interpreter','None')
axis on
set(gca, 'xtick', [0,size(pre_map2_sort,2),size([pre_map2_sort,cellMapTemp2_Sort],2),size(cellMapTemp2_Sort_Norm,2)]);
set(gca, 'xticklabels', TICKText2);
clim([0 1])

end

        