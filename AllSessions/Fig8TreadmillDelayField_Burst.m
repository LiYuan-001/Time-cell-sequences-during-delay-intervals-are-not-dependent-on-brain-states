% time field in this code is firstly defined by whether the within session
% shuffled stability passed 5 percentile
% stability is definsed as median value of all trial by trial correlation
% in one session
% then extract field by peak rate
function Fig8TreadmillDelayField_Burst(inFile,AnalyzeSes)
close all

p.savePlot = 1;
p.writeToFile = 1;

p.peak = 0; % unit Hz
p.sdThreshold = 2;
p.fieldTreshold = 0.2;
% p.timeGap = 1; % if field is seperated by 1 bin, concatenate time field
    
% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\DelayTimeField_Burst');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end
    
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    % load DelayFire_Burst map
    delayFile = fullfile(savedir2, 'Fig8DelayTimeMap_Burst.mat');
    load(delayFile);
    % load correlation values
    stabilityFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayFireStability_Burst_2Session.mat');
    load(stabilityFile);

    % initiate the data
    DelayField_Burst_2Session.rat = sessInfo(i).animal;
    DelayField_Burst_2Session.day = sessInfo(i).day;
    DelayField_Burst_2Session.timeBin = DelayFire_Burst.timeBin;
    DelayField_Burst_2Session.gaussSigma = DelayFire_Burst.gaussSigma;

    Fig8DelayTimeMap_2Session_Burst.rat = sessInfo(i).animal;
    Fig8DelayTimeMap_2Session_Burst.day = sessInfo(i).day;
    Fig8DelayTimeMap_2Session_Burst.timeBin = DelayFire_Burst.timeBin;
    Fig8DelayTimeMap_2Session_Burst.gaussSigma = DelayFire_Burst.gaussSigma;
     
    TList = DelayFire_Burst.tList;
    clusterNum = length(DelayFire_Burst.tList);
    Fig8DelayTimeMap_2Session_Burst.TList = TList;
    
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    map_Burst_1 = [];
    map_Burst_2 = [];
    map_SingleSpk_1 = [];
    map_SingleSpk_2 = [];
    for j = 1:length(sessDirs)
        sessName = sessDirs{j}(1:end-2);
        trialNumTemp = size(DelayFire_Burst.(sessDirs{j}).delayTstart1,2);
        if contains(sessDirs{j},'_1')
            trialNum_2Session.(sessName)(1) = trialNumTemp;
            trialInd = [1:trialNum_2Session.(sessName)(1)];
        elseif contains(sessDirs{j},'_2')
            trialNum_2Session.(sessName)(2) = trialNumTemp;
            trialInd = [trialNum_2Session.(sessName)(1)+1:trialNum_2Session.(sessName)(1)+trialNumTemp];
        end
        binCount = size(DelayFire_Burst.(sessDirs{j}).Burst.spikeRate1_Smooth{1},2);
        prebinCount = size(DelayFire_Burst.(sessDirs{j}).Burst.pre_spikeRate1_Smooth{1},2);
        % -----------------------------------------------------------------
        for k = 1:clusterNum
            spikeRate1_Smooth = DelayFire_Burst.(sessDirs{j}).Burst.spikeRate1_Smooth{k};
            map_Burst_1.(sessName)(trialInd,1:binCount,k) = spikeRate1_Smooth;
            spikeRate2_Smooth = DelayFire_Burst.(sessDirs{j}).Burst.spikeRate2_Smooth{k};
            map_Burst_2.(sessName)(trialInd,1:binCount,k) = spikeRate2_Smooth;
            
            spikeRate1_Smooth = DelayFire_Burst.(sessDirs{j}).SingleSpk.spikeRate1_Smooth{k};
            map_SingleSpk_1.(sessName)(trialInd,1:binCount,k) = spikeRate1_Smooth;
            spikeRate2_Smooth = DelayFire_Burst.(sessDirs{j}).SingleSpk.spikeRate2_Smooth{k};
            map_SingleSpk_2.(sessName)(trialInd,1:binCount,k) = spikeRate2_Smooth;
%             premap_1.(sessName)(trialInd,1:prebinCount,k) = DelayFire_Burst.(sessDirs{j}).pre_spikeRate1_Smooth{k};
%             premap_2.(sessName)(trialInd,1:prebinCount,k) = DelayFire_Burst.(sessDirs{j}).pre_spikeRate2_Smooth{k};
        end
    end
    clear trialNum_2Session
    
    if length(sessDirs) == 8
        sessName2 = {'on10','off10','on30','off30'};
    elseif contains(sessDirs{1},'on')
        sessName2 = {'on10','on30'};
    else
        sessName2 = {'off10','off30'};
    end
    
    for j = 1:length(sessName2)
        
        %% burst
        shuffleSig_1 = DelayFireStability_Burst_2Session.(sessName2{j}).Burst.delayCorr_Def1_Trial_Shuffle_Label99;
        timeField_1 = nan(1,clusterNum);
        fieldLabel_1 = cell(1,clusterNum);
        
        shuffleSig_2 = DelayFireStability_Burst_2Session.(sessName2{j}).Burst.delayCorr_Def2_Trial_Shuffle_Label99;
        timeField_2 = nan(1,clusterNum);
        fieldLabel_2 = cell(1,clusterNum);
        
        for k = 1:clusterNum
            
            % get map for 2 sessions combined
            % def 1
            spikeRate1_2session = map_Burst_1.(sessName2{j})(:,:,k);
            mapMean_Def1_2Session = nanmean(spikeRate1_2session,1);
            avgRate1 = mean(mapMean_Def1_2Session);
            peakRate_1 = max(mapMean_Def1_2Session);
            std_Def1_2Session = std(mapMean_Def1_2Session);
            
            % def 2
            spikeRate2_2session = map_Burst_2.(sessName2{j})(:,:,k);
            mapMean_Def2_2Session = nanmean(spikeRate2_2session,1);
            avgRate2 = mean(mapMean_Def2_2Session);
            peakRate_2 = max(mapMean_Def2_2Session);
            std_Def2_2Session = std(mapMean_Def2_2Session);
            
%             % get pre map for 2 sessions combined
%             % def 1
%             pre_spikeRate1_2session = premap_1.(sessName2{j})(:,:,k);
%             pre_mapMean_Def1_2Session = nanmean(pre_spikeRate1_2session,1);
%             
%             % def 2
%             pre_spikeRate2_2session = premap_2.(sessName2{j})(:,:,k);
%             pre_mapMean_Def2_2Session = nanmean(pre_spikeRate2_2session,1);
            
            % put map into the .mat file
            Fig8DelayTimeMap_2Session_Burst.(sessName2{j}).Burst.spikeRate1_Smooth{k} = spikeRate1_2session;
            Fig8DelayTimeMap_2Session_Burst.(sessName2{j}).Burst.spikeRate1_Combined_Smooth{k} = mapMean_Def1_2Session;
            
            Fig8DelayTimeMap_2Session_Burst.(sessName2{j}).Burst.spikeRate2_Smooth{k} = spikeRate2_2session;
            Fig8DelayTimeMap_2Session_Burst.(sessName2{j}).Burst.spikeRate2_Combined_Smooth{k} = mapMean_Def2_2Session;
            
%             Fig8DelayTimeMap_2Session_Assembly.(sessName2{j}).pre_spikeRate1_Smooth{k} = pre_spikeRate1_2session;
%             Fig8DelayTimeMap_2Session_Assembly.(sessName2{j}).pre_spikeRate1_Combined_Smooth{k} = pre_mapMean_Def1_2Session;
%             
%             Fig8DelayTimeMap_2Session_Assembly.(sessName2{j}).pre_spikeRate2_Smooth{k} = pre_spikeRate2_2session;
%             Fig8DelayTimeMap_2Session_Assembly.(sessName2{j}).pre_spikeRate2_Combined_Smooth{k} = pre_mapMean_Def2_2Session;
            
            
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
            
%             if shuffleSig_2(k) == 1 && (peakRate_2 >= (avgRate2 + p.sdThreshold*std_Def2_2Session))
                % calculate information/spike, info/sec, coherence, sparseness
            if shuffleSig_2(k) == 1
                [nFields_Def2,FieldBinX_Def2,fieldLabelTemp_Def2] = timeFieldSearch(mapMean_Def2_2Session,p);
            else
                nFields_Def2 = 0;
                fieldLabelTemp_Def2 = zeros(length(mapMean_Def2_2Session),1);
            end
            timeField_2(k) = nFields_Def2;
            fieldLabel_2{k} = fieldLabelTemp_Def2;
            
            if j == 1                
                h = figure(k);
                h.Position = [100 100 1600 900];
                h = figure(k+clusterNum);
                h.Position = [100 100 1600 900];
            end
            
            figure(k)
            subplot(2,length(sessName2),j)
            imagesc([0,size(spikeRate1_2session,2)],[20-size(spikeRate1_2session,1)+1 size(spikeRate1_2session,1)],spikeRate1_2session/max(max(spikeRate1_2session)))
            hold on
            imagesc([0,size(mapMean_Def1_2Session,2)],22,mapMean_Def1_2Session/max(mapMean_Def1_2Session))
            colormap(jet)
            if nFields_Def1 ==1
                plot(FieldBinX_Def1-1.5,ones(1,length(FieldBinX_Def1))*23,'r','LineWidth',3)
            end
            
            if j == 1
                TITLE1 = sprintf('%s%d%s%d%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-Burst-',sessName2{j});
            else
                TITLE1 = sessName2{j};
            end
            title({TITLE1},'Interpreter','None')
            axis off
            text2 = sprintf('%s%d','Shuffled significance: ',shuffleSig_1(k));
            text1 = sprintf('%s%2.2f%s','Peak rate : ',peakRate_1,' Hz');
            text3 = sprintf('%s%d','TimeField: ',nFields_Def1);
            text([1,1,1],[24,25,26],{text1,text2,text3},'FontSize',10,'Interpreter','None')
            ylim([-2 30])
            
            
            figure(k+clusterNum)
            subplot(2,length(sessName2),j)
            imagesc([0,size(spikeRate2_2session,2)],[20-size(spikeRate2_2session,1)+1 size(spikeRate2_2session,1)],spikeRate2_2session/max(max(spikeRate2_2session)))
            hold on
            imagesc([0,size(mapMean_Def2_2Session,2)],22,mapMean_Def2_2Session/max(mapMean_Def2_2Session))
            colormap(jet)
            if nFields_Def2 ==1
                plot(FieldBinX_Def2-1.5,ones(1,length(FieldBinX_Def2))*23,'r','LineWidth',3)
            end
            
            if j == 1
                TITLE1 = sprintf('%s%d%s%d%s%s','DelayZone2-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-Burst-',sessName2{j});
            else
                TITLE1 = sessName2{j};
            end
            title({TITLE1},'Interpreter','None')
            axis off
            text2 = sprintf('%s%d','Shuffled significance: ',shuffleSig_2(k));
            text1 = sprintf('%s%2.2f%s','Peak rate : ',peakRate_2,' Hz');
            text3 = sprintf('%s%d','TimeField: ',nFields_Def2);
            text([1,1,1],[24,25,26],{text1,text2,text3},'FontSize',10,'Interpreter','None')
            ylim([-2 30])
            
        end 
        DelayField_Burst_2Session.(sessName2{j}).Burst.timeField_1 = timeField_1;        
        DelayField_Burst_2Session.(sessName2{j}).Burst.fieldLabel_1 = fieldLabel_1;     
        DelayField_Burst_2Session.(sessName2{j}).Burst.timeField_2 = timeField_2;        
        DelayField_Burst_2Session.(sessName2{j}).Burst.fieldLabel_2 = fieldLabel_2;
        
        %% single spike
        shuffleSig_1 = DelayFireStability_Burst_2Session.(sessName2{j}).SingleSpk.delayCorr_Def1_Trial_Shuffle_Label99;
        timeField_1 = nan(1,clusterNum);
        fieldLabel_1 = cell(1,clusterNum);
        
        shuffleSig_2 = DelayFireStability_Burst_2Session.(sessName2{j}).SingleSpk.delayCorr_Def2_Trial_Shuffle_Label99;
        timeField_2 = nan(1,clusterNum);
        fieldLabel_2 = cell(1,clusterNum);
        
        for k = 1:clusterNum
            
            % get map for 2 sessions combined
            % def 1
            spikeRate1_2session = map_SingleSpk_1.(sessName2{j})(:,:,k);
            mapMean_Def1_2Session = nanmean(spikeRate1_2session,1);
            avgRate1 = mean(mapMean_Def1_2Session);
            peakRate_1 = max(mapMean_Def1_2Session);
            std_Def1_2Session = std(mapMean_Def1_2Session);
            
            % def 2
            spikeRate2_2session = map_SingleSpk_2.(sessName2{j})(:,:,k);
            mapMean_Def2_2Session = nanmean(spikeRate2_2session,1);
            avgRate2 = mean(mapMean_Def2_2Session);
            peakRate_2 = max(mapMean_Def2_2Session);
            std_Def2_2Session = std(mapMean_Def2_2Session);
            
%             % get pre map for 2 sessions combined
%             % def 1
%             pre_spikeRate1_2session = premap_1.(sessName2{j})(:,:,k);
%             pre_mapMean_Def1_2Session = nanmean(pre_spikeRate1_2session,1);
%             
%             % def 2
%             pre_spikeRate2_2session = premap_2.(sessName2{j})(:,:,k);
%             pre_mapMean_Def2_2Session = nanmean(pre_spikeRate2_2session,1);
            
            % put map into the .mat file
            Fig8DelayTimeMap_2Session_Burst.(sessName2{j}).SingleSpk.spikeRate1_Smooth{k} = spikeRate1_2session;
            Fig8DelayTimeMap_2Session_Burst.(sessName2{j}).SingleSpk.spikeRate1_Combined_Smooth{k} = mapMean_Def1_2Session;
            
            Fig8DelayTimeMap_2Session_Burst.(sessName2{j}).SingleSpk.spikeRate2_Smooth{k} = spikeRate2_2session;
            Fig8DelayTimeMap_2Session_Burst.(sessName2{j}).SingleSpk.spikeRate2_Combined_Smooth{k} = mapMean_Def2_2Session;
            
%             Fig8DelayTimeMap_2Session_Assembly.(sessName2{j}).pre_spikeRate1_Smooth{k} = pre_spikeRate1_2session;
%             Fig8DelayTimeMap_2Session_Assembly.(sessName2{j}).pre_spikeRate1_Combined_Smooth{k} = pre_mapMean_Def1_2Session;
%             
%             Fig8DelayTimeMap_2Session_Assembly.(sessName2{j}).pre_spikeRate2_Smooth{k} = pre_spikeRate2_2session;
%             Fig8DelayTimeMap_2Session_Assembly.(sessName2{j}).pre_spikeRate2_Combined_Smooth{k} = pre_mapMean_Def2_2Session;
            
            
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
            
%             if shuffleSig_2(k) == 1 && (peakRate_2 >= (avgRate2 + p.sdThreshold*std_Def2_2Session))
                % calculate information/spike, info/sec, coherence, sparseness
            if shuffleSig_2(k) == 1
                [nFields_Def2,FieldBinX_Def2,fieldLabelTemp_Def2] = timeFieldSearch(mapMean_Def2_2Session,p);
            else
                nFields_Def2 = 0;
                fieldLabelTemp_Def2 = zeros(length(mapMean_Def2_2Session),1);
            end
            timeField_2(k) = nFields_Def2;
            fieldLabel_2{k} = fieldLabelTemp_Def2;
            
%             if j == 1                
%                 h = figure(k);
%                 h.Position = [100 100 1600 600];
%                 h = figure(k+clusterNum);
%                 h.Position = [100 100 1600 600];
%             end
            
            figure(k)
            subplot(2,length(sessName2),j+4)
            imagesc([0,size(spikeRate1_2session,2)],[20-size(spikeRate1_2session,1)+1 size(spikeRate1_2session,1)],spikeRate1_2session/max(max(spikeRate1_2session)))
            hold on
            imagesc([0,size(mapMean_Def1_2Session,2)],22,mapMean_Def1_2Session/max(mapMean_Def1_2Session))
            colormap(jet)
            if nFields_Def1 ==1
                plot(FieldBinX_Def1-1.5,ones(1,length(FieldBinX_Def1))*23,'r','LineWidth',3)
            end
            
            if j == 1
                TITLE1 = sprintf('%s%d%s%d%s%s','DelayZone1-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-SingleSpk-',sessName2{j});
            else
                TITLE1 = sessName2{j};
            end
            title({TITLE1},'Interpreter','None')
            axis off
            text2 = sprintf('%s%d','Shuffled significance: ',shuffleSig_1(k));
            text1 = sprintf('%s%2.2f%s','Peak rate : ',peakRate_1,' Hz');
            text3 = sprintf('%s%d','TimeField: ',nFields_Def1);
            text([1,1,1],[24,25,26],{text1,text2,text3},'FontSize',10,'Interpreter','None')
            ylim([-2 30])
            
            
            figure(k+clusterNum)
            subplot(2,length(sessName2),j+4)
            imagesc([0,size(spikeRate2_2session,2)],[20-size(spikeRate2_2session,1)+1 size(spikeRate2_2session,1)],spikeRate2_2session/max(max(spikeRate2_2session)))
            hold on
            imagesc([0,size(mapMean_Def2_2Session,2)],22,mapMean_Def2_2Session/max(mapMean_Def2_2Session))
            colormap(jet)
            if nFields_Def2 ==1
                plot(FieldBinX_Def2-1.5,ones(1,length(FieldBinX_Def2))*23,'r','LineWidth',3)
            end
            
            if j == 1
                TITLE1 = sprintf('%s%d%s%d%s%s','DelayZone2-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-SingleSpk-',sessName2{j});
            else
                TITLE1 = sessName2{j};
            end
            title({TITLE1},'Interpreter','None')
            axis off
            text2 = sprintf('%s%d','Shuffled significance: ',shuffleSig_2(k));
            text1 = sprintf('%s%2.2f%s','Peak rate : ',peakRate_2,' Hz');
            text3 = sprintf('%s%d','TimeField: ',nFields_Def2);
            text([1,1,1],[24,25,26],{text1,text2,text3},'FontSize',10,'Interpreter','None')
            ylim([-2 30])
            
        end 
        DelayField_Burst_2Session.(sessName2{j}).SingleSpk.timeField_1 = timeField_1;        
        DelayField_Burst_2Session.(sessName2{j}).SingleSpk.fieldLabel_1 = fieldLabel_1;     
        DelayField_Burst_2Session.(sessName2{j}).SingleSpk.timeField_2 = timeField_2;        
        DelayField_Burst_2Session.(sessName2{j}).SingleSpk.fieldLabel_2 = fieldLabel_2;
    end
    
    if p.savePlot == 1
        for k = 1:clusterNum
            figure(k)
            figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',TList{k},'-DelayTimeField_Burst_Onset-barrier');
            print(figName,'-dpng','-r300');
            figure(k+clusterNum)
            figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',TList{k},'-DelayTimeField__Burst_Onset-entrance');
            print(figName,'-dpng','-r300');
        end
    end
    close all
    if p.writeToFile == 1
        save(fullfile(savedir2,'DelayField_Burst_2Session.mat'), 'DelayField_Burst_2Session');
        save(fullfile(savedir2,'Fig8DelayTimeMap_2Session_Burst.mat'), 'Fig8DelayTimeMap_2Session_Burst');
    end
    clear DelayField_Assembly_2Session Fig8DelayTimeMap_2Session_Assembly map_1 map_2
    fprintf('Finished analysis for session %d\n',i);
end

end
