% it needs to be stable in each session
% and stable in 2 session combined
%
function Fig8TreadmillTimeCellQuant_Ratio_Session(inFile,AnalyzeSes)

close all

p.savePlot = 1;
p.writeToFile = 1;
p.avgRateThres = 0.5;

% % Read in input information
sessInfo = SessInfoImport(inFile);

if p.savePlot
    % directory for plot figures
    % generate a folder for each rat eah day under the current folder
    savedir = sprintf('%s%s',cd,'\Figures\Quantification Figures\Time cell num-dist quant');
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end
end

% get each phase names (no delay etc)
if length(sessInfo(1).sessDirs) == 8
    sessDirs = {'on10','off10','on30','off30'};
elseif contains(sessInfo(1).sessDirs{1},'on')
    sessDirs = {'on10','on30'};
else
    sessDirs = {'off10','off30'};
end

for j = 1:length(sessDirs)
    delay_on.(sessDirs{j}) = [];
    timeCell.(sessDirs{j}) = [];
    timeLimited.(sessDirs{j}) = [];
    persistent_on.(sessDirs{j}) = [];
    peak5.(sessDirs{j}) = [];
    cellNum.(sessDirs{j}) = [];
    
    delay_on_animal.(sessDirs{j}) = [];
    timeCell_animal.(sessDirs{j}) = [];
    timeLimited_animal.(sessDirs{j}) = [];
    persistent_on_animal.(sessDirs{j}) = [];
    peak5_animal.(sessDirs{j}) = [];
    cellNum_animal.(sessDirs{j}) = [];
    
    delay_on_all.(sessDirs{j}) = [];
    timeCell_all.(sessDirs{j}) = [];
    timeLimited_all.(sessDirs{j}) = [];
    persistent_on_all.(sessDirs{j}) = [];
    peak5_all.(sessDirs{j}) = [];
    cellNum_all.(sessDirs{j}) = [];
end
  
ratID =[];
for i = AnalyzeSes(1:end)
    
    avgRate = [];    
    delay_onIdx = [];
    rate_Delay = [];
    
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load time field file
    timeFieldFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayField_TrialbyTrial_2Session.mat');
    load(timeFieldFile);
    % load delayFire map
    delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap_2Session.mat');
    load(delayFile);
    % load correlation values
    stabilityFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayFireStability_TrialbyTrial_2Session.mat');
    load(stabilityFile);
    % load arm rate file
    armRateFile = fullfile(sessInfo(i).mainDir,'Cell Property', 'Fig8TreadmillArmRate.mat');
    load(armRateFile);    
    
    clusterNum = length(SpikeProp.AvgRate.Fig8Rate);
    idx = find(SpikeProp.AvgRate.Fig8Rate>0.1 & SpikeProp.AvgRate.Fig8Rate<5);
    avgRate = [avgRate,SpikeProp.AvgRate.Fig8Rate];
        
    ratID = [ratID;sessInfo(i).animal];
    for j = 1:length(sessDirs)
        fieldLabel_Def1.(sessDirs{j}) = [];
        fieldLabel_Def2.(sessDirs{j}) = [];
        timeMap_Def1.(sessDirs{j}) = [];
        timeMap_Def2.(sessDirs{j}) = [];
        pre_Map_Def1.(sessDirs{j}) = [];
        post_Map_Def1.(sessDirs{j}) = [];
        corr_Def1.(sessDirs{j}) = [];
        corr_Def2.(sessDirs{j}) = [];
        endField_1.(sessDirs{j}) = [];
        endField_2.(sessDirs{j}) = [];
        corrLabel_Def1.(sessDirs{j}) = [];
        corrLabel_Def2.(sessDirs{j}) = [];
        peakIdx.(sessDirs{j}) = [];
        
        fieldLabel_Def1.(sessDirs{j}) = [fieldLabel_Def1.(sessDirs{j}),DelayField_TrialbyTrial_2Session.(sessDirs{j}).timeField_1(idx)];
        fieldLabel_Def2.(sessDirs{j}) = [fieldLabel_Def2.(sessDirs{j}),DelayField_TrialbyTrial_2Session.(sessDirs{j}).timeField_2(idx)];
        
        % average correlation value of the 2 blocks
        corr_Def1.(sessDirs{j}) = [corr_Def1.(sessDirs{j}),DelayFireStability_TrialbyTrial_2Session.(sessDirs{j}).delayCorr_Def1_Trial(idx)]; 
        corr_Def2.(sessDirs{j}) = [corr_Def2.(sessDirs{j}),DelayFireStability_TrialbyTrial_2Session.(sessDirs{j}).delayCorr_Def2_Trial(idx)]; 
        
        corrLabel_Def1.(sessDirs{j}) = [corrLabel_Def1.(sessDirs{j}),DelayFireStability_TrialbyTrial_2Session.(sessDirs{j}).delayCorr_Def1_Trial_Shuffle_Label95(idx)];
        corrLabel_Def2.(sessDirs{j}) = [corrLabel_Def2.(sessDirs{j}),DelayFireStability_TrialbyTrial_2Session.(sessDirs{j}).delayCorr_Def2_Trial_Shuffle_Label95(idx)];
        
        sess1 = strcat(sessDirs{j},'_1');
        sess2 = strcat(sessDirs{j},'_2');

%         blockFieldTemp = DelayField_TrialbyTrial_Shuffle.(sess1).timeField_1(idx) + DelayField_TrialbyTrial_Shuffle.(sess2).timeField_1(idx);
%         block_field_Def1.(sessDirs{j}) = [block_field_Def1.(sessDirs{j}),blockFieldTemp];
%         blockFieldTemp = DelayField_TrialbyTrial_Shuffle.(sess1).timeField_2(idx) + DelayField_TrialbyTrial_Shuffle.(sess2).timeField_2(idx);
%         block_field_Def2.(sessDirs{j}) = [block_field_Def2.(sessDirs{j}),blockFieldTemp];
        
        for k = idx
            timeMap_Def1.(sessDirs{j}) = [timeMap_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate1_Combined_Smooth{k}];
            timeMap_Def2.(sessDirs{j}) = [timeMap_Def2.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate2_Combined_Smooth{k}];
 
            pre_Map_Def1.(sessDirs{j}) = [pre_Map_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).pre_spikeRate1_Combined_Smooth{k}];
            post_Map_Def1.(sessDirs{j}) = [post_Map_Def1.(sessDirs{j});Fig8DelayTimeMap_2Session.(sessDirs{j}).post_mapMean_Def1_2Session{k}(1:20)];
            
            if sum(DelayField_TrialbyTrial_2Session.(sessDirs{j}).fieldLabel_1{k}(end-6:end))>0
                endField_1.(sessDirs{j}) = [endField_1.(sessDirs{j}),1];
            else
                endField_1.(sessDirs{j}) = [endField_1.(sessDirs{j}),0];
            end
            
            if sum(DelayField_TrialbyTrial_2Session.(sessDirs{j}).fieldLabel_2{k}(end-6:end))>0
                endField_2.(sessDirs{j}) = [endField_2.(sessDirs{j}),1];
            else
                endField_2.(sessDirs{j}) = [endField_2.(sessDirs{j}),0];
            end
            
        end 
        [~,peakInd] = max(timeMap_Def1.(sessDirs{j}),[],2);
        peakIdx.(sessDirs{j}) = peakInd';
        rate_Delay = [rate_Delay;Fig8TreadmillArmRate.(sessDirs{j}).rateDelay(idx,1)'];  
    end 
 
    delay_onIdx = (sum(rate_Delay>p.avgRateThres,1))>0;    
    
    for j = 1:length(sessDirs)
        delay_onIdx_temp = sum(rate_Delay(j,:)>p.avgRateThres,1)>0;
        
        delay_on.(sessDirs{j}) = [delay_on.(sessDirs{j});sum(delay_onIdx_temp)];
        
        timeCell.(sessDirs{j}) = [timeCell.(sessDirs{j});sum(fieldLabel_Def1.(sessDirs{j})>0 & ...
            corr_Def1.(sessDirs{j})>0 & delay_onIdx)];
        
        timeLimited.(sessDirs{j}) = [timeLimited.(sessDirs{j});sum(fieldLabel_Def1.(sessDirs{j})>0 & ...
            corr_Def1.(sessDirs{j})>0 & delay_onIdx & ~endField_1.(sessDirs{j}))];
        
        persistent_on.(sessDirs{j}) = [persistent_on.(sessDirs{j});sum(fieldLabel_Def1.(sessDirs{j})>0 & ...
            corr_Def1.(sessDirs{j})>0 & delay_onIdx & endField_1.(sessDirs{j}))];
        
        peak5.(sessDirs{j}) = [peak5.(sessDirs{j});sum(fieldLabel_Def1.(sessDirs{j})>0 & ...
            corr_Def1.(sessDirs{j})>0 & delay_onIdx & ~endField_1.(sessDirs{j}) & peakIdx.(sessDirs{j})<=ceil(5/0.15))];
        
        cellNum.(sessDirs{j}) = [cellNum.(sessDirs{j});length(delay_onIdx_temp)];
    end 
end

close all
