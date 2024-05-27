% Read in sequences and form input for checkmotseq.m
% run after code function checkmotseq_SeqReadin(inFile,AnalyzeSes) and rt=rtfun2_LiModified(narr)
% Nov-19-2021, Li Yuan, UCSD
function checkTemp_TimeCell_Batch(inFile,AnalyzeSes)


p.writeToFile = 1;
p.cellNum = 16;
p.shuffle = 100;

% Read in input information
sessInfo = SessInfoImport(inFile);

rtFile = 'W:\Li Yuan\Codes\Fig8MazeTreadmill\Cell Property\ReactivationCode-Leibold\matlab code\Results\rt.mat';
load(rtFile);
  
zmat_TimeCell.LL = [];
zmat_TimeCell.LR = [];
zmat_TimeCell.RR = [];
zmat_TimeCell.RL = [];

zmat_NonTimeCell.LL = [];
zmat_NonTimeCell.LR = [];
zmat_NonTimeCell.RR = [];
zmat_NonTimeCell.RL = [];

for i = AnalyzeSes(1:end)
    close all

    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    % load reactivation file
    seqFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_RegionSeq.mat');
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load time field file
    timeFieldFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayTimeField_Trial_Shuffle.mat');
    load(timeFieldFile);
    
    % get valid cell ind
    % rate < 5 Hz in whole session
    % time cell / non-time cell
    clusterNum = length(SpikeProp.max_AvgRate);
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5;  
    rateFig8 = SpikeProp.AvgRate.Fig8Rate;
    rateFig8Label = rateFig8>0.1 & rateFig8<5;
    
    % DelayField_TrialbyTrial_Shuffle.(SessDirs{j}).timeField_1(idx);
    % label time cell
    timeCellLabel = (DelayField_TrialbyTrial_Shuffle.on10_1.timeField_1+DelayField_TrialbyTrial_Shuffle.on10_2.timeField_1+...
        DelayField_TrialbyTrial_Shuffle.off10_1.timeField_1+DelayField_TrialbyTrial_Shuffle.off10_2.timeField_1+...
        DelayField_TrialbyTrial_Shuffle.on30_1.timeField_1+DelayField_TrialbyTrial_Shuffle.on30_2.timeField_1+...
        DelayField_TrialbyTrial_Shuffle.off30_1.timeField_1+DelayField_TrialbyTrial_Shuffle.off30_2.timeField_1)>0;
    timeCellLabel = timeCellLabel&rateFig8Label;
    nonTimeCell = ~timeCellLabel&rateLabel;
    
    if sum(rateLabel) >= p.cellNum
         
        load(seqFile);
        
%         DelayPopFire_RegionSeq.rat = sessInfo(i).animal;
%         DelayPopFire_RegionSeq.day = sessInfo(i).day;
%         DelayPopFire_RegionSeq.clusterNum = clusterNum;
%         DelayPopFire_RegionSeq.PyrCellNum = sum(rateLabel);
%              
%         on_eventTsp2Start.delay 
%         on_eventTsp2Start.reward 
%                   
%         off_eventTsp2Start.delay
%         off_eventTsp2Start.reward
% 
%         sleep_eventTsp2Start 
          
        % calculate template then test template
        
        leftInd = DelayPopFire_RegionSeq.off_eventTsp2Start.delayLeftInd;
        rightInd = ~DelayPopFire_RegionSeq.off_eventTsp2Start.delayLeftInd;
        firstSpkTime_Left = DelayPopFire_RegionSeq.off_eventTsAllspks.delay(leftInd==1,:);
        firstSpkTime_Right = DelayPopFire_RegionSeq.off_eventTsAllspks.delay(rightInd==1,:);
        
        firstSpkTime_Left_TimeCell = firstSpkTime_Left(:,timeCellLabel);
        firstSpkTime_Right_TimeCell = firstSpkTime_Right(:,timeCellLabel);        

        firstSpkTime_Left_Non_TimeCell = firstSpkTime_Left(:,nonTimeCell);
        firstSpkTime_Right_Non_TimeCell = firstSpkTime_Right(:,nonTimeCell);  
               
        firstSpkTime_Left_2 = [];
        firstSpkTime_Right_2 = [];

        seq_Left_TimeCell = [];
        mm = 0;
        
        for k = 1:size(firstSpkTime_Left_TimeCell,1)
            [tsTemp,seqTemp] = sort(firstSpkTime_Left_TimeCell(k,:));
            seqTemp = seqTemp(~isnan(tsTemp));
            if length(seqTemp)>=3
                mm = mm+1;
                seq_Left_TimeCell{mm,1} = seqTemp;
                firstSpkTime_Left_2(mm,:) = firstSpkTime_Left_TimeCell(k,:);
            end
        end

        seq_Right_TimeCell = [];
        mm = 0;
        for k = 1:size(firstSpkTime_Right_TimeCell,1)
            [tsTemp,seqTemp] = sort(firstSpkTime_Right_TimeCell(k,:));
            seqTemp = seqTemp(~isnan(tsTemp));
            if length(seqTemp)>=3
                mm = mm+1;
                seq_Right_TimeCell{mm,1} = seqTemp;
                firstSpkTime_Right_2(mm,:) = firstSpkTime_Right_TimeCell(k,:);
            end
        end
        
        seqNum_Left_TimeCell = length(seq_Left_TimeCell);
        seqNum_Right_TimeCell = length(seq_Right_TimeCell);
        
        [~,template_Left_TimeCell] = sort(mean(firstSpkTime_Left_2,1,'omitnan'));
        [~,template_Right_TimeCell] = sort(mean(firstSpkTime_Right_2,1,'omitnan'));

        % non-time cell
        firstSpkTime_Left_3 = [];
        firstSpkTime_Right_3 = [];

        seq_Left_Non_TimeCell = [];
        mm = 0;
        for k = 1:size(firstSpkTime_Left_Non_TimeCell,1)
            [tsTemp,seqTemp] = sort(firstSpkTime_Left_Non_TimeCell(k,:));
            seqTemp = seqTemp(~isnan(tsTemp));
            if length(seqTemp)>=3
                mm = mm+1;
                seq_Left_Non_TimeCell{mm,1} = seqTemp;
                firstSpkTime_Left_3(mm,:) = firstSpkTime_Left_Non_TimeCell(k,:);
            end
        end

        seq_Right_Non_TimeCell = [];
        mm = 0;
        for k = 1:size(firstSpkTime_Right_Non_TimeCell,1)
            [tsTemp,seqTemp] = sort(firstSpkTime_Right_Non_TimeCell(k,:));
            seqTemp = seqTemp(~isnan(tsTemp));
            if length(seqTemp)>=3
                mm = mm+1;
                seq_Right_Non_TimeCell{mm,1} = seqTemp;
                firstSpkTime_Right_3(mm,:) = firstSpkTime_Right_Non_TimeCell(k,:);
            end
        end
        
        seqNum_Left_Non_TimeCell = length(seq_Left_Non_TimeCell);
        seqNum_Right_Non_TimeCell = length(seq_Right_Non_TimeCell);
        
        [~,template_Left_Non_TimeCell] = sort(mean(firstSpkTime_Left_3,1,'omitnan'));
        [~,template_Right_Non_TimeCell] = sort(mean(firstSpkTime_Right_3,1,'omitnan'));
        
        
        
         % sleep session
        firstSpkTime_Sleep = DelayPopFire_RegionSeq.sleep_eventTsp2Start;
        firstSpkTime_Sleep = firstSpkTime_Sleep(:,timeCellLabel);
        seqNum_Sleep = size(firstSpkTime_Sleep,1);
        seq_Sleep = cell(seqNum_Sleep,1);
        for k = 1:size(firstSpkTime_Sleep,1)
            [tsTemp,seqTemp] = sort(firstSpkTime_Sleep(k,:));
            seqTemp = seqTemp(~isnan(tsTemp));
            seq_Sleep{k} = seqTemp;
        end
        
        
        
        if seqNum_Left_TimeCell>5 && seqNum_Right_TimeCell>5
            % Left template matching
            [rval_LL,len_LL]=checktempseq(seq_Left_TimeCell, template_Left_TimeCell);
            [rval_LR,len_LR]=checktempseq(seq_Right_TimeCell, template_Left_TimeCell);
            [rval_LSleep,LSleep]=checktempseq(seq_Sleep, template_Left_TimeCell);
            
            % Right template matching
            [rval_RR,len_RR]=checktempseq(seq_Right_TimeCell, template_Right_TimeCell);
            [rval_RL,len_RL]=checktempseq(seq_Left_TimeCell, template_Right_TimeCell);
            [rval_RSleep,RSleep]=checktempseq(seq_Sleep, template_Right_TimeCell);
            
            % test left template significance
            % Ltemp to L significance
            Nt_TimeCell.delay.LL = r_Normalize(rval_LL,len_LL,rt);
            % Left temp to R significance
            Nt_TimeCell.delay.LR = r_Normalize(rval_LR,len_LR,rt);
            Nt_TimeCell.delay.LSleep = r_Normalize(rval_LSleep,LSleep,rt);
            
            % test left template significance
            % Rtemp to R significance
            Nt_TimeCell.delay.RR = r_Normalize(rval_RR,len_RR,rt);
            % Right temp to R significance
            Nt_TimeCell.delay.RL = r_Normalize(rval_RL,len_RL,rt);
            Nt_TimeCell.delay.RSleep = r_Normalize(rval_RSleep,RSleep,rt);
            
            figure
            Violin(Nt_TimeCell.delay.LL.zval,1);
            Violin(Nt_TimeCell.delay.LR.zval,2);
            Violin(Nt_TimeCell.delay.LSleep.zval,3);
            [h,pvAL]=ttest2(Nt_TimeCell.delay.LL.zval,Nt_TimeCell.delay.LR.zval);
            text(1.3,2,strcat('p = ',num2str(round(pvAL,3))));
            set(gca,'XTick',[1,2,3],'XTickLabel',{'L','R','Sleep'})
            TITLE1 = sprintf('%d%s%d',sessInfo(i).animal,' Day-',sessInfo(i).day);
            TITLE2 = 'Time cell Left turn template';
            title({TITLE1;TITLE2});
            
            figure
            Violin(Nt_TimeCell.delay.RR.zval,1);
            Violin(Nt_TimeCell.delay.RL.zval,2);
            Violin(Nt_TimeCell.delay.RSleep.zval,3);
            
            [h,pvAL]=ttest2(Nt_TimeCell.delay.RR.zval,Nt_TimeCell.delay.RL.zval);
            text(1.3,2,strcat('p = ',num2str(round(pvAL,3))));
            set(gca,'XTick',[1,2,3],'XTickLabel',{'R','L','Sleep'})
            TITLE1 = sprintf('%d%s%d',sessInfo(i).animal,' Day-',sessInfo(i).day);
            TITLE2 = 'Time cell Right turn template';
            title({TITLE1;TITLE2});
        end
        
        
        if seqNum_Left_Non_TimeCell>5 && seqNum_Right_Non_TimeCell>=5
            [rval_LL_NonTime,len_LL_NonTime]=checktempseq(seq_Left_Non_TimeCell, template_Left_Non_TimeCell);
            [rval_LR_NonTime,len_LR_NonTime]=checktempseq(seq_Right_Non_TimeCell, template_Left_Non_TimeCell);
            [rval_LSleep_NonTime,LSleep_NonTime]=checktempseq(seq_Sleep, template_Left_Non_TimeCell);
            
            [rval_RR_NonTime,len_RR_NonTime]=checktempseq(seq_Right_Non_TimeCell, template_Right_Non_TimeCell);
            [rval_RL_NonTime,len_RL_NonTime]=checktempseq(seq_Left_Non_TimeCell, template_Right_Non_TimeCell);
            [rval_RSleep_NonTime,RSleep_NonTime]=checktempseq(seq_Sleep, template_Right_Non_TimeCell);
            
            
            Nt_NonTimeCell.delay.LL = r_Normalize(rval_LL_NonTime,len_LL_NonTime,rt);
            Nt_NonTimeCell.delay.LR = r_Normalize(rval_LR_NonTime,len_LR_NonTime,rt);
            Nt_NonTimeCell.delay.LSleep = r_Normalize(rval_LSleep_NonTime,LSleep_NonTime,rt);
            
            Nt_NonTimeCell.delay.RR = r_Normalize(rval_RR_NonTime,len_RR_NonTime,rt);
            Nt_NonTimeCell.delay.RL = r_Normalize(rval_RL_NonTime,len_RL_NonTime,rt);
            Nt_NonTimeCell.delay.RSleep = r_Normalize(rval_RSleep_NonTime,RSleep_NonTime,rt);
            

            figure
            Violin(Nt_NonTimeCell.delay.LL.zval,1);
            Violin(Nt_NonTimeCell.delay.LR.zval,2);
            Violin(Nt_NonTimeCell.delay.LSleep.zval,3);
            
            [h,pvAL]=ttest2(Nt_NonTimeCell.delay.LL.zval,Nt_NonTimeCell.delay.LR.zval);
            text(1.3,2,strcat('p = ',num2str(round(pvAL,3))));
            set(gca,'XTick',[1,2,3],'XTickLabel',{'L','R','Sleep'})
            TITLE1 = sprintf('%d%s%d',sessInfo(i).animal,' Day-',sessInfo(i).day);
            TITLE2 = 'Non-Time cell Left turn template';
            title({TITLE1;TITLE2});
            
            figure
            Violin(Nt_NonTimeCell.delay.RR.zval,1);
            Violin(Nt_NonTimeCell.delay.RL.zval,2);
            Violin(Nt_NonTimeCell.delay.RSleep.zval,3);
            
            [h,pvAL]=ttest2(Nt_NonTimeCell.delay.RR.zval,Nt_NonTimeCell.delay.RL.zval);
            text(1.3,2,strcat('p = ',num2str(round(pvAL,3))));
            set(gca,'XTick',[1,2,3],'XTickLabel',{'R','L','Sleep'})
            TITLE1 = sprintf('%d%s%d',sessInfo(i).animal,' Day-',sessInfo(i).day);
            TITLE2 = 'Non-Time cell Right turn template';
            title({TITLE1;TITLE2});
            
            
            zmat_TimeCell.LL = [zmat_TimeCell.LL;Nt_TimeCell.delay.LL.zval];
            zmat_TimeCell.LR = [zmat_TimeCell.LR;Nt_TimeCell.delay.LR.zval];
            zmat_TimeCell.RR = [zmat_TimeCell.RR;Nt_TimeCell.delay.RR.zval];
            zmat_TimeCell.RL = [zmat_TimeCell.RL;Nt_TimeCell.delay.RL.zval];
  
            zmat_NonTimeCell.LL = [zmat_NonTimeCell.LL;Nt_NonTimeCell.delay.LL.zval];
            zmat_NonTimeCell.LR = [zmat_NonTimeCell.LR;Nt_NonTimeCell.delay.LR.zval];
            zmat_NonTimeCell.RR = [zmat_NonTimeCell.RR;Nt_NonTimeCell.delay.RR.zval];
            zmat_NonTimeCell.RL = [zmat_NonTimeCell.RL;Nt_NonTimeCell.delay.RL.zval];
            

            
        end                 
    fprintf('Finished analysis for session %d\n',i);
    clear DelayPopFire_RegionSeq on_eventTsp2Start off_eventTsp2Start sleep_eventTsp2Start
    end
end

figure
Violin(zmat_TimeCell.LL, 1);
Violin(zmat_TimeCell.LR, 2);
Violin(zmat_TimeCell.RR, 4);
Violin(zmat_TimeCell.RL, 5);
title('Time cell')

figure
Violin(zmat_NonTimeCell.LL, 1);
Violin(zmat_NonTimeCell.LR, 2);
Violin(zmat_NonTimeCell.RR, 4);
Violin(zmat_NonTimeCell.RL, 5);
title('Non-Time cell')

end