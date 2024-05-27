% Read in sequences and form input for checkmotseq.m
% run after code function checkmotseq_SeqReadin(inFile,AnalyzeSes) and rt=rtfun2_LiModified(narr)
% Nov-19-2021, Li Yuan, UCSD
function checkTemp_Batch(inFile,AnalyzeSes)


p.writeToFile = 1;
p.cellNum = 20;
p.shuffle = 100;

% Read in input information
sessInfo = SessInfoImport(inFile);

rtFile = 'W:\Li Yuan\Codes\Fig8MazeTreadmill\Cell Property\ReactivationCode-Leibold\matlab code\Results\rt.mat';
load(rtFile);
  
zmat.LL = [];
zmat.LR = [];
zmat.RR = [];
zmat.RL = [];

for i = AnalyzeSes(1:end)
    close all

    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    % load reactivation file
    seqFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_RegionSeq.mat');
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    
    % get valid cell ind
    % rate < 5 Hz in whole session
    % time cell / non-time cell
    clusterNum = length(SpikeProp.max_AvgRate);
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    
    if sum(rateLabel) >= p.cellNum
        load(seqFile);
        % calculate template then test template
        
        leftInd = DelayPopFire_RegionSeq.off_eventTsp2Start.delayLeftInd;
        rightInd = ~DelayPopFire_RegionSeq.off_eventTsp2Start.delayLeftInd;
        firstSpkTime_Left = DelayPopFire_RegionSeq.off_eventTsp2Start.delay(leftInd==1,:);
        firstSpkTime_Right = DelayPopFire_RegionSeq.off_eventTsp2Start.delay(rightInd==1,:);
        
        firstSpkTime_Left = firstSpkTime_Left(:,rateLabel);
        firstSpkTime_Right = firstSpkTime_Right(:,rateLabel);
        

        seqNum_Left = sum(leftInd);
        seqNum_Right = sum(rightInd);
                        
        firstSpkTime_Left = firstSpkTime_Left(:,sum(~isnan(firstSpkTime_Left),1)>=seqNum_Left/5);
        firstSpkTime_Right = firstSpkTime_Right(:,sum(~isnan(firstSpkTime_Right),1)>=seqNum_Right/5);
        
        
        seq_Left = cell(seqNum_Left,1);
        for k = 1:size(firstSpkTime_Left,1)
            [tsTemp,seqTemp] = sort(firstSpkTime_Left(k,:));
            seqTemp = seqTemp(~isnan(tsTemp));
            seq_Left{k} = seqTemp;
        end

        seq_Right = cell(seqNum_Right,1);
        for k = 1:size(firstSpkTime_Right,1)
            [tsTemp,seqTemp] = sort(firstSpkTime_Right(k,:));
            seqTemp = seqTemp(~isnan(tsTemp));
            seq_Right{k} = seqTemp;
        end
        
         % sleep session
        firstSpkTime_Sleep = DelayPopFire_RegionSeq.sleep_eventTsp2Start;
        seqNum_Sleep = size(firstSpkTime_Sleep,1);
        seq_Sleep = cell(seqNum_Sleep,1);
        for k = 1:size(firstSpkTime_Sleep,1)
            [tsTemp,seqTemp] = sort(firstSpkTime_Sleep(k,:));
            seqTemp = seqTemp(~isnan(tsTemp));
            seq_Sleep{k} = seqTemp;
        end
        
        
        if seqNum_Left>=10 && seqNum_Right>=10
            
        rval_LL = [];
        len_LL = [];
            for m = 1:seqNum_Left
                seqID = 1:seqNum_Left;
                seqID = find(seqID~=m);                
                [~,template_Left] = sort(mean(firstSpkTime_Left(seqID,:),1,'omitnan'));                
                % Left template matching
                [rval_LL(m),len_LL(m)]=checktempseq(seq_Left(m), template_Left);
                [rval_LR,len_LR]=checktempseq(seq_Right, template_Left);
                [rval_LSleep,LSleep]=checktempseq(seq_Sleep, template_Left);
            end
            
            rval_RR = [];
            len_RR = [];
            for m = 1:seqNum_Right
                seqID = 1:seqNum_Right;
                seqID = find(seqID~=m);
                [~,template_Right] = sort(mean(firstSpkTime_Right(seqID,:),1,'omitnan'));    
                % Right template matching
                [rval_RR(m),len_RR(m)]=checktempseq(seq_Right(m), template_Right);
                [rval_RL,len_RL]=checktempseq(seq_Left, template_Right);
                [rval_RSleep,RSleep]=checktempseq(seq_Sleep, template_Right);
            end
%             % shuffle Left
%             nrep = seqNum_Left;
%             [rvalR_Shuffle,rvalL_Shuffle,lenR_Shuffle,lenL_Shuffle]=checktempshuffle(seq_Left, template_Left, template_Right, 100);
            
            % test left template significance
            % Ltemp to L significance
            Nt.delay.LL = r_Normalize(rval_LL,len_LL,rt);
            % Left temp to R significance
            Nt.delay.LR = r_Normalize(rval_LR,len_LR,rt);
            Nt.delay.LSleep = r_Normalize(rval_LSleep,LSleep,rt);
            
            
            
            figure
            Violin(Nt.delay.LL.zval,1);
            Violin(Nt.delay.LR.zval,2);
            Violin(Nt.delay.LSleep.zval,3);
            [h,pvAL]=ttest2(Nt.delay.LL.zval,Nt.delay.LR.zval);
            text(1.3,2,strcat('p = ',num2str(round(pvAL,3))));
            set(gca,'XTick',[1,2,3],'XTickLabel',{'L','R','Sleep'})
            TITLE1 = sprintf('%d%s%d',sessInfo(i).animal,' Day-',sessInfo(i).day);
            TITLE2 = 'Left turn template';
            title({TITLE1;TITLE2});
            
            
            % test left template significance
            % Rtemp to R significance
            Nt.delay.RR = r_Normalize(rval_RR,len_RR,rt);
            % Right temp to R significance
            Nt.delay.RL = r_Normalize(rval_RL,len_RL,rt);
            Nt.delay.RSleep = r_Normalize(rval_RSleep,RSleep,rt);
            
            
            figure
            Violin(Nt.delay.RR.zval,1);
            Violin(Nt.delay.RL.zval,2);
            Violin(Nt.delay.RSleep.zval,3);
            
            [h,pvAL]=ttest2(Nt.delay.RR.zval,Nt.delay.RL.zval);
            text(1.3,2,strcat('p = ',num2str(round(pvAL,3))));
            set(gca,'XTick',[1,2,3],'XTickLabel',{'R','L','Sleep'})
            TITLE1 = sprintf('%d%s%d',sessInfo(i).animal,' Day-',sessInfo(i).day);
            TITLE2 = 'Right turn template';
            title({TITLE1;TITLE2});
            
%             zmat.LL = [zmat.LL,Nt.delay.LL.zval];
%             zmat.LR = [zmat.LR,Nt.delay.LR.zval];
%             zmat.RR = [zmat.RR,Nt.delay.RR.zval];
%             zmat.RL = [zmat.RL,Nt.delay.RL.zval];

           
        end                 
    fprintf('Finished analysis for session %d\n',i);
    clear DelayPopFire_RegionSeq on_eventTsp2Start off_eventTsp2Start sleep_eventTsp2Start
    end
end

figure
Violin(zmat.LL, 1);
Violin(zmat.LR, 2);
Violin(zmat.RR, 4);
Violin(zmat.RL, 5);

end