% Read in sequences and form input for checkmotseq.m
% run after code function checkmotseq_SeqReadin(inFile,AnalyzeSes) and rt=rtfun2_LiModified(narr)
% Nov-19-2021, Li Yuan, UCSD
function checkmotseq_Batch(inFile,AnalyzeSes)

close all

p.writeToFile = 0;

% Read in input information
sessInfo = SessInfoImport(inFile);

rtFile = 'W:\Li Yuan\Codes\Fig8MazeTreadmill\Cell Property\ReactivationCode-Leibold\matlab code\Results\rt2.mat';
load(rtFile);
    
for i = AnalyzeSes(1:end)
    
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    % load reactivation file
    seqFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_RegionSeq.mat');

    if isfile(seqFile)
        load(seqFile);
        
%         DelayPopFire_RegionSeq.rat = sessInfo(i).animal;
%         DelayPopFire_RegionSeq.day = sessInfo(i).day;
%         DelayPopFire_RegionSeq.clusterNum = clusterNum;
%         DelayPopFire_RegionSeq.PyrCellNum = sum(rateLabel);
%         
        
%         on_eventTsp2Start.delay 
%         on_eventTsp2Start.reward 
%         
%         on_eventTsp2Start.delayLeft 
%         on_eventTsp2Start.delayRight 
%         on_eventTsp2Start.rewardLeft 
%         on_eventTsp2Start.rewardRight 
%                 
%         off_eventTsp2Start.delay
%         off_eventTsp2Start.reward
% 
%         off_eventTsp2Start.delayLeft 
%         off_eventTsp2Start.delayRight 
%         off_eventTsp2Start.rewardLeft
%         off_eventTsp2Start.rewardRight
%         sleep_eventTsp2Start 
            
        
        % checkmotseq then testmot
        % off delay as whole
        leftInd = DelayPopFire_RegionSeq.off_eventTsp2Start.delayLeftInd;
        rightInd = ~DelayPopFire_RegionSeq.off_eventTsp2Start.delayLeftInd;
        firstSpkTime_Left = DelayPopFire_RegionSeq.off_eventTsp2Start.delay(leftInd==1,:);
        firstSpkTime_Right = DelayPopFire_RegionSeq.off_eventTsp2Start.delay(rightInd==1,:);
        
        firstSpkTime = [firstSpkTime_Right;firstSpkTime_Left];
%         firstSpkTime = [DelayPopFire_RegionSeq.off_eventTsAllspks.delayRight;DelayPopFire_RegionSeq.off_eventTsAllspks.delayLeft];
        rightNum = size(firstSpkTime_Right,1);
        leftNum = size(firstSpkTime_Left,1);
        turnDir1 = [zeros(rightNum,1);ones(leftNum,1)];
        
        seqNum = size(firstSpkTime,1);
        seq1 = cell(seqNum,1);
        for k = 1:size(firstSpkTime,1)
            [tsTemp,seqTemp] = sort(firstSpkTime(k,:));
            seqTemp = seqTemp(~isnan(tsTemp));
            seq1{k} = seqTemp;
        end
        % put right turn seq first and concatenate left turn seq
        [turnDir2,seqOrder] = sort(turnDir1);
        seq2 = seq1(seqOrder);
        
        [rval.offDelay,len.offDelay]=checkmotseq_LiModified(seq2);
        N.delay = testmot_LiModified(rval.offDelay,len.offDelay,rt);

%         a=[];
%         heatMat = zeros(size(N.delay.zmat));
%         for kk=1:seqNum-1
%             for kkk = kk+1:seqNum
%                 a = [a,N.delay.zmat(kk,kkk)];
%                 if turnDir2(kk)+turnDir2(kkk) == 0
%                     heatMat(kk,kkk) = 1;
%                     heatMat(kkk,kk) = 1;
%                 elseif turnDir2(kk)+turnDir2(kkk) == 1
%                     heatMat(kk,kkk) = 2;
%                     heatMat(kkk,kk) = 2;
%                 elseif turnDir2(kk)+turnDir2(kkk) == 2
%                     heatMat(kk,kkk) = 3;
%                     heatMat(kkk,kk) = 3;
%                 end
%             end     
%         end
%         figure
%         Violin(a,1);
%         figure
%         temp = zeros(size(heatMat));
%         temp((N.delay.zmat)>N.delay.sigZ) = heatMat((N.delay.zmat)>N.delay.sigZ);
%         heatmap(temp,'ColorLimits',[0 3]);  
% 
%         cmap = [0.5, 0.5, 0.5; ...   % gray
%             1, 0, 0; ...       % red
%             0, 1, 0; ...       % Green
%             0,0,1];            % blue
        
%         colormap(cmap)

        % plot social network
        figure
        cmap = [1, 0, 0;0,0,1];
        weight = zeros(size(N.delay.zmat));
        % find >95% and <5%
        weight((N.delay.zmat)>N.delay.sigZ) = (N.delay.zmat((N.delay.zmat)>N.delay.sigZ));
        counts = sum(weight~=0);
        counts(counts==0) = 1;
        G = graph(weight);
        LWidths = 2*G.Edges.Weight/max(G.Edges.Weight);
        h = plot(G,'NodeLabel',[1:seqNum],'LineWidth',LWidths);
        for nn = 1:seqNum
            highlight(h,nn, 'MarkerSize', counts(nn),'NodeColor',cmap(turnDir2(nn)+1,:))           
        end
        title('Delay event social network')
        
%         % plot chord diagram
%         groups = 1:seqNum;
%         connections = [];
%         sigcount = 0;
%         for kk=1:seqNum-1
%             for kkk = kk+1:seqNum
%                 if (N.delay.zmat(kk,kkk))>N.delay.sigZ
%                     sigcount = sigcount+1;
%                     connections(sigcount,1) = kk;
%                     connections(sigcount,2) = kkk;
%                     connections(sigcount,3) = (N.delay.zmat(kk,kkk));
%                     if turnDir2(kk)+turnDir2(kkk) == 0
%                         connections(sigcount,4) = 1;
%                     elseif turnDir2(kk)+turnDir2(kkk) == 1
%                         connections(sigcount,4) = 2;
%                     elseif turnDir2(kk)+turnDir2(kkk) == 2
%                         connections(sigcount,4) = 3;
%                     end
%                 end
%             end     
%         end
%         connectionTable = array2table(connections);
%         p = chordPlot(groups, connectionTable);
%         for nn = 1:seqNum
%             highlight(p,nn, 'MarkerSize', 2*counts(nn),'NodeColor',cmap(turnDir2(nn)+1,:))
%         end
%         title('Delay event chord diagram')
        
        % delay and reward
        % off delay as whole
        firstSpkTime = [DelayPopFire_RegionSeq.off_eventTsp2Start.delay;DelayPopFire_RegionSeq.off_eventTsp2Start.reward;DelayPopFire_RegionSeq.on_eventTsAllspks.reward];
%         firstSpkTime = [DelayPopFire_RegionSeq.off_eventTsAllspks.delayRight;DelayPopFire_RegionSeq.off_eventTsAllspks.delayLeft];
        
        seqNum = size(firstSpkTime,1);
        seqNum_Delay = size(DelayPopFire_RegionSeq.off_eventTsp2Start.delay,1);
        seqNum_Reward = size([DelayPopFire_RegionSeq.off_eventTsp2Start.reward;DelayPopFire_RegionSeq.on_eventTsAllspks.reward],1);
        
        seq1 = cell(seqNum,1);
        for k = 1:size(firstSpkTime,1)
            [tsTemp,seqTemp] = sort(firstSpkTime(k,:));
            seqTemp = seqTemp(~isnan(tsTemp));
            seq1{k} = seqTemp;
        end

        [rval.DelayRewrad,len.DelayRewrad]=checkmotseq_LiModified(seq1);
        N.DelayRewrad = testmot_LiModified(rval.DelayRewrad,len.DelayRewrad,rt);

%         a=[];
%         heatMat = zeros(size(N.DelayRewrad.zmat));
%         for kk=1:seqNum-1
%             for kkk = kk+1:seqNum
%                 a = [a,N.DelayRewrad.zmat(kk,kkk)];
%                 if kk<=seqNum_Delay && kkk<=seqNum_Delay
%                     heatMat(kk,kkk) = 1;
%                     heatMat(kkk,kk) = 1;
%                 elseif kk<=seqNum_Delay && kkk>seqNum_Delay
%                     heatMat(kk,kkk) = 2;
%                     heatMat(kkk,kk) = 2;
%                 elseif kk>seqNum_Delay && kkk>seqNum_Delay
%                     heatMat(kk,kkk) = 3;
%                     heatMat(kkk,kk) = 3;
%                 end
%             end     
%         end
%         figure
%         Violin(a,1);
%         figure
%         temp = zeros(size(heatMat));
%         temp((N.DelayRewrad.zmat)>N.DelayRewrad.sigZ) = heatMat((N.DelayRewrad.zmat)>N.DelayRewrad.sigZ);
%         heatmap(temp,'ColorLimits',[0 3]);  
% 
%         cmap = [0.5, 0.5, 0.5; ...   % gray
%             1, 0, 0; ...       % red
%             0, 1, 0; ...       % Green
%             0,0,1];            % blue
%         
%         colormap(cmap)

        % plot social network
        figure
        cmap = [1, 0, 0;0,0,1];
        weight = zeros(size(N.DelayRewrad.zmat));
        % find >95% and <5%
        weight((N.DelayRewrad.zmat)>N.DelayRewrad.sigZ) = (N.DelayRewrad.zmat((N.DelayRewrad.zmat)>N.DelayRewrad.sigZ));
        counts = sum(weight~=0);
        counts(counts==0) = 1;
        G = graph(weight);
        LWidths = 2*G.Edges.Weight/max(G.Edges.Weight);
        h = plot(G,'NodeLabel',[1:seqNum],'LineWidth',LWidths);
        for nn = 1:seqNum_Delay
            highlight(h,nn, 'MarkerSize', counts(nn),'NodeColor','r')           
        end
        for nn = seqNum_Delay+1:seqNum_Delay+seqNum_Reward
            highlight(h,nn, 'MarkerSize', counts(nn),'NodeColor','b')           
        end
        title('Delay and reward event social network')
        
%         % plot chord diagram
%         groups = 1:seqNum;
%         connections = [];
%         sigcount = 0;
%         for kk=1:seqNum-1
%             for kkk = kk+1:seqNum
%                 if (N.DelayRewrad.zmat(kk,kkk))>N.DelayRewrad.sigZ
%                     sigcount = sigcount+1;
%                     connections(sigcount,1) = kk;
%                     connections(sigcount,2) = kkk;
%                     connections(sigcount,3) = (N.DelayRewrad.zmat(kk,kkk));
%                     if kk<=seqNum_Delay && kkk<=seqNum_Delay
%                         connections(sigcount,4) = 1;
%                     elseif kk<=seqNum_Delay && kkk>seqNum_Delay
%                         connections(sigcount,4) = 2;
%                     elseif kk>seqNum_Delay && kkk>seqNum_Delay
%                         connections(sigcount,4) = 3;
%                     end
%                 end
%             end     
%         end
%         connectionTable = array2table(connections);
%         p = chordPlot(groups, connectionTable);
%         for nn = 1:seqNum_Delay
%             highlight(p,nn, 'MarkerSize', 2*counts(nn),'NodeColor','r')
%         end
%         for nn = seqNum_Delay+1:seqNum_Delay+seqNum_Reward
%             highlight(p,nn, 'MarkerSize', 2*counts(nn),'NodeColor','b')          
%         end
%         title('Delay reward event chord diagram')
%         
        
        
        
        % sleep session
        firstSpkTime = [DelayPopFire_RegionSeq.sleep1_eventTsp2Start;DelayPopFire_RegionSeq.sleep2_eventTsp2Start]; 
        seqNum = size(firstSpkTime,1);
        seqNum_Sleep1 = size(DelayPopFire_RegionSeq.sleep1_eventTsp2Start,1);
        seqNum_Sleep2 = size(DelayPopFire_RegionSeq.sleep2_eventTsp2Start,1);
        
        seq1 = cell(seqNum,1);
        for k = 1:size(firstSpkTime,1)
            [tsTemp,seqTemp] = sort(firstSpkTime(k,:));
            seqTemp = seqTemp(~isnan(tsTemp));
            seq1{k} = seqTemp;
        end
        [rval.sleep,len.sleep]=checkmotseq_LiModified(seq1);
        N.sleep = testmot_LiModified(rval.sleep,len.sleep,rt);
        a=[];
        for kk=1:seqNum-1
            for kkk = kk+1:seqNum
                a = [a,N.sleep.zmat(kk,kkk)];
            end
        end
        
        figure
        Violin(a,2);
        
        % plot social network
        figure
        cmap = [1, 0, 0;0,0,1];
        weight = zeros(size(N.sleep.zmat));
        % find >95% and <5%
        weight((N.sleep.zmat)>N.sleep.sigZ) = (N.sleep.zmat((N.sleep.zmat)>N.sleep.sigZ));
        counts = sum(weight~=0);
        counts(counts==0) = 1;
        G = graph(weight);
        LWidths = G.Edges.Weight/max(G.Edges.Weight);
        h = plot(G,'NodeLabel',[1:seqNum],'LineWidth',LWidths);
        for nn = 1:seqNum_Sleep1
            highlight(h,nn, 'MarkerSize', counts(nn)./max(counts)*15,'NodeColor','r')
        end
        for nn = seqNum_Sleep1+1:seqNum_Sleep1+seqNum_Sleep2
            highlight(h,nn, 'MarkerSize', counts(nn)./max(counts)*15,'NodeColor','b')
        end
        title('Sleep event social network')
        
%         % plot chord diagram
%         figure
%         groups = 1:seqNum;
%         connections = [];
%         sigcount = 0;
%         for kk=1:seqNum-1
%             for kkk = kk+1:seqNum
%                 if (N.sleep.zmat(kk,kkk))>N.sleep.sigZ
%                     sigcount = sigcount+1;
%                     connections(sigcount,1) = kk;
%                     connections(sigcount,2) = kkk;
%                     connections(sigcount,3) = (N.sleep.zmat(kk,kkk));
%                     connections(sigcount,4) = 3;
%                 end
%             end     
%         end
%         connectionTable = array2table(connections);
%         p = chordPlot(groups, connectionTable);
%         for nn = 1:seqNum_Sleep1
%             highlight(p,nn, 'MarkerSize', counts(nn)./max(counts)*10,'NodeColor','r')
%         end
%         for nn = 1:seqNum_Sleep1+1:seqNum_Sleep1+seqNum_Sleep2
%             highlight(p,nn, 'MarkerSize', counts(nn)./max(counts)*10,'NodeColor','b')
%         end
%  
%         title('Sleep event chord diagram')
        
    end
    
    fprintf('Finished analysis for session %d\n',i);
    clear DelayPopFire_RegionSeq on_eventTsp2Start off_eventTsp2Start sleep_eventTsp2Start
end

    

end