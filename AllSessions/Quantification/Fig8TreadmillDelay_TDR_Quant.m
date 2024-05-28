% Quantify theta length and break timein delay area
% Li Yuan, Feb-22-22, UCSD
%
function Fig8TreadmillDelay_TDR_Quant(inFile,AnalyzeSes)
close all

p.savePlot = 0;
p.writeToFile = 0;
% Read in input information
sessInfo = SessInfoImport(inFile);
sessDirs = {'on10','off10','on30','off30'};

for i = AnalyzeSes(1:end)
    
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    eegCh = (sessInfo(i).EEGch2);
    
    tdrFile = fullfile(mainDir,'Cell Property', 'ThetaDeltaRatioAll.mat');
    load(tdrFile);
    eegCh_matrix = nan(length(ThetaDeltaRatioAll.EEGch),1);
    for m = 1:length(eegCh_matrix)
        eegCh_matrix(m) = str2num(ThetaDeltaRatioAll.EEGch{m});
    end
    [~,eegCh_matrix_Ind] = ismember(eegCh,eegCh_matrix);
    
    for j = 1:length(sessDirs)
        if i == AnalyzeSes(1)
            tdr.(sessDirs{j}) = [];
        end
        tdr_meanTemp = [];
        
        % two blocks
        for n = 1:2
            sess{n} = strcat(sessDirs{j},'_',num2str(n));
            
            pathZoneFile = fullfile(mainDir,sess{n}, 'PathZone.mat');
            load(pathZoneFile);
            delayFile = fullfile(mainDir,sess{n}, 'Fig8DelayZonePos.mat');
            load(delayFile);
            
            timeStamp = ThetaDeltaRatioAll.(sess{n}).timeStamp;
            % get from every channel being chosen
            tdrTemp = nanmean(ThetaDeltaRatioAll.(sess{n}).tdr_Smooth(eegCh_matrix_Ind,:),1);
            
            delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
            delayTend1 = Fig8DelayZonePos.delayPos1.endT;
            delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
            delayTend2 = Fig8DelayZonePos.delayPos2.endT;
            
            trialNum = size(delayTstart1,2);
            if contains(sess{n},'10')
                maxT = 10;
                delayTend1_2 = delayTstart1+maxT;
                delayTend2_2 = delayTstart2+maxT;
            elseif contains(sess{n},'30')
                maxT = 29;
                delayTend1_2 = delayTstart1+maxT;
                delayTend2_2 = delayTstart2+maxT;
            else
                error('Delay time is wrong')
            end
            timeGap = timeStamp(2) - timeStamp(1);
            indLength = round(maxT./timeGap);
            
            [~,startInd] = min(abs(delayTstart1 - timeStamp));
            
            % get for every trial
            for k = 1:length(startInd)
                tdr_meanTemp = [tdr_meanTemp;tdrTemp(startInd(k):startInd(k)+indLength)];
            end
        end
        tdr.(sessDirs{j}) = [tdr.(sessDirs{j});mean(tdr_meanTemp,1)];
    end
end

% figure
% Violin(thetaLength_on10,1,'ViolinColor',[0.85,0.33,0.1]);
% Violin(thetaLength_off10,2,'ViolinColor',[0.3,0.74,0.93]);
% Violin(thetaLength_on30,3,'ViolinColor',[0.64,0.08,0.18]);
% Violin(thetaLength_off30,4,'ViolinColor',[0,0.44,0.74]);
%
% figure
% Violin(thetaBreak_on10,1,'ViolinColor',[0.85,0.33,0.1]);
% Violin(thetaBreak_off10,2,'ViolinColor',[0.3,0.74,0.93]);
% Violin(thetaBreak_on30,3,'ViolinColor',[0.64,0.08,0.18]);
% Violin(thetaBreak_off30,4,'ViolinColor',[0,0.44,0.74]);

h = figure;
h.Position = [100,100,1200,900];
% meanVal = mean(tdr.on10,1);
% ste = std(tdr.on10,1)./sqrt(size(tdr.on10,1)-1);
% seqLen = size(tdr.on10,2);
% plot(meanVal,'r')
% hold on
% shade([(meanVal-ste)',(meanVal+ste)'],'FillType',[1,2;2,1],'FillColor','r');

stdshade(tdr.on10,0.5,[1,0,0]);
hold on
stdshade(tdr.off10,0.5,[0,0,0]);
ylim([0 5])


h = figure;
h.Position = [100,100,1200,900];
stdshade(tdr.on30,0.5,[1,0,0]);
hold on
stdshade(tdr.off30,0.5,[0,0,0]);
ylim([0 5])

end