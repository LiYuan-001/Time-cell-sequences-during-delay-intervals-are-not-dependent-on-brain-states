% This code is to define delay zone in two ywas
% 1. Rat totally enters delay box as start
% 2. Rat reaches the delay area as start
function Fig8TreadDelayZoneDef(inFile,AnalyzeSes)
%   
%  Li Yuan, UCSD, 15-Nov-2019
%  plot general map for fig8 maze nd open box session
% -------------------------------------------------------------------------
% set parameters for analysis

p.savePlot = 1;
p.writeToFile = 1;

% -------------------------------------------------------------------------

% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
   close all
    if p.savePlot
        % Li Yuan, UCSD, 15-Nov-2019
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day);
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end

    display(['About to do session ' sessInfo(i).mainDir]);
    % get each phase names (no delay etc)
    SessDirs = sessInfo(i).sessDirs;
    mainDir = sessInfo(i).mainDir;
    
    if prod(~isempty(char(SessDirs{1}))) 
        
        for j = 1:length(SessDirs)
            % load analyzed positions
            posFile = fullfile(mainDir,SessDirs{j}, 'locInfo.mat');
            locInfo = load(posFile);
            pathFile = fullfile(mainDir,SessDirs{j}, 'pathData.mat');
            pathData = load(pathFile);
            % method 1
            % count start when animal totally enters delay area
            % detect aniaml on Arm A25, then use first Ypos < 30 as start
            
            % method 2
            % count start when animal approaches delay area
            % detect animal on Arm A25, use start of consecutive A25 pos as start

            % find time period in each segment on maze
            % label: N1 N2 N3 N4 N5 N6 A23 A34 A45 A56 A16 A12 A25
            % IDseq: 1  2  3  4  5  6   7   8   9  10  11  12  13
            % regToUse = {'center','choice','base', 'return', 'reward', 'choiceEnd'};
            % regLocs.center = {'A25'};
            % regLocs.choice = {'N5'};
            % regLocs.base = {'A56','A45'};
            % regLocs.return = {'A16','N1','A12','A23','N3','A34'};
            % regLocs.reward = {'N6','N4'};
            % regLocs.choiceEnd = {'N2'};
            
            % choice time segment and matching speed
            % find animal on the center arm
            % position label could go wrong if animal moves a lot
            % duration correction
            % -----------------------------------------------------------------
            segLabel = zeros(length(locInfo.IDSeq),1);
            % find first N2 and Last N25 for one pass
            segLabel(locInfo.IDSeq == 12 | locInfo.IDSeq == 7) = 1;
            segLabel(locInfo.IDSeq == 2) = 3;
            segLabel(locInfo.IDSeq == 13) = 6;
            segLabel(locInfo.IDSeq == 5) = 10;
            segLabel = [0;segLabel];
            
            diffLabel = diff(segLabel);
            startidx = [];
            endidx = [];
            k = 1;
            while k < length(diffLabel)-1
                if diffLabel(k) == 1 && diffLabel(k+1) == 2
                    startidx = [startidx;k+1];
                    m = 1;
                    while (k+m) < length(diffLabel) && diffLabel(k+m)~= 4
                        m = m+1;
                    end
                    endidx = [endidx;k+m-1];
                    k = k+m;
                else
                    k = k+1;
                end
            end
            % double check start and end
            if length(startidx) ~= length(endidx) || any(endidx < startidx)
                error('arm duration detection is wrong')
            end
            tempIdx = ((endidx-startidx)==0);
            startidx(tempIdx)=[];
            endidx(tempIdx)=[];
            % -----------------------------------------------------------------
            % method 1
            % delay start from delay barrier
            delayPos1 = posFind1(pathData,locInfo,startidx,endidx);
            % method 2
            % delay start from delay entrance
            delayPos2 = posFind2(pathData,locInfo,startidx,endidx);           
            
            % for .mat files which is going to be write down
            Fig8DelayZonePos.mainDir = mainDir;
            Fig8DelayZonePos.session = SessDirs{j};
            Fig8DelayZonePos.delayPos1 = delayPos1;
            Fig8DelayZonePos.delayPos2 = delayPos2;            
         
            if p.writeToFile
                save(fullfile(sessInfo(i).mainDir,SessDirs{j},'Fig8DelayZonePos.mat'), 'Fig8DelayZonePos');
            end
            
            % plot delay positions to confirm the detection is correct
            delayTstart1 = delayPos1.startT;
            delayTstart2 = delayPos2.startT;
            if contains(SessDirs{j},'10')
                maxT = 10;
                delayTend1_2 = delayTstart1+maxT;
                delayTend2_2 = delayTstart2+maxT;
            elseif contains(SessDirs{j},'30')
                maxT = 25;
                delayTend1_2 = delayTstart1+maxT;
                delayTend2_2 = delayTstart2+maxT;
            else
                error('Delay time is wrong')
            end
            figure(1)
            subplot(2,length(SessDirs)/2,j)
            [~, startInd] = min(abs(pathData.t-delayTstart1));
            [~, endInd] = min(abs(pathData.t-delayTend1_2));
            plot(pathData.x,pathData.y,'b')
            hold on
            axis tight
            axis off
            for k = 1:length(startInd)
                plot(pathData.x(startInd(k):endInd(k)),pathData.y(startInd(k):endInd(k)),'r')
            end
            set(gca,'YDir','Reverse')
            if j == 1
                TITLE1 = sprintf('%s%d%s%d%s%s','DelayZone1-',sessInfo(i).animal,'-Day-',sessInfo(i).day,'-',SessDirs{j},'-DelayDef');
            else
                TITLE1 = SessDirs{j};
            end
            title({TITLE1},'Interpreter','None')
            
            
            clear Fig8DelayZonePos
        end
   
        if p.savePlot == 1
            figure(1)
            figName = sprintf('%s%s%d%s%d%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-DelayZoneDefinition');
            print(figName,'-dpng','-r300');
        end
    close all

    end
    fprintf('Finished ratemap analysis for session %d\n',i);
end 

end
