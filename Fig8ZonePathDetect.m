% This code is to define each trial start and end in each sub zone
% 
function Fig8ZonePathDetect(inFile,AnalyzeSes)
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

    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        % create the folder to save tracking images
        savedir = fullfile(sessInfo(i).mainDir,'TrackingImages');
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
               
            PathZone.rat = sessInfo(i).animal;
            PathZone.day = sessInfo(i).day;
            PathZone.sesID = SessDirs{j};
            
            % load analyzed positions
            posFile = fullfile(mainDir,SessDirs{j}, 'locInfo.mat');
            locInfo = load(posFile);
            pathFile = fullfile(mainDir,SessDirs{j}, 'pathData.mat');
            pathData = load(pathFile);
            
            [PathZone.posStartInd,PathZone.posEndInd,PathZone.posStartT,PathZone.posEndT] = armTime(locInfo,pathData);
         
            if p.writeToFile
                save(fullfile(sessInfo(i).mainDir,SessDirs{j},'PathZone.mat'), 'PathZone');
            end
            clear PathZone
            
            if p.savePlot == 1
                TITLE = sprintf('%s%d%s%d%s%s%s%s','Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',SessDirs{j},'-PathSubZoneDetection');
                title(TITLE,'Interpreter','None');
                figName = sprintf('%s%s%d%s%d%s%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',SessDirs{j},'-PathSubZoneTrialTime');
                print(figName,'-dpng','-r300');                                            
                close all
            end
        end


    end
    fprintf('Finished ratemap analysis for session %d\n',i);
end 

end

% get start and end of each run on center, choice, reward, return, base
% area respectively
% Li Yuan, 16-Oct-2020, UCSD
function [posStartInd,posEndInd,posStartT,posEndT] = armTime(locInfo,pathData)
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

%% find position idx on center arm
% range [2,13,5)
segLabel = zeros(length(locInfo.IDSeq),1);

segLabel(locInfo.IDSeq == 12 | locInfo.IDSeq == 7) = 1;
segLabel(locInfo.IDSeq == 2) = 3;
segLabel(locInfo.IDSeq == 13) = 6;
segLabel(locInfo.IDSeq == 5) = 10;
segLabel = [0;segLabel];

diffLabel = diff(segLabel);
% idx on IDSeq;
startidx = [];
endidx = [];
n = 1;
while n < length(diffLabel)-1
    if diffLabel(n) == 1 && diffLabel(n+1) == 2
        startidx = [startidx;n+1];
        m = 1;
        while (n+m) < length(diffLabel) && diffLabel(n+m)~= 4
            m = m+1;
        end
        endidx = [endidx;n+m-1];
        n = n+m;
    else
        n = n+1;
    end
end
if length(startidx) ~= length(endidx) || any((endidx-startidx) < 0)
    error('arm duration detection is wrong')
end
tempIdx = ((endidx-startidx)==0);
startidx(tempIdx)=[];
endidx(tempIdx)=[];

posStartInd.Center = find(ismember(pathData.t,locInfo.tInt(startidx,1))==1);
posEndInd.Center =  find(ismember(pathData.t,locInfo.tInt(endidx,2))==1);
posStartT.Center = pathData.t(posStartInd.Center);
posEndT.Center =  pathData.t(posEndInd.Center);
%% find position idx on choice arm
% range [5,9,4)
segLabel = zeros(length(locInfo.IDSeq),1);
segLabel(locInfo.IDSeq == 13) = 1;
segLabel(locInfo.IDSeq == 5) = 3;
segLabel(locInfo.IDSeq == 10 | locInfo.IDSeq == 9) = 6;
segLabel(locInfo.IDSeq == 4 | locInfo.IDSeq == 6) = 10;
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

if length(startidx) ~= length(endidx) || any((endidx-startidx) < 0)
    error('arm duration detection is wrong')
end
tempIdx = ((endidx-startidx)==0);
startidx(tempIdx)=[];
endidx(tempIdx)=[];

posStartInd.Choice = find(ismember(pathData.t,locInfo.tInt(startidx,1))==1);
posEndInd.Choice =  find(ismember(pathData.t,locInfo.tInt(endidx,2))==1);
posStartT.Choice = pathData.t(posStartInd.Choice);
posEndT.Choice =  pathData.t(posEndInd.Choice);
%% return arm position idx
% range [8/11,8/11]
segLabel = zeros(length(locInfo.IDSeq),1);
segLabel(locInfo.IDSeq == 4 | locInfo.IDSeq == 6) = 1;
segLabel(locInfo.IDSeq == 8 | locInfo.IDSeq == 11) = 3;
segLabel(locInfo.IDSeq == 3 | locInfo.IDSeq == 1) = 6;
segLabel = [0;segLabel];

diffLabel = diff(segLabel);
startidx = [];
endidx = [];

k = 1;
while k < length(diffLabel)-1
    if diffLabel(k) == 1 && diffLabel(k+1) == 2
        startidx = [startidx;k+1];
        m = 1;
        while (k+m) < length(diffLabel) && diffLabel(k+m)~= 3
            m = m+1;
        end
        endidx = [endidx;k+m];
        k = k+m;
    else
        k = k+1;
    end
end

% double check start and end
if length(startidx) ~= length(endidx) || any((endidx-startidx) < 0)
    error('arm duration detection is wrong')
end
tempIdx = ((endidx-startidx)==0);
startidx(tempIdx)=[];
endidx(tempIdx)=[];

posStartInd.Return = find(ismember(pathData.t,locInfo.tInt(startidx,1))==1);
posEndInd.Return =  find(ismember(pathData.t,locInfo.tInt(endidx,1))==1)-1;
posStartT.Return = pathData.t(posStartInd.Return);
posEndT.Return =  pathData.t(posEndInd.Return);
%% base arm position idx
segLabel = zeros(length(locInfo.IDSeq),1);
segLabel(locInfo.IDSeq == 8 | locInfo.IDSeq == 11) = 1;
segLabel(locInfo.IDSeq == 3 | locInfo.IDSeq == 1) = 3;
segLabel(locInfo.IDSeq == 7 | locInfo.IDSeq == 12) = 6;
segLabel(locInfo.IDSeq == 2) = 10;
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
if length(startidx) ~= length(endidx) || any((endidx-startidx) < 0)
    error('arm duration detection is wrong')
end
tempIdx = ((endidx-startidx)==0);
startidx(tempIdx)=[];
endidx(tempIdx)=[];

posStartInd.Base = find(ismember(pathData.t,locInfo.tInt(startidx,1))==1);
posEndInd.Base =  find(ismember(pathData.t,locInfo.tInt(endidx,2))==1);
posStartT.Base = pathData.t(posStartInd.Base);
posEndT.Base =  pathData.t(posEndInd.Base);
%% Reward position idx
% reward is after choice and before next return
posStartInd.Reward = posEndInd.Choice + 1;
posEndInd.Reward = posStartInd.Return(2:end) -1;
posEndInd.Reward(end+1) = length(pathData.t);

if any((posEndInd.Reward-posStartInd.Reward)<0)
    error('Reward time assignment wrong')
end
posStartT.Reward = pathData.t(posStartInd.Reward);
posEndT.Reward =  pathData.t(posEndInd.Reward);

% segLabel = zeros(length(locInfo.IDSeq),1);
% segLabel(locInfo.IDSeq == 9 | locInfo.IDSeq == 10) = 1;
% segLabel(locInfo.IDSeq == 4 | locInfo.IDSeq == 6) = 3;
% segLabel(locInfo.IDSeq == 8 | locInfo.IDSeq == 11) = 6;
% segLabel = [segLabel;6];
% 
% diffLabel = diff(segLabel);
% startidx = [];
% endidx = [];
% 
% k = 1;
% while k < length(diffLabel)-1
%     if diffLabel(k) == 1 && diffLabel(k+1) == 2
%         startidx = [startidx;k+1];
%         m = 1;
%         while (k+m) < length(diffLabel) && diffLabel(k+m)~= 3
%             m = m+1;
%         end
%         endidx = [endidx;k+m];
%         k = k+m;
%     else
%         k = k+1;
%     end
% end
% 
% startidx = startidx+1;
% endidx = endidx+1;
% if endidx(end) == length(locInfo.IDSeq)+1
%     endidx(end) = length(locInfo.IDSeq);
% end
%     
% % double check start and end
% if length(startidx) ~= length(endidx) || any((endidx-startidx) < 0)
%     error('arm duration detection is wrong')
% end
% tempIdx = ((endidx-startidx)==0);
% startidx(tempIdx)=[];
% endidx(tempIdx)=[];
% % % always session start from reward position
% % if startidx(1) > 2
% %     startidx = [1;startidx];
% %     endidx = [2;endidx];
% % end
% 
% posStartInd.Reward = find(ismember(pathData.t,locInfo.tInt(startidx,1))==1);
% posEndInd.Reward =  find(ismember(pathData.t,locInfo.tInt(endidx,1))==1)-1;
% posStartT.Reward = pathData.t(posStartInd.Reward);
% posEndT.Reward =  pathData.t(posEndInd.Reward);
% 
% if posEndInd.Reward(1) < posEndInd.Choice(1)
%     posStartInd.Reward(1) = [];
%     posEndInd.Reward(1) = [];
%     posStartT.Reward(1) = [];
%     posEndT.Reward(1) = [];
% end
% 
% tempInd = [];
% j = 0;
% for i = 1:length(posStartInd.Center)-1
%     j = j+1;
%     if posStartInd.Reward(i+1)<posStartInd.Return(j)
%         tempInd = [tempInd,i+1];
%         j = j-1;
%     end
% end
% 
% posStartInd.Reward(tempInd) = [];
% posEndInd.Reward(tempInd-1) = [];
% posStartT.Reward(tempInd) = [];
% posEndT.Reward(tempInd-1) = [];
%     
% 
% %%
% if length(unique([length(posStartInd.Reward),length(posStartInd.Return),length(posStartInd.Base),length(posStartInd.Choice),...
%         length(posStartInd.Center)]))>1   
%     error('Arm detection is wrong')
% end

% plot path to make sure
figure
for i = 1:length(posStartInd.Center)
    plot(pathData.x(posStartInd.Center(i):posEndInd.Center(i)),pathData.y(posStartInd.Center(i):posEndInd.Center(i)),'r');
    hold on
    plot(pathData.x(posStartInd.Choice(i):posEndInd.Choice(i)),pathData.y(posStartInd.Choice(i):posEndInd.Choice(i)),'g');
    plot(pathData.x(posStartInd.Reward(i):posEndInd.Reward(i)),pathData.y(posStartInd.Reward(i):posEndInd.Reward(i)),'k');
    plot(pathData.x(posStartInd.Return(i):posEndInd.Return(i)),pathData.y(posStartInd.Return(i):posEndInd.Return(i)),'y');
    plot(pathData.x(posStartInd.Base(i):posEndInd.Base(i)),pathData.y(posStartInd.Base(i):posEndInd.Base(i)),'b');
end
set (gca,'YDir','reverse')
legend('Center','Choice','Reward','Return','Base','Location','NorthEastOutside')
end
