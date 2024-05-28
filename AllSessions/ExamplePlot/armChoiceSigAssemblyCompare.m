% plot significant off 30 stem choice assemblies
% % Read in input information
% stem detected ID: 104422      104431      104441      104442      104451      104633      107911      107932      107942
% delay detected ID:  104423      104433      104434      104441      104445      104454      104456      107939      107941      107946
sessInfo = SessInfoImport('V:\LiYuan\Codes\Fig8MazeTreadmill_V2\Fig8Treadmill_OnOff.xlsx');

figure


i = 6;
stemFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_Stem-25ms.mat');
load(stemFile);
delayFile =  fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_DelayOnOff_LR-25ms.mat');
load(delayFile);

stemAssembly = CellAssembly_Stem.DelayOff.AssmblWght(:,2);
delayAssembly = CellAssembly_DelayLR.DelayOff.AssmblWght(:,3);
sim = dot(stemAssembly,delayAssembly)/(norm(stemAssembly)*norm(delayAssembly));
subplot(8,2,1)
stem(delayAssembly,'k')
view([90 -90])
ylim([-1 1])

subplot(8,2,2)
stem(stemAssembly,'b')
view([90 -90])
ylim([-1 1])
TITLE = sprintf('%s%1.2f','Sim: ',sim);
title(TITLE);

i = 7;
stemFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_Stem-25ms.mat');
load(stemFile);
delayFile =  fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_DelayOnOff_LR-25ms.mat');
load(delayFile);

stemAssembly = CellAssembly_Stem.DelayOff.AssmblWght(:,1);
delayAssembly = CellAssembly_DelayLR.DelayOff.AssmblWght(:,[3,4]);

subplot(8,3,4)
stem(delayAssembly(:,1),'k')
view([90 -90])
ylim([-1 1])
sim = dot(stemAssembly,delayAssembly(:,1))/(norm(stemAssembly)*norm(delayAssembly(:,1)));
TITLE = sprintf('%s%1.2f','Sim: ',sim);
title(TITLE);

subplot(8,3,5)
stem(delayAssembly(:,2),'k')
view([90 -90])
ylim([-1 1])
sim = dot(stemAssembly,delayAssembly(:,2))/(norm(stemAssembly)*norm(delayAssembly(:,2)));
TITLE = sprintf('%s%1.2f','Sim: ',sim);
title(TITLE);

subplot(8,3,6)
stem(stemAssembly,'b')
view([90 -90])
ylim([-1 1])
