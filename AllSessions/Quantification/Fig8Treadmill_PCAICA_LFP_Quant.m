
function Fig8Treadmill_PCAICA_LFP_Quant(inFile,AnalyzeSes)
close all

on_eventLFP_onShare = [];
on_eventLFP_offShare = [];
on_eventLFP_onSpec = [];
on_eventLFP_offSpec = [];

off_eventLFP_onShare = [];
off_eventLFP_offShare = [];
off_eventLFP_onSpec = [];
off_eventLFP_offSpec = [];

% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
       
    assemblyLFPFile = fullfile(sessInfo(i).mainDir,'Cell Property','Assembly_LFPShape.mat');
    load(assemblyLFPFile);
    
    %% on blocks
    if ~isempty(Assembly_LFPShape.on.eventLFP_onShare)
        for k = 1:length(Assembly_LFPShape.on.eventLFP_onShare)
            on_eventLFP_onShare = [on_eventLFP_onShare;Assembly_LFPShape.on.eventLFP_onShare{k}];
        end       
    end
    
    if ~isempty(Assembly_LFPShape.on.eventLFP_offShare)
        for k = 1:length(Assembly_LFPShape.on.eventLFP_offShare)
            on_eventLFP_offShare = [on_eventLFP_offShare;Assembly_LFPShape.on.eventLFP_offShare{k}];
        end       
    end
    
    if ~isempty(Assembly_LFPShape.on.eventLFP_onSpec)
        for k = 1:length(Assembly_LFPShape.on.eventLFP_onSpec)
            on_eventLFP_onSpec = [on_eventLFP_onSpec;Assembly_LFPShape.on.eventLFP_onSpec{k}];
        end       
    end
    
    if ~isempty(Assembly_LFPShape.on.eventLFP_offSpec)
        for k = 1:length(Assembly_LFPShape.on.eventLFP_offSpec)
            on_eventLFP_offSpec = [on_eventLFP_offSpec;Assembly_LFPShape.on.eventLFP_offSpec{k}];
        end       
    end
    
    %% off blocks
    if ~isempty(Assembly_LFPShape.off.eventLFP_onShare)
        for k = 1:length(Assembly_LFPShape.off.eventLFP_onShare)
            off_eventLFP_onShare = [off_eventLFP_onShare;Assembly_LFPShape.off.eventLFP_onShare{k}];
        end       
    end
    
    if ~isempty(Assembly_LFPShape.off.eventLFP_offShare)
        for k = 1:length(Assembly_LFPShape.off.eventLFP_offShare)
            off_eventLFP_offShare = [off_eventLFP_offShare;Assembly_LFPShape.off.eventLFP_offShare{k}];
        end       
    end
    
    if ~isempty(Assembly_LFPShape.off.eventLFP_onSpec)
        for k = 1:length(Assembly_LFPShape.off.eventLFP_onSpec)
            off_eventLFP_onSpec = [off_eventLFP_onSpec;Assembly_LFPShape.off.eventLFP_onSpec{k}];
        end       
    end
    
    if ~isempty(Assembly_LFPShape.off.eventLFP_offSpec)
        for k = 1:length(Assembly_LFPShape.off.eventLFP_offSpec)
            off_eventLFP_offSpec = [off_eventLFP_offSpec;Assembly_LFPShape.off.eventLFP_offSpec{k}];
        end       
    end
 
end

h = figure(1);
h.Position = [100,100,1200,800];

subplot(2,4,1)
stdshade(on_eventLFP_onShare,0.5,[1,0,0]);
TITLE1 = 'On share assembly';
TITLE2 = sprintf('%s%d','On block, n = ',size(on_eventLFP_onShare,1));
title({TITLE1;TITLE2})
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

subplot(2,4,2)
stdshade(on_eventLFP_offShare,0.5,[1,0,0]);
TITLE1 = 'Off share assembly';
TITLE2 = sprintf('%s%d','On block, n = ',size(on_eventLFP_offShare,1));
title({TITLE1;TITLE2})
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

subplot(2,4,3)
stdshade(on_eventLFP_onSpec,0.5,[1,0,0]);
TITLE1 = 'On spec assembly';
TITLE2 = sprintf('%s%d','On block, n = ',size(on_eventLFP_onSpec,1));
title({TITLE1;TITLE2})
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

subplot(2,4,4)
stdshade(on_eventLFP_offSpec,0.5,[1,0,0]);
TITLE1 = 'Off spec assembly';
TITLE2 = sprintf('%s%d','On block, n = ',size(on_eventLFP_offSpec,1));
title({TITLE1;TITLE2})
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])


subplot(2,4,5)
stdshade(off_eventLFP_onShare,0.5,[0,0,0]);
TITLE2 = sprintf('%s%d','Off block, n = ',size(off_eventLFP_onShare,1));
title(TITLE2)
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

subplot(2,4,6)
stdshade(off_eventLFP_offShare,0.5,[0,0,0]);
TITLE2 = sprintf('%s%d','Off block, n = ',size(off_eventLFP_offShare,1));
title(TITLE2)
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

subplot(2,4,7)
stdshade(off_eventLFP_onSpec,0.5,[0,0,0]);
TITLE2 = sprintf('%s%d','Off block, n = ',size(off_eventLFP_onSpec,1));
title(TITLE2)
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

subplot(2,4,8)
stdshade(off_eventLFP_offSpec,0.5,[0,0,0]);
TITLE2 = sprintf('%s%d','Off block, n = ',size(off_eventLFP_offSpec,1));
title(TITLE2)
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

figure
subplot(1,2,1)
stdshade([on_eventLFP_onShare;on_eventLFP_onSpec],0.5,[1,0,0]);
TITLE1 = 'On assembly';
TITLE2 = sprintf('%s%d','On block, n = ',size(on_eventLFP_onShare,1)+size(on_eventLFP_onSpec,1));
title({TITLE1;TITLE2})
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

subplot(1,2,2)
stdshade([off_eventLFP_offShare;off_eventLFP_offSpec],0.5,[0,0,0]);
TITLE1 = 'Off assembly';
TITLE2 = sprintf('%s%d','Off block, n = ',size(off_eventLFP_offShare,1)+size(off_eventLFP_offSpec,1));
title({TITLE1;TITLE2})
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

end