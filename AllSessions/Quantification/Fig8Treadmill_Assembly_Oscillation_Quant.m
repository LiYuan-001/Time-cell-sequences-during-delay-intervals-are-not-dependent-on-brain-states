
function Fig8Treadmill_Assembly_Oscillation_Quant(inFile,AnalyzeSes)
close all

AssemblySimilarityOn = [];
AssemblySimilarityOff = [];
assembly_ID_on = [];
assembly_ID_off = [];

on_eventLFP_onShare.swr = [];
on_eventLFP_offShare.swr = [];
on_eventLFP_onSpec.swr = [];
on_eventLFP_offSpec.swr = [];

off_eventLFP_onShare.swr = [];
off_eventLFP_offShare.swr = [];
off_eventLFP_onSpec.swr = [];
off_eventLFP_offSpec.swr = [];


on_eventLFP_onShare.strength_swr = [];
on_eventLFP_offShare.strength_swr = [];
on_eventLFP_onSpec.strength_swr = [];
on_eventLFP_offSpec.strength_swr = [];

off_eventLFP_onShare.strength_swr = [];
off_eventLFP_offShare.strength_swr = [];
off_eventLFP_onSpec.strength_swr = [];
off_eventLFP_offSpec.strength_swr = [];

on_eventLFP.swr = [];
on_eventLFP.strength_swr = [];
off_eventLFP.swr = [];
off_eventLFP.strength_swr = [];


% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
       
    assemblyLFPFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8_Assembly_OnOff_SWR_Rate-25ms.mat');
    load(assemblyLFPFile);
    Delay_OnOff_File = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_DelayOnOff_LR-25ms.mat');
    load(Delay_OnOff_File);
%     %% on blocks
%     if ~isempty(Fig8_Assembly_OnOff_SWR_Rate.on)
%         % get for all on blocks swr and strength
%         swrTemp = Fig8_Assembly_OnOff_SWR_Rate.on.on10_1.Assembly.rateDelay(:,4) + Fig8_Assembly_OnOff_SWR_Rate.on.on10_2.Assembly.rateDelay(:,4) + ...
%             Fig8_Assembly_OnOff_SWR_Rate.on.on30_1.Assembly.rateDelay(:,4) + Fig8_Assembly_OnOff_SWR_Rate.on.on30_2.Assembly.rateDelay(:,4);
%         eventNum = Fig8_Assembly_OnOff_SWR_Rate.on.on10_1.Assembly.rateDelay(:,3) + Fig8_Assembly_OnOff_SWR_Rate.on.on10_2.Assembly.rateDelay(:,3) + ...
%             Fig8_Assembly_OnOff_SWR_Rate.on.on30_1.Assembly.rateDelay(:,3) + Fig8_Assembly_OnOff_SWR_Rate.on.on30_2.Assembly.rateDelay(:,3);
%         strengthTemp = Fig8_Assembly_OnOff_SWR_Rate.on.on10_1.Assembly.rateDelay(:,5) + Fig8_Assembly_OnOff_SWR_Rate.on.on10_2.Assembly.rateDelay(:,5) + ...
%             Fig8_Assembly_OnOff_SWR_Rate.on.on30_1.Assembly.rateDelay(:,5) + Fig8_Assembly_OnOff_SWR_Rate.on.on30_2.Assembly.rateDelay(:,5);
%         
%         on_eventLFP.swr = [on_eventLFP.swr;swrTemp./eventNum];
%         on_eventLFP.strength_swr = [on_eventLFP.strength_swr;strengthTemp./swrTemp];
%     end
    
    
    %% off blocks
    if ~isempty(Fig8_Assembly_OnOff_SWR_Rate.off)
         % get for all on blocks swr and strength
        swrTemp = Fig8_Assembly_OnOff_SWR_Rate.off.off10_1.Assembly.rateDelay(:,4) + Fig8_Assembly_OnOff_SWR_Rate.off.off10_2.Assembly.rateDelay(:,4) + ...
            Fig8_Assembly_OnOff_SWR_Rate.off.off30_1.Assembly.rateDelay(:,4) + Fig8_Assembly_OnOff_SWR_Rate.off.off30_2.Assembly.rateDelay(:,4);
        eventNum = Fig8_Assembly_OnOff_SWR_Rate.off.off10_1.Assembly.rateDelay(:,3) + Fig8_Assembly_OnOff_SWR_Rate.off.off10_2.Assembly.rateDelay(:,3) + ...
            Fig8_Assembly_OnOff_SWR_Rate.off.off30_1.Assembly.rateDelay(:,3) + Fig8_Assembly_OnOff_SWR_Rate.off.off30_2.Assembly.rateDelay(:,3);
        strengthTemp = Fig8_Assembly_OnOff_SWR_Rate.off.off10_1.Assembly.rateDelay(:,5) + Fig8_Assembly_OnOff_SWR_Rate.off.off10_2.Assembly.rateDelay(:,5) + ...
            Fig8_Assembly_OnOff_SWR_Rate.off.off30_1.Assembly.rateDelay(:,5) + Fig8_Assembly_OnOff_SWR_Rate.off.off30_2.Assembly.rateDelay(:,5);
        
        off_eventLFP.swr = [off_eventLFP.swr;swrTemp./eventNum];
        off_eventLFP.strength_swr = [off_eventLFP.strength_swr;strengthTemp./swrTemp];
        
        swrTemp = Fig8_Assembly_OnOff_SWR_Rate.on.off10_1.Assembly.rateDelay(:,4) + Fig8_Assembly_OnOff_SWR_Rate.on.off10_2.Assembly.rateDelay(:,4) + ...
            Fig8_Assembly_OnOff_SWR_Rate.on.off30_1.Assembly.rateDelay(:,4) + Fig8_Assembly_OnOff_SWR_Rate.on.off30_2.Assembly.rateDelay(:,4);
        eventNum = Fig8_Assembly_OnOff_SWR_Rate.on.off10_1.Assembly.rateDelay(:,3) + Fig8_Assembly_OnOff_SWR_Rate.on.off10_2.Assembly.rateDelay(:,3) + ...
            Fig8_Assembly_OnOff_SWR_Rate.on.off30_1.Assembly.rateDelay(:,3) + Fig8_Assembly_OnOff_SWR_Rate.on.off30_2.Assembly.rateDelay(:,3);
        strengthTemp = Fig8_Assembly_OnOff_SWR_Rate.on.off10_1.Assembly.rateDelay(:,5) + Fig8_Assembly_OnOff_SWR_Rate.on.off10_2.Assembly.rateDelay(:,5) + ...
            Fig8_Assembly_OnOff_SWR_Rate.on.off30_1.Assembly.rateDelay(:,5) + Fig8_Assembly_OnOff_SWR_Rate.on.off30_2.Assembly.rateDelay(:,5);
        
        on_eventLFP.swr = [on_eventLFP.swr;swrTemp./eventNum];
        on_eventLFP.strength_swr = [on_eventLFP.strength_swr;strengthTemp./swrTemp];
    end
 
    for m = 1:CellAssembly_DelayLR.DelayOn.patNum
        sim_Temp = [];
        pattern_On = CellAssembly_DelayLR.DelayOn.AssmblWght(:,m);
        for j = 1:CellAssembly_DelayLR.DelayOff.patNum
            pattern_Off = CellAssembly_DelayLR.DelayOff.AssmblWght(:,j);
            sim_Temp(j) = dot(pattern_On,pattern_Off)/(norm(pattern_On)*norm(pattern_Off));
        end
        AssemblySimilarityOn = [AssemblySimilarityOn,max(sim_Temp)];
        assembly_ID_on = [assembly_ID_on,CellAssembly_DelayLR.rat*100+CellAssembly_DelayLR.day*10+m];
    end
    
    for m = 1:CellAssembly_DelayLR.DelayOff.patNum
        sim_Temp = [];
        pattern_Off = CellAssembly_DelayLR.DelayOff.AssmblWght(:,m);
        for j = 1:CellAssembly_DelayLR.DelayOn.patNum
            pattern_On = CellAssembly_DelayLR.DelayOn.AssmblWght(:,j);
            sim_Temp(j) = dot(pattern_Off,pattern_On)/(norm(pattern_Off)*norm(pattern_On));
        end
        AssemblySimilarityOff = [AssemblySimilarityOff,max(sim_Temp)];
        assembly_ID_off = [assembly_ID_off,CellAssembly_DelayLR.rat*100+CellAssembly_DelayLR.day*10+m];
    end
    
end

on_sim = AssemblySimilarityOn > 0.5;
off_sim = AssemblySimilarityOff > 0.5;

h = figure(1);
Violin(on_eventLFP.swr(on_sim),1,'ShowData',false);
Violin(on_eventLFP.swr(~on_sim),2,'ShowData',false);
Violin(off_eventLFP.swr(off_sim),3,'ShowData',false);
Violin(off_eventLFP.swr(~off_sim),4,'ShowData',false);
set(gca,'XTick',[1,2,3,4],'XtickLabel',{'On share','On Spec','Off share','Off Spec'})

h = figure(2);
Violin(on_eventLFP.strength_swr,1);
Violin(off_eventLFP.strength_swr,2);

subplot(2,4,1)
stdshade(on_eventLFP_onShare(:,4),0.5,[1,0,0]);
TITLE1 = 'On share assembly';
TITLE2 = sprintf('%s%d','On block, n = ',size(on_eventLFP_onShare,1));
title({TITLE1;TITLE2})
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

subplot(2,4,2)
stdshade(on_eventLFP_offShare(:,4),0.5,[1,0,0]);
TITLE1 = 'Off share assembly';
TITLE2 = sprintf('%s%d','On block, n = ',size(on_eventLFP_offShare,1));
title({TITLE1;TITLE2})
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

subplot(2,4,3)
stdshade(on_eventLFP_onSpec(:,4),0.5,[1,0,0]);
TITLE1 = 'On spec assembly';
TITLE2 = sprintf('%s%d','On block, n = ',size(on_eventLFP_onSpec,1));
title({TITLE1;TITLE2})
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

subplot(2,4,4)
stdshade(on_eventLFP_offSpec(:,4),0.5,[1,0,0]);
TITLE1 = 'Off spec assembly';
TITLE2 = sprintf('%s%d','On block, n = ',size(on_eventLFP_offSpec,1));
title({TITLE1;TITLE2})
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])


subplot(2,4,5)
stdshade(off_eventLFP_onShare(:,4),0.5,[0,0,0]);
TITLE2 = sprintf('%s%d','Off block, n = ',size(off_eventLFP_onShare,1));
title(TITLE2)
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

subplot(2,4,6)
stdshade(off_eventLFP_offShare(:,4),0.5,[0,0,0]);
TITLE2 = sprintf('%s%d','Off block, n = ',size(off_eventLFP_offShare,1));
title(TITLE2)
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

subplot(2,4,7)
stdshade(off_eventLFP_onSpec(:,4),0.5,[0,0,0]);
TITLE2 = sprintf('%s%d','Off block, n = ',size(off_eventLFP_onSpec,1));
title(TITLE2)
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

subplot(2,4,8)
stdshade(off_eventLFP_offSpec(:,4),0.5,[0,0,0]);
TITLE2 = sprintf('%s%d','Off block, n = ',size(off_eventLFP_offSpec,1));
title(TITLE2)
set(gca,'XTick',[1 500 1000],'XtickLabel',[-250 0 250])
xlabel('ms');
ylim([-50 50])

end