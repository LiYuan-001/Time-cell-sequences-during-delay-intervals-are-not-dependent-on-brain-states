function plotTimeCellMap_2Block(mapSession1_Sort_1,mapSession2_Sort_1,mapSession1_Sort_Norm_1,mapSession2_Sort_Norm_1,...
       mapSession1_Sort_2,mapSession2_Sort_2,mapSession1_Sort_Norm_2,mapSession2_Sort_Norm_2,sessInfo,block1,block2)

   
   pos1 = [0.07 0.6 0.18 0.3];
   pos2 = [0.28 0.6 0.18 0.3];
   pos3 = [0.07 0.2 0.18 0.3];
   pos4 = [0.28 0.2 0.18 0.3];
   pos5 = [0.54 0.6 0.18 0.3];
   pos6 = [0.75 0.6 0.18 0.3];
   pos7 = [0.54 0.2 0.18 0.3];
   pos8 = [0.75 0.2 0.18 0.3];
      
subplot('Position',pos1)
imagesc(mapSession1_Sort_1)
colormap(jet)
TITLE1 = sprintf('%s%d%s%d%s%s%s%s%s%s','DelayZone1-',sessInfo.animal,' Day-',sessInfo.day,'-',block1);
title({TITLE1},'Interpreter','None')
xlabel('Time (Sec)')
ylabel('Cells')
xlim([0 size(mapSession1_Sort_1,2)])
xTick= [0 size(mapSession1_Sort_1,2)];
xTickLabel = [0 10];
set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);

subplot('Position',pos2)
imagesc(mapSession2_Sort_1)
colormap(jet)
TITLE1 = sprintf('%s%d%s%d%s%s%s%s%s%s','DelayZone1-',sessInfo.animal,' Day-',sessInfo.day,'-',block2);
title({TITLE1},'Interpreter','None')
xlabel('Time (Sec)')
ylabel('Cells')
xlim([0 size(mapSession1_Sort_1,2)])
xTick= [0 size(mapSession1_Sort_1,2)];
xTickLabel = [0 10];
set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);

subplot('Position',pos3)
imagesc(mapSession1_Sort_Norm_1)
colormap(jet)
TITLE1 = sprintf('%s%d%s%d%s%s%s%s%s%s','DelayZone1-',sessInfo.animal,' Day-',sessInfo.day,'-',block1);
TITLE2 = 'Normalized to each cell peak';
title({TITLE1;TITLE2},'Interpreter','None')
xlabel('Time (Sec)')
ylabel('Cells')
xlim([0 size(mapSession1_Sort_1,2)])
xTick= [0 size(mapSession1_Sort_1,2)];
xTickLabel = [0 10];
set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);

subplot('Position',pos4)
imagesc(mapSession2_Sort_Norm_1)
colormap(jet)
TITLE1 = sprintf('%s%d%s%d%s%s%s%s%s%s','DelayZone1-',sessInfo.animal,' Day-',sessInfo.day,'-',block2);
TITLE2 = 'Normalized to each cell peak';
title({TITLE1;TITLE2},'Interpreter','None')
xlabel('Time (Sec)')
ylabel('Cells')
xlim([0 size(mapSession1_Sort_1,2)])
xTick= [0 size(mapSession1_Sort_1,2)];
xTickLabel = [0 10];
set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);

annotationString = sprintf('%s%s','Cell picked sorted by ',block1);
annotation('textbox',[0.025, 0.025, 0.45, 0.95], 'string',annotationString, 'Interpreter','None','FontSize',14,'LineWidth',2)

subplot('Position',pos5)
imagesc(mapSession1_Sort_2)
colormap(jet)
TITLE1 = sprintf('%s%d%s%d%s%s%s%s%s%s','DelayZone1-',sessInfo.animal,' Day-',sessInfo.day,'-',block2);
title({TITLE1},'Interpreter','None')
xlabel('Time (Sec)')
ylabel('Cells')
xlim([0 size(mapSession1_Sort_2,2)])
xTick= [0 size(mapSession1_Sort_2,2)];
xTickLabel = [0 10];
set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);

subplot('Position',pos6)
imagesc(mapSession2_Sort_2)
colormap(jet)
TITLE1 = sprintf('%s%d%s%d%s%s%s%s%s%s','DelayZone1-',sessInfo.animal,' Day-',sessInfo.day,'-',block1);
title({TITLE1},'Interpreter','None')
xlabel('Time (Sec)')
ylabel('Cells')
xlim([0 size(mapSession1_Sort_2,2)])
xTick= [0 size(mapSession1_Sort_2,2)];
xTickLabel = [0 10];
set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);

subplot('Position',pos7)
imagesc(mapSession1_Sort_Norm_2)
colormap(jet)
TITLE1 = sprintf('%s%d%s%d%s%s%s%s%s%s','DelayZone1-',sessInfo.animal,' Day-',sessInfo.day,'-',block2);
TITLE2 = 'Normalized to each cell peak';
title({TITLE1;TITLE2},'Interpreter','None')
xlabel('Time (Sec)')
ylabel('Cells')
xlim([0 size(mapSession1_Sort_2,2)])
xTick= [0 size(mapSession1_Sort_2,2)];
xTickLabel = [0 10];
set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);

subplot('Position',pos8)
imagesc(mapSession2_Sort_Norm_2)
colormap(jet)
TITLE1 = sprintf('%s%d%s%d%s%s%s%s%s%s','DelayZone1-',sessInfo.animal,' Day-',sessInfo.day,'-',block1);
TITLE2 = 'Normalized to each cell peak';
title({TITLE1;TITLE2},'Interpreter','None')
xlabel('Time (Sec)')
ylabel('Cells')
xlim([0 size(mapSession1_Sort_2,2)])
xTick= [0 size(mapSession1_Sort_2,2)];
xTickLabel = [0 10];
set(gca, 'XTick', xTick, 'XTickLabel', xTickLabel, 'TickLength', [0, 0]);
annotationString = sprintf('%s%s','Cell picked sorted by ',block2);
annotation('textbox',[0.5, 0.025, 0.45, 0.95], 'string',annotationString, 'Interpreter','None','FontSize',14,'LineWidth',2)

end