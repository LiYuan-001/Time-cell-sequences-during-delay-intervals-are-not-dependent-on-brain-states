close all
% order: on10 off10 on30 off30
delayActive_mean = [32.83,35.61,30.58,35.93];
delayActive_sem = [3.39,5.07,5.15,2.98];
mazeActive_mean = [46.54,42.63,42.50,43.90];
mazeActive_sem = [5.64,5.28,6.31,5.64];
figure
subplot(1,3,1)
bar(1:4,mazeActive_mean)
hold on
errorbar(1:4,mazeActive_mean,mazeActive_sem,mazeActive_sem)
bar(6:9,delayActive_mean)
errorbar(6:9,delayActive_mean,delayActive_sem,delayActive_sem)
set(gca, 'XTick', [1:4,6:9], 'XTickLabel', {'on10','off10','on30','off30'});
title('OnOff group maze and delay active cells %')
ylim([0 100])


% order: on10 on30 off10 off3019.07

timeCell_mean = [20.07,14.99,19.07,16.42];
timeCell_sem = [3.38,4.37,1.87,1.87];
subplot(1,3,2)
bar(1:4,timeCell_mean)
hold on
errorbar(1:4,timeCell_mean,timeCell_sem,timeCell_sem)
set(gca, 'XTick', [1:4], 'XTickLabel', {'on10','on30','off10','off30'});
title('Time cells %')
ylim([0 30])

% ctr group
timeCell_Ctr_mean = [13.96,9.14,13.60,10.98];
timeCell_Ctr_sem = [2.73,3.15,2.71,3.09];
subplot(1,3,3)
bar(1:4,timeCell_Ctr_mean)
hold on
errorbar(1:4,timeCell_Ctr_mean,timeCell_Ctr_sem,timeCell_Ctr_sem)
set(gca, 'XTick', [1:4], 'XTickLabel', {'on10','on30','off10','off30'});
title('Ctr Time cells %')
ylim([0 30])