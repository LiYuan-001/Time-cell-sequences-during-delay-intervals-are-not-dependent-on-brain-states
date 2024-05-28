%data order is single cell: on10,on30,off10,off30; all assembly: on10,off10,on30,off30; share assembly:on10,on30,off10,off30;
% specific assembly: on10,on30,off10,off30;
delaySelective = [10,11,9,11,4,4,8,9,3,3,4,3,1,1,4,6];
delayAll = [193,193,193,193,61,64,95,105,43,43,43,43,18,21,52,62];
delay_pout = myBinomTest(delaySelective,delayAll,0.05,'one');

stemSelective = [12,14,18,18,2,3,6,10,1,2,4,2,1,1,2,8];
stemAll = [147,147,147,147,24,24,27,30,12,14,14,14,12,10,13,16];
stem_pout = myBinomTest(stemSelective,stemAll,0.05,'one');

choiceSelective = [77,81,88,86,6,10,17,20,3,6,7,11,3,4,10,9];
choiceAll = [148,148,148,148,10,15,28,32,5,10,11,17,5,5,17,15];
choice_pout = myBinomTest(choiceSelective,choiceAll,0.05,'one');

figure
bar(100*delaySelective./delayAll);
title('Delay selective: single cell, all assemblies, shared assemblies, specific assemblies')
set(gca,'XTick',[1:16],'XTickLabel',{'on10','on30','off10','off30'});
 
figure
bar(100*stemSelective./stemAll);
title('Stem selective: single cell, all assemblies, shared assemblies, specific assemblies')
set(gca,'XTick',[1:16],'XTickLabel',{'on10','on30','off10','off30'});

figure
bar(100*choiceSelective./choiceAll);
title('Choice selective: single cell, all assemblies, shared assemblies, specific assemblies')
set(gca,'XTick',[1:16],'XTickLabel',{'on10','on30','off10','off30'});