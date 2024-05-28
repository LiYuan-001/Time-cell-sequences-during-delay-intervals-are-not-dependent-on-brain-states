function Fig8Treadmill_Task_Accuracy(inFile,AnalyzeSes)

close all

p.savePlot = 0;
p.writeToFile = 0;

p.avgRateThres = 1;

% % Read in input information
sessInfo = SessInfoImport(inFile);

% if p.savePlot
%     % directory for plot figures
%     % generate a folder for each rat eah day under the current folder
%     savedir = sprintf('%s%s',cd,'\Figures\Quantification Figures\Time cell num-dist quant');
%     if ~exist(savedir, 'dir')
%         mkdir(savedir);
%     end
% end

for i = AnalyzeSes(1:end)

    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;    
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    for j = 1:length(sessDirs)
        if i == AnalyzeSes(1)
            Accuracy.(sessDirs{j}) = [];
        end
        
        % load turning directions and correctness
        tInfoFile = fullfile(mainDir,sessDirs{j}, 'trialInfo.mat');
        tInfo = load(tInfoFile);
        accuracyTemp = sum(tInfo.success)/length(tInfo.success);
        Accuracy.(sessDirs{j}) = [Accuracy.(sessDirs{j}),accuracyTemp];
    end  
end

% on 10, off 10, on 30, off 30
Accuracy_on10 = 100*mean([Accuracy.on10_1;Accuracy.on10_2]);
Accuracy_off10 = 100*mean([Accuracy.off10_1;Accuracy.off10_2]);
Accuracy_on30 = 100*mean([Accuracy.on30_1;Accuracy.on30_2]);
Accuracy_off30 = 100*mean([Accuracy.off30_1;Accuracy.off30_2]);

figure(1)
% find peak and sort
plot([Accuracy_on10;Accuracy_off10;Accuracy_on30;Accuracy_off30],'d-','Color',[0.8,0.8,0.8])
hold on
errorbar(1,mean(Accuracy_on10),std(Accuracy_on10)/sqrt(length(Accuracy_on10)),'ro')
errorbar(2,mean(Accuracy_off10),std(Accuracy_off10)/sqrt(length(Accuracy_off10)),'ko')
errorbar(3,mean(Accuracy_on30),std(Accuracy_on30)/sqrt(length(Accuracy_on30)),'ro')
errorbar(4,mean(Accuracy_off30),std(Accuracy_off30)/sqrt(length(Accuracy_off30)),'ko')
% boxplot(Accuracy_on10,'Position',1);
% boxplot(Accuracy_off10,'Position',2);
% boxplot(Accuracy_on30,'Position',3);
% boxplot(Accuracy_off30,'Position',4);
xlim([0 5])
ylim([0 100])
title('Turning accuracy')
set(gca, 'XTick', [1,2,3,4], 'XTickLabel', {'on10','off10','on30','off30'});
end