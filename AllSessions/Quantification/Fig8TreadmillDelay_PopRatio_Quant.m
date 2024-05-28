% calculate the population activation ratio
%
% Li Yuan, 28-Mar-2022
function Fig8TreadmillDelay_PopRatio_Quant(inFile,AnalyzeSes)
close all

p.cellNumThres = 25;

% Read in input information
sessInfo = SessInfoImport(inFile);
validDay = 0;

for i = AnalyzeSes(1:end)
    
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');

    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);

    % load reactivation file
    popFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayPopFire_CellRatio_Delay.mat');
    load(popFile);
    
    % get valid cell ind
    % time cell / non-time cell
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;
    clusterNum = length(rateCluster);
    pyrNum = sum(rateLabel);
    
    if i == AnalyzeSes(1)
        % get each phase names (no delay etc)
        sessDirs = sessInfo(i).sessDirs;
        for j = 1:length(sessDirs)
            popRate.(sessDirs{j}) = [];
            cellRatio.(sessDirs{j}) = [];
        end
    end
    
    if pyrNum >= p.cellNumThres 
        validDay = validDay + 1;
        % get processed data from each subfolder
        for j = 1:length(sessDirs)
            sizeTemp = size(DelayPopFire_CellRatio_Delay.(sessDirs{j}).ratePopMat);
            popRateTemp = DelayPopFire_CellRatio_Delay.(sessDirs{j}).ratePopMat'./pyrNum;
            ratioTemp = DelayPopFire_CellRatio_Delay.(sessDirs{j}).cellPopMat'./pyrNum;
            popRate.(sessDirs{j}) = [popRate.(sessDirs{j}),reshape(popRateTemp,[1,sizeTemp(1)*sizeTemp(2)])];
            cellRatio.(sessDirs{j}) = [cellRatio.(sessDirs{j}),reshape(ratioTemp,[1,sizeTemp(1)*sizeTemp(2)])];
        end
    end   
end

validDay
ratioRegion = 0:0.03:1;
cellRatioDist.on10 = histcounts([cellRatio.on10_1,cellRatio.on10_2],ratioRegion);
cellRatioDist.off10 = histcounts([cellRatio.off10_1,cellRatio.off10_2],ratioRegion);
cellRatioDist.on30 = histcounts([cellRatio.on30_1,cellRatio.on30_2],ratioRegion);
cellRatioDist.off30 = histcounts([cellRatio.off30_1,cellRatio.off30_2],ratioRegion);

figure
plot(ratioRegion(2:end),100*cellRatioDist.on10/(sum(cellRatioDist.on10)));
hold on
plot(ratioRegion(2:end),100*cellRatioDist.off10/(sum(cellRatioDist.off10)));
plot(ratioRegion(2:end),100*cellRatioDist.on30/(sum(cellRatioDist.on30)));
plot(ratioRegion(2:end),100*cellRatioDist.off30/(sum(cellRatioDist.off30)));

on_Bin = (cellRatioDist.on10/sum(cellRatioDist.on10) + cellRatioDist.on30/sum(cellRatioDist.on30))/2;
off_Bin = (cellRatioDist.off10/sum(cellRatioDist.off10) + cellRatioDist.off30/sum(cellRatioDist.off30))/2;

figure
bar(ratioRegion(2:end),100*(on_Bin-off_Bin));
ylim([-6 6])
title('Occupation difference')

figure
bar(ratioRegion(2:11),(on_Bin(1:10)-off_Bin(1:10))./off_Bin(1:10));
title('Fold difference')

ratioRegion = 0:0.1:3;
popRateDist.on10 = histcounts([popRate.on10_1,popRate.on10_2],ratioRegion);
popRateDist.off10 = histcounts([popRate.off10_1,popRate.off10_2],ratioRegion);
popRateDist.on30 = histcounts([popRate.on30_1,popRate.on30_2],ratioRegion);
popRateDist.off30 = histcounts([popRate.off30_1,popRate.off30_2],ratioRegion);

figure
plot(ratioRegion(2:end),popRateDist.on10/(sum(popRateDist.on10)));
hold on
plot(ratioRegion(2:end),popRateDist.off10/(sum(popRateDist.off10)));
plot(ratioRegion(2:end),popRateDist.on30/(sum(popRateDist.on30)));
plot(ratioRegion(2:end),popRateDist.off30/(sum(popRateDist.off30)));

figure
plot(ratioRegion(2:end),(popRateDist.on10/(sum(popRateDist.on10)))-(popRateDist.off10/(sum(popRateDist.off10))));
hold on
plot(ratioRegion(2:end),(popRateDist.on30/(sum(popRateDist.on30)))-(popRateDist.off30/(sum(popRateDist.off30))));

end
