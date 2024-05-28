function Fig8Treadmill_DelayPopRate_Quant(inFile,AnalyzeSes)
close all
% Read in input information
sessInfo = SessInfoImport(inFile);

sessDirs2 = {'on10','off10','on30','off30'};

delayRate.on10 = [];
delayRate.off10 = [];
delayRate.on30 = [];
delayRate.off30 = [];

delayRatio.on10 = [];
delayRatio.off10 = [];
delayRatio.on30 = [];
delayRatio.off30 = [];

for i = AnalyzeSes(1:end)

    fileName = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8TreadmillDelayRatePop.mat');
    load(fileName);
    delayRate.on10 = [delayRate.on10,Fig8TreadmillDelayRatePop.on10.delayRate];
    delayRate.off10 = [delayRate.off10,Fig8TreadmillDelayRatePop.off10.delayRate];
    delayRate.on30 = [delayRate.on30,Fig8TreadmillDelayRatePop.on30.delayRate];
    delayRate.off30 = [delayRate.off30,Fig8TreadmillDelayRatePop.off30.delayRate];

    delayRatio.on10 = [delayRatio.on10,Fig8TreadmillDelayRatePop.on10.delayRatio];
    delayRatio.off10 = [delayRatio.off10,Fig8TreadmillDelayRatePop.off10.delayRatio];
    delayRatio.on30 = [delayRatio.on30,Fig8TreadmillDelayRatePop.on30.delayRatio];
    delayRatio.off30 = [delayRatio.off30,Fig8TreadmillDelayRatePop.off30.delayRatio];

end

h = figure;
% h.Position = [100,100,1200,900];
stdshade(delayRate.on10',0.5,[1,0,0]);
hold on
stdshade(delayRate.off10',0.5,[0,0,0]);


h = figure;
stdshade(delayRate.on30',0.5,[1,0,0]);
hold on
stdshade(delayRate.off30',0.5,[0,0,0]);

h = figure;
stdshade(delayRatio.on10',0.5,[1,0,0]);
hold on
stdshade(delayRatio.off10',0.5,[0,0,0]);


h = figure;
stdshade(delayRatio.on30',0.5,[1,0,0]);
hold on
stdshade(delayRatio.off30',0.5,[0,0,0]);

end
