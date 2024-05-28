
function Fig8Treadmill_OnOffRateDiff_PCAICA(inFile,AnalyzeSes)
close all
figureCount = 0;
ArmPCAChoice.rate = [];

sessDirs = {'on10','off10','on30','off30'};
assemblyNum_All = 0;

for j = 1:length(sessDirs)
    PCA_Rate.(sessDirs{j}) = [];
    PCA_Strength.(sessDirs{j}) = [];
end
    
% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    
    
    Fig8ArmRateFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8TreadmillArmRate_SameDelay_DelayAssembly.mat');
    if exist(Fig8ArmRateFile)
        load(Fig8ArmRateFile);
    
        for j = 1:length(sessDirs)
            % rate comparison
            PCA_Rate.(sessDirs{j}) = [PCA_Rate.(sessDirs{j}),Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).rateDelay];
            PCA_Strength.(sessDirs{j}) = [PCA_Strength.(sessDirs{j}),Fig8TreadmillArmRate_SameDelay_DelayAssembly.(sessDirs{j}).eventStrengthDelay];
        end
    end
end

figure
Violin(PCA_Rate.on10,1)
Violin(PCA_Rate.off10,2)
Violin(PCA_Rate.on30,3)
Violin(PCA_Rate.off30,4)

delayActive_on10 = PCA_Rate.on10>0.05;
delayActive_off10 = PCA_Rate.off10>0.05;
delayActive_on30 = PCA_Rate.on30>0.05;
delayActive_off30 = PCA_Rate.off30>0.05;

figure
% plot wenn of the delay on cell distribution
setListData = {find(delayActive_on10==1); find(delayActive_off10==1); find(delayActive_on30==1); find(delayActive_off30==1)};
setLabels = ["on 10"; "off 10"; "on 30"; "off 30"];
h = vennEulerDiagram(setListData, setLabels, 'drawProportional', true);

figure
plot([PCA_Rate.on10-PCA_Rate.off10],[PCA_Rate.on30-PCA_Rate.off30],'ko')
% fit off session
x = [PCA_Rate.on10-PCA_Rate.off10];
y1 = [PCA_Rate.on30-PCA_Rate.off30];
[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'k-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
pText = sprintf('%s%1.2f%s%1.3f','r2 = ',r_Val,' pVal = ',p_val);
text(0,-1,pText,'Color',[0,0,0]);


figure
plot([PCA_Strength.on10-PCA_Strength.off10],[PCA_Strength.on30-PCA_Strength.off30],'ko')
% fit off session
x = [PCA_Strength.on10-PCA_Strength.off10];
y1 = [PCA_Strength.on30-PCA_Strength.off30];
xNaN = isnan(x);
x = x(~xNaN);
y1 = y1(~xNaN);
yNaN = isnan(y1);
x = x(~yNaN);
y1 = y1(~yNaN);

[P,S] = polyfit(x,y1,1);
yfit = P(1)*x+P(2);
hold on;
plot(x,yfit,'k-');
lm = fitlm(x,y1);
p_val = lm.coefTest; 
r_Val = lm.Rsquared.Ordinary;
pText = sprintf('%s%1.2f%s%1.3f','r2 = ',r_Val,' pVal = ',p_val);
text(50,-50,pText,'Color',[0,0,0]);

end