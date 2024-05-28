
function Fig8Treadmill_ThetaLock_Assembly_Quant(inFile,AnalyzeSes)
close all
p.numT = 30;
% for theta lock figure plot
MinH = 0;
MaxH = 0.4;
xDetail = 0:360;
MaxSpike = 3000;    % max spike number for raster
MaxT = 0.4;

sessDirs = {'on10','on30'};

for j = 1:length(sessDirs)
    assembly_ppcT.(sessDirs{j})= [];
    assembly_rT.(sessDirs{j})= [];
    assembly_angT.(sessDirs{j})= [];
    assembly_numT.(sessDirs{j})= [];
    assembly_pT.(sessDirs{j})= [];
    
    ppcT.(sessDirs{j})= [];
    rT.(sessDirs{j})= [];
    angT.(sessDirs{j})= [];
    numT.(sessDirs{j})= [];
    pT.(sessDirs{j})= [];
end

% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    
    LockFile = fullfile(sessInfo(i).mainDir,'Cell Property','ThetaLock_Assembly_2Session.mat');
    load(LockFile);
    LockFile = fullfile(sessInfo(i).mainDir,'Cell Property','ThetaLock_2Session.mat');
    load(LockFile);
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.AvgRate.Fig8Rate;
    rateLabel = rateCluster<5 & rateCluster>0.1;    
    
    if ThetaLock_Assembly_2Session.AssemblyNum > 0
        for j = 1:length(sessDirs)
            assembly_ppcT.(sessDirs{j}) = [assembly_ppcT.(sessDirs{j}),ThetaLock_Assembly_2Session.(sessDirs{j}).thetaLock_delay.ppcT];
            assembly_rT.(sessDirs{j}) = [assembly_rT.(sessDirs{j}),ThetaLock_Assembly_2Session.(sessDirs{j}).thetaLock_delay.rT];
            assembly_angT.(sessDirs{j}) = [assembly_angT.(sessDirs{j}),ThetaLock_Assembly_2Session.(sessDirs{j}).thetaLock_delay.angT];
            assembly_numT.(sessDirs{j}) = [assembly_numT.(sessDirs{j}),ThetaLock_Assembly_2Session.(sessDirs{j}).thetaLock_delay.numT];
            assembly_pT.(sessDirs{j})= [assembly_pT.(sessDirs{j}),ThetaLock_Assembly_2Session.(sessDirs{j}).thetaLock_delay.pT];
            
            ppcT.(sessDirs{j})= [ppcT.(sessDirs{j}),ThetaLock_2Session.(sessDirs{j}).thetaLock_delay.ppcT(rateLabel)];
            rT.(sessDirs{j})= [rT.(sessDirs{j}),ThetaLock_2Session.(sessDirs{j}).thetaLock_delay.rT(rateLabel)];
            angT.(sessDirs{j})= [angT.(sessDirs{j}),ThetaLock_2Session.(sessDirs{j}).thetaLock_delay.angT(rateLabel)];
            numT.(sessDirs{j})= [numT.(sessDirs{j}),ThetaLock_2Session.(sessDirs{j}).thetaLock_delay.numT(rateLabel)];
            pT.(sessDirs{j})= [pT.(sessDirs{j}),ThetaLock_2Session.(sessDirs{j}).thetaLock_delay.pT(rateLabel)];
        end
    end
end

validInd_assembly_on10 =  assembly_numT.on10 >= p.numT;
validInd_assembly_on30 =  assembly_numT.on30 >= p.numT;
validInd_on10 =  numT.on10 >= p.numT;
validInd_on30 =  numT.on30 >= p.numT;

% correct the angle to remove the effect from circMean
assembly_angT.on10(assembly_angT.on10<0) = assembly_angT.on10(assembly_angT.on10<0) + 360;
assembly_angT.on30(assembly_angT.on30<0) = assembly_angT.on30(assembly_angT.on30<0) + 360;
angT.on10(angT.on10<0) = angT.on10(angT.on10<0) + 360;
angT.on30(angT.on30<0) = angT.on30(angT.on30<0) + 360;


figure
subplot(1,3,1)
bar(1,100*sum(assembly_pT.on10(validInd_assembly_on10)<=0.05)/sum(validInd_assembly_on10));
hold on
bar(2,100*sum(assembly_pT.on30(validInd_assembly_on30)<=0.05)/sum(validInd_assembly_on30));
bar(4,100*sum(pT.on10(validInd_on10)<=0.05)/sum(validInd_on10));
bar(5,100*sum(pT.on30(validInd_on30)<=0.05)/sum(validInd_on30));
ylim([0 100])
title('lock %')
set(gca, 'XTick', [1,2,4,5],'XTickLabel',{'Assembly on10','on30','Single on10','on30'});


subplot(1,3,2)
Violin(assembly_ppcT.on10(validInd_assembly_on10),1,'ShowData',false);
hold on
Violin(assembly_ppcT.on30(validInd_assembly_on30),2,'ShowData',false);
Violin(ppcT.on10(validInd_on10),4,'ShowData',false);
Violin(ppcT.on30(validInd_on30),5,'ShowData',false);
ylim([0 0.5])
title('ppc')

subplot(1,3,3)
Violin(assembly_rT.on10(validInd_assembly_on10),1,'ShowData',false);
hold on
Violin(assembly_rT.on30(validInd_assembly_on30),2,'ShowData',false);
Violin(rT.on10(validInd_on10),4,'ShowData',false);
Violin(rT.on30(validInd_on30),5,'ShowData',false);
ylim([0 0.8])
title('resultant length')

figure
subplot(2,4,1)
% polar_lim(spikeRadPhaseT,ones(length(ts),1)*MaxT,'o')
cellInd = (validInd_assembly_on10 & assembly_pT.on10<=0.05);
Num = sum(cellInd);
[tout,rout]=rose(circ_ang2rad(assembly_angT.on10(cellInd)),360/30);
polar_lim(tout,rout./Num,MaxT);
title('Assembly on10')

subplot(2,4,5)
x=15:30:345;
nT_asmb_on10 = hist(assembly_angT.on10(cellInd),x)/Num;
bar([x x+360],[nT_asmb_on10 nT_asmb_on10],1,'r')
xlabel('Theta phase (deg)'); ylabel('Normalized counts')
ylim([0 0.5])
set(gca, 'XTick', [0,180,360,540,720]);
xlim([0 720])

subplot(2,4,2)
% polar_lim(spikeRadPhaseT,ones(length(ts),1)*MaxT,'o')
cellInd = (validInd_assembly_on30 & assembly_pT.on30<=0.05);
Num = sum(cellInd);
[tout,rout]=rose(circ_ang2rad(assembly_angT.on30(cellInd)),360/30);
polar_lim(tout,rout./Num,MaxT);
title('Assembly on30')

subplot(2,4,6)
x=15:30:345;
nT_asmb_on30 = hist(assembly_angT.on30(cellInd),x)/Num;
bar([x x+360],[nT_asmb_on30 nT_asmb_on30],1,'r')
xlabel('Theta phase (deg)'); ylabel('Normalized counts')
ylim([0 0.5])
set(gca, 'XTick', [0,180,360,540,720]);
xlim([0 720])

subplot(2,4,3)
% polar_lim(spikeRadPhaseT,ones(length(ts),1)*MaxT,'o')
cellInd = (validInd_on10 & pT.on10<=0.05);
Num = sum(cellInd);
[tout,rout]=rose(circ_ang2rad(angT.on10(cellInd)),360/30);
polar_lim(tout,rout./Num,MaxT);
title('single cell on10')

subplot(2,4,7)
x=15:30:345;
nT_on10 = hist(angT.on10(cellInd),x)/Num;
bar([x x+360],[nT_on10 nT_on10],1,'r')
xlabel('Theta phase (deg)'); ylabel('Normalized counts')
ylim([0 0.5])
set(gca, 'XTick', [0,180,360,540,720]);
xlim([0 720])

subplot(2,4,4)
% polar_lim(spikeRadPhaseT,ones(length(ts),1)*MaxT,'o')
cellInd = (validInd_on30 & pT.on30<=0.05);
Num = sum(cellInd);
[tout,rout]=rose(circ_ang2rad(angT.on30(cellInd)),360/30);
polar_lim(tout,rout./Num,MaxT);
title('single cell on30')

subplot(2,4,8)
x=15:30:345;
nT_on30 = hist(angT.on30(cellInd),x)/Num;
bar([x x+360],[nT_on30 nT_on30],1,'r')
xlabel('Theta phase (deg)'); ylabel('Normalized counts')
ylim([0 0.5])
set(gca, 'XTick', [0,180,360,540,720]);
xlim([0 720])

[h1,p1] = kstest2(nT_asmb_on10,nT_on10)
[h2,p2] = kstest2(nT_asmb_on30,nT_on30)
end