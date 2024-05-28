% calculate cell assemblies throught the sleep-Fig8-sleep session
% treat it as a unit and plot rate map
% Li Yuan, UCSD, Aug-30-2022
% 
function Fig8Treadmill_PCAICA_Delay_LR_Quant(inFile,AnalyzeSes)

close all
p.savePlot = 0;
p.writeToFile = 1;


% % Read in input information
sessInfo = SessInfoImport(inFile);
strengthNames = {'on_L','on_R','off_L','off_R'};
for n = 1:length(strengthNames)
    on_L_strength.(strengthNames{n}) = [];
    on_R_strength.(strengthNames{n}) = [];
    off_L_strength.(strengthNames{n}) = [];
    off_R_strength.(strengthNames{n}) = [];
    
    on_L_patSim.(strengthNames{n}) = [];
    on_R_patSim.(strengthNames{n}) = [];
    off_L_patSim.(strengthNames{n}) = [];
    off_R_patSim.(strengthNames{n}) = [];
    
    on_L_commonInd.(strengthNames{n}) = [];
    on_R_commonInd.(strengthNames{n}) = [];
    off_L_commonInd.(strengthNames{n}) = [];
    off_R_commonInd.(strengthNames{n}) = [];
end

on_strength.on = [];
on_strength.off = [];
off_strength.on = [];
off_strength.off = [];

on_patSim.on = [];
on_patSim.off = [];
off_patSim.on = [];
off_patSim.off = [];

on_commonInd.on = [];
on_commonInd.off = [];
off_commonInd.on = [];
off_commonInd.off = [];

on_cellNum = [];
off_cellNum = [];

for i = AnalyzeSes(1:end)
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');    
    % load assembly
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_DelayOnOff_LR-100ms.mat');
    load(assemblyFile);
    
    binTime = CellAssembly_DelayLR.DelayOn_L.binTime;
    
    bin_on_L_ind = zeros(1,length(binTime));
    bin_on_R_ind = zeros(1,length(binTime));
    bin_off_L_ind = zeros(1,length(binTime));
    bin_off_R_ind = zeros(1,length(binTime));
    
    sessDirs = {'on10_1','on10_2','on30_1','on30_2'};    
    for j = 1:length(sessDirs) 
        % delay area spike and bin for assembly detection
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        map_1D = load(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'));

        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        % delayTend1 = Fig8DelayZonePos.delayPos1.endT;        
        trialNum = size(delayTstart1,2);
        
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
        else
            error('Delay time is wrong')
        end       
        
        for m = 1:trialNum

            startT = delayTstart1(m);
            endT = delayTend1_2(m);
            
            if map_1D.ratesByECLR.ECLR(m) == 1 || map_1D.ratesByECLR.ECLR(m) == 2
%                 posLabel = 2;
                turnLabel = 1;
                ind = (binTime>=startT & binTime<endT);
                bin_on_L_ind = bin_on_L_ind + ind;
            elseif map_1D.ratesByECLR.ECLR(m) == 3 || map_1D.ratesByECLR.ECLR(m) == 4
%                 posLabel = 1;
                turnLabel = 0;
                ind = (binTime>=startT & binTime<endT);
                bin_on_R_ind = bin_on_R_ind + ind;
            else
                error('Turning label ERROR');
            end
        end        
    end
    
    sessDirs = {'off10_1','off10_2','off30_1','off30_2'};    
    for j = 1:length(sessDirs) 
        % delay area spike and bin for assembly detection
        % load analyzed positions
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        map_1D = load(fullfile(sessInfo(i).mainDir,sessDirs{j},'ratesByECLR.mat'));

        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        % delayTend1 = Fig8DelayZonePos.delayPos1.endT;        
        trialNum = size(delayTstart1,2);
        
        if contains(sessDirs{j},'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
        elseif contains(sessDirs{j},'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
        else
            error('Delay time is wrong')
        end       
        
        for m = 1:trialNum

            startT = delayTstart1(m);
            endT = delayTend1_2(m);
            
            if map_1D.ratesByECLR.ECLR(m) == 1 || map_1D.ratesByECLR.ECLR(m) == 2
%                 posLabel = 2;
                turnLabel = 1;
                ind = (binTime>=startT & binTime<endT);
                bin_off_L_ind = bin_off_L_ind + ind;
            elseif map_1D.ratesByECLR.ECLR(m) == 3 || map_1D.ratesByECLR.ECLR(m) == 4
%                 posLabel = 1;
                turnLabel = 0;
                ind = (binTime>=startT & binTime<endT);
                bin_off_R_ind = bin_off_R_ind + ind;
            else
                error('Turning label ERROR');
            end
        end        
    end
    
    for k = 1:CellAssembly_DelayLR.DelayOn_L.patNum
        on_L_strength.on_L = [on_L_strength.on_L,sum(CellAssembly_DelayLR.DelayOn_L.AssmblStrength(bin_on_L_ind>0)>=5)];
        on_L_strength.on_R = [on_L_strength.on_R,sum(CellAssembly_DelayLR.DelayOn_L.AssmblStrength(bin_on_R_ind>0)>=5)];
        on_L_strength.off_L = [on_L_strength.off_L,sum(CellAssembly_DelayLR.DelayOn_L.AssmblStrength(bin_off_L_ind>0)>=5)];
        on_L_strength.off_R = [on_L_strength.off_R,sum(CellAssembly_DelayLR.DelayOn_L.AssmblStrength(bin_off_R_ind>0)>=5)];
        
        simVal = zeros(length(CellAssembly_DelayLR.DelayOn_R.patNum),1);
        for m = 1:CellAssembly_DelayLR.DelayOn_R.patNum
            pattern_1 = CellAssembly_DelayLR.DelayOn_L.AssmblWght(:,k);
            pattern_2 = CellAssembly_DelayLR.DelayOn_R.AssmblWght(:,m);
            simVal(m) = dot(pattern_1,pattern_2)/(norm(pattern_1)*norm(pattern_2)); 
        end
        on_L_patSim.on_R = [on_L_patSim.on_R,max(simVal)];
        on_L_commonInd.on_R = [on_L_commonInd.on_R,max(simVal)>=0.5];
        
        simVal = zeros(length(CellAssembly_DelayLR.DelayOff_L.patNum),1);
        for m = 1:CellAssembly_DelayLR.DelayOff_L.patNum
            pattern_1 = CellAssembly_DelayLR.DelayOn_L.AssmblWght(:,k);
            pattern_2 = CellAssembly_DelayLR.DelayOff_L.AssmblWght(:,m);
            simVal(m) = dot(pattern_1,pattern_2)/(norm(pattern_1)*norm(pattern_2)); 
        end
        on_L_patSim.off_L = [on_L_patSim.off_L,max(simVal)];
        on_L_commonInd.off_L = [on_L_commonInd.off_L,max(simVal)>=0.5];
        
        
        simVal = zeros(length(CellAssembly_DelayLR.DelayOff_R.patNum),1);
        for m = 1:CellAssembly_DelayLR.DelayOff_R.patNum
            pattern_1 = CellAssembly_DelayLR.DelayOn_L.AssmblWght(:,k);
            pattern_2 = CellAssembly_DelayLR.DelayOff_R.AssmblWght(:,m);
            simVal(m) = dot(pattern_1,pattern_2)/(norm(pattern_1)*norm(pattern_2)); 
        end
        on_L_patSim.off_R = [on_L_patSim.off_R,max(simVal)];
        on_L_commonInd.off_R = [on_L_commonInd.off_R,max(simVal)>=0.5];
    end
    
    for k = 1:CellAssembly_DelayLR.DelayOn_R.patNum
        on_R_strength.on_L = [on_R_strength.on_L,sum(CellAssembly_DelayLR.DelayOn_R.AssmblStrength(bin_on_L_ind>0)>=5)];
        on_R_strength.on_R = [on_R_strength.on_R,sum(CellAssembly_DelayLR.DelayOn_R.AssmblStrength(bin_on_R_ind>0)>=5)];
        on_R_strength.off_L = [on_R_strength.off_L,sum(CellAssembly_DelayLR.DelayOn_R.AssmblStrength(bin_off_L_ind>0)>=5)];
        on_R_strength.off_R = [on_R_strength.off_R,sum(CellAssembly_DelayLR.DelayOn_R.AssmblStrength(bin_off_R_ind>0)>=5)];
        
        
        simVal = zeros(length(CellAssembly_DelayLR.DelayOn_L.patNum),1);
        for m = 1:CellAssembly_DelayLR.DelayOn_L.patNum
            pattern_1 = CellAssembly_DelayLR.DelayOn_R.AssmblWght(:,k);
            pattern_2 = CellAssembly_DelayLR.DelayOn_L.AssmblWght(:,m);
            simVal(m) = dot(pattern_1,pattern_2)/(norm(pattern_1)*norm(pattern_2)); 
        end
        on_R_patSim.on_L = [on_R_patSim.on_L,max(simVal)];
        on_R_commonInd.on_L = [on_R_commonInd.on_L,max(simVal)>=0.5];
        
        
        simVal = zeros(length(CellAssembly_DelayLR.DelayOff_L.patNum),1);
        for m = 1:CellAssembly_DelayLR.DelayOff_L.patNum
            pattern_1 = CellAssembly_DelayLR.DelayOn_R.AssmblWght(:,k);
            pattern_2 = CellAssembly_DelayLR.DelayOff_L.AssmblWght(:,m);
            simVal(m) = dot(pattern_1,pattern_2)/(norm(pattern_1)*norm(pattern_2)); 
        end
        on_R_patSim.off_L = [on_R_patSim.off_L,max(simVal)];
        on_R_commonInd.off_L = [on_R_commonInd.off_L,max(simVal)>=0.5];
        
        
        simVal = zeros(length(CellAssembly_DelayLR.DelayOff_R.patNum),1);
        for m = 1:CellAssembly_DelayLR.DelayOff_R.patNum
            pattern_1 = CellAssembly_DelayLR.DelayOn_R.AssmblWght(:,k);
            pattern_2 = CellAssembly_DelayLR.DelayOff_R.AssmblWght(:,m);
            simVal(m) = dot(pattern_1,pattern_2)/(norm(pattern_1)*norm(pattern_2)); 
        end
        on_R_patSim.off_R = [on_R_patSim.off_R,max(simVal)];
        on_R_commonInd.off_R = [on_R_commonInd.off_R,max(simVal)>=0.5];
    end
    
    for k = 1:CellAssembly_DelayLR.DelayOff_L.patNum
        off_L_strength.on_L = [off_L_strength.on_L,mean(CellAssembly_DelayLR.DelayOff_L.AssmblStrength(bin_on_L_ind>0))];
        off_L_strength.on_R = [off_L_strength.on_R,mean(CellAssembly_DelayLR.DelayOff_L.AssmblStrength(bin_on_R_ind>0))];
        off_L_strength.off_L = [off_L_strength.off_L,mean(CellAssembly_DelayLR.DelayOff_L.AssmblStrength(bin_off_L_ind>0))];
        off_L_strength.off_R = [off_L_strength.off_R,mean(CellAssembly_DelayLR.DelayOff_L.AssmblStrength(bin_off_R_ind>0))];
        
        simVal = zeros(length(CellAssembly_DelayLR.DelayOn_L.patNum),1);
        for m = 1:CellAssembly_DelayLR.DelayOn_L.patNum
            pattern_1 = CellAssembly_DelayLR.DelayOff_L.AssmblWght(:,k);
            pattern_2 = CellAssembly_DelayLR.DelayOn_L.AssmblWght(:,m);
            simVal(m) = dot(pattern_1,pattern_2)/(norm(pattern_1)*norm(pattern_2)); 
        end
        off_L_patSim.on_L = [off_L_patSim.on_L,max(simVal)];
        off_L_commonInd.on_L = [off_L_commonInd.on_L,max(simVal)>=0.5];
        
        
        simVal = zeros(length(CellAssembly_DelayLR.DelayOn_R.patNum),1);
        for m = 1:CellAssembly_DelayLR.DelayOn_R.patNum
            pattern_1 = CellAssembly_DelayLR.DelayOff_L.AssmblWght(:,k);
            pattern_2 = CellAssembly_DelayLR.DelayOn_R.AssmblWght(:,m);
            simVal(m) = dot(pattern_1,pattern_2)/(norm(pattern_1)*norm(pattern_2)); 
        end
        off_L_patSim.on_R = [off_L_patSim.on_R,max(simVal)];
        off_L_commonInd.on_R = [off_L_commonInd.on_R,max(simVal)>=0.5];
        
        simVal = zeros(length(CellAssembly_DelayLR.DelayOff_R.patNum),1);
        for m = 1:CellAssembly_DelayLR.DelayOff_R.patNum
            pattern_1 = CellAssembly_DelayLR.DelayOff_L.AssmblWght(:,k);
            pattern_2 = CellAssembly_DelayLR.DelayOff_R.AssmblWght(:,m);
            simVal(m) = dot(pattern_1,pattern_2)/(norm(pattern_1)*norm(pattern_2)); 
        end
        off_L_patSim.off_R = [off_L_patSim.off_R,max(simVal)];
        off_L_commonInd.off_R = [off_L_commonInd.off_R,max(simVal)>=0.5];
    end
    
    for k = 1:CellAssembly_DelayLR.DelayOff_R.patNum
        off_R_strength.on_L = [off_R_strength.on_L,mean(CellAssembly_DelayLR.DelayOff_R.AssmblStrength(bin_on_L_ind>0))];
        off_R_strength.on_R = [off_R_strength.on_R,mean(CellAssembly_DelayLR.DelayOff_R.AssmblStrength(bin_on_R_ind>0))];
        off_R_strength.off_L = [off_R_strength.off_L,mean(CellAssembly_DelayLR.DelayOff_R.AssmblStrength(bin_off_L_ind>0))];
        off_R_strength.off_R = [off_R_strength.off_R,mean(CellAssembly_DelayLR.DelayOff_R.AssmblStrength(bin_off_R_ind>0))];
        
        simVal = zeros(length(CellAssembly_DelayLR.DelayOn_L.patNum),1);
        for m = 1:CellAssembly_DelayLR.DelayOn_L.patNum
            pattern_1 = CellAssembly_DelayLR.DelayOff_R.AssmblWght(:,k);
            pattern_2 = CellAssembly_DelayLR.DelayOn_L.AssmblWght(:,m);
            simVal(m) = dot(pattern_1,pattern_2)/(norm(pattern_1)*norm(pattern_2)); 
        end
        off_R_patSim.on_L = [off_R_patSim.on_L,max(simVal)];
        off_R_commonInd.on_L = [off_R_commonInd.on_L,max(simVal)>=0.5];
        
        
        simVal = zeros(length(CellAssembly_DelayLR.DelayOn_R.patNum),1);
        for m = 1:CellAssembly_DelayLR.DelayOn_R.patNum
            pattern_1 = CellAssembly_DelayLR.DelayOff_R.AssmblWght(:,k);
            pattern_2 = CellAssembly_DelayLR.DelayOn_R.AssmblWght(:,m);
            simVal(m) = dot(pattern_1,pattern_2)/(norm(pattern_1)*norm(pattern_2)); 
        end
        off_R_patSim.on_R = [off_R_patSim.on_R,max(simVal)];
        off_R_commonInd.on_R = [off_R_commonInd.on_R,max(simVal)>=0.5];
        
        simVal = zeros(length(CellAssembly_DelayLR.DelayOff_L.patNum),1);
        for m = 1:CellAssembly_DelayLR.DelayOff_L.patNum
            pattern_1 = CellAssembly_DelayLR.DelayOff_R.AssmblWght(:,k);
            pattern_2 = CellAssembly_DelayLR.DelayOff_L.AssmblWght(:,m);
            simVal(m) = dot(pattern_1,pattern_2)/(norm(pattern_1)*norm(pattern_2)); 
        end
        off_R_patSim.off_L = [off_R_patSim.off_L,max(simVal)];
        off_R_commonInd.off_L = [off_R_commonInd.off_L,max(simVal)>=0.5];
    end
   
    for k = 1:CellAssembly_DelayLR.DelayOn.patNum
        on_strength.on = [on_strength.on,mean(CellAssembly_DelayLR.DelayOn.AssmblStrength((bin_on_L_ind+bin_on_R_ind)>0))];
        on_strength.off = [on_strength.off,mean(CellAssembly_DelayLR.DelayOn.AssmblStrength((bin_off_L_ind+bin_off_R_ind)>0))];
        on_cellNum = [on_cellNum,length(CellAssembly_DelayLR.DelayOn.AssmblPtrnCellIDs{k})];
        
        simVal = zeros(length(CellAssembly_DelayLR.DelayOff.patNum),1);
        for m = 1:CellAssembly_DelayLR.DelayOff.patNum
            pattern_1 = CellAssembly_DelayLR.DelayOn.AssmblWght(:,k);
            pattern_2 = CellAssembly_DelayLR.DelayOff.AssmblWght(:,m);
            simVal(m) = dot(pattern_1,pattern_2)/(norm(pattern_1)*norm(pattern_2)); 
        end
        on_patSim.off = [on_patSim.off,max(simVal)];
        on_commonInd.off = [on_commonInd.off,max(simVal)>=0.5];
    end
    
    for k = 1:CellAssembly_DelayLR.DelayOff.patNum
        off_strength.on = [off_strength.on,mean(CellAssembly_DelayLR.DelayOff.AssmblStrength((bin_on_L_ind+bin_on_R_ind)>0))];
        off_strength.off = [off_strength.off,mean(CellAssembly_DelayLR.DelayOff.AssmblStrength((bin_off_L_ind+bin_off_R_ind)>0))];
        off_cellNum = [off_cellNum,length(CellAssembly_DelayLR.DelayOff.AssmblPtrnCellIDs{k})];
        
        simVal = zeros(length(CellAssembly_DelayLR.DelayOn.patNum),1);
        for m = 1:CellAssembly_DelayLR.DelayOn.patNum
            pattern_1 = CellAssembly_DelayLR.DelayOff.AssmblWght(:,k);
            pattern_2 = CellAssembly_DelayLR.DelayOn.AssmblWght(:,m);
            simVal(m) = dot(pattern_1,pattern_2)/(norm(pattern_1)*norm(pattern_2)); 
        end
        off_patSim.on = [off_patSim.on,max(simVal)];
        off_commonInd.on = [off_commonInd.on,max(simVal)>=0.5];
    end
end

figure
Violin(on_patSim.off,1);
Violin(off_patSim.on,2);

figure
boxplot(on_cellNum,'Position',1);
hold on
boxplot(off_cellNum,'Position',2);
xlim([0 3])

figure
subplot(2,1,1)
% Violin(on_L_strength.on_L,1);
Violin(on_L_patSim.on_R,1);
Violin(on_R_patSim.on_L,2);
xlim([0 3])

subplot(2,1,2)
Violin(off_L_patSim.off_R,1);
Violin(off_R_patSim.off_L,2);
xlim([0 3])

figure
subplot(4,1,3)
Violin(off_L_patSim.on_L,1);
Violin(off_L_patSim.on_R,2);
% Violin(Off_R_patSim.off_L,3);
Violin(off_L_patSim.off_R,4);
xlim([0 5])

subplot(4,1,4)
Violin(off_R_patSim.on_L,1);
Violin(off_R_patSim.on_R,2);
Violin(off_R_patSim.off_L,3);
% Violin(Off_R_patSim.off_R,4);
xlim([0 5])


figure
subplot(4,2,1)
% Violin(on_L_strength.on_L(on_L_commonInd.on_R==0),1);
% Violin(on_L_strength.on_R(on_L_commonInd.on_R==0),2);
% Violin(on_L_strength.on_L(on_L_commonInd.on_R==1),3);
% Violin(on_L_strength.on_R(on_L_commonInd.on_R==1),4);

plot([on_L_strength.on_L(on_L_commonInd.on_R==0);on_L_strength.on_R(on_L_commonInd.on_R==0)],'o-');
subplot(4,2,2)
plot([on_L_strength.on_L(on_L_commonInd.on_R==1);on_L_strength.on_R(on_L_commonInd.on_R==1)],'o-');

% subplot(4,1,2)
% Violin(on_R_strength.on_L(on_R_commonInd.on_L==0),1);
% Violin(on_R_strength.on_R(on_R_commonInd.on_L==0),2);
% Violin(on_R_strength.on_L(on_R_commonInd.on_L==1),3);
% Violin(on_R_strength.on_R(on_R_commonInd.on_L==1),4);

% Violin(on_R_strength.on_L,1);
% Violin(on_R_strength.on_R,2);
% Violin(on_R_strength.off_L,3);
% Violin(on_R_strength.off_R,4);

subplot(4,2,3)
% Violin(off_L_strength.on_L,1);
% Violin(off_L_strength.on_R,2);
% Violin(off_L_strength.off_L,3);
% Violin(off_L_strength.off_R,4);
plot([on_R_strength.on_L(on_R_commonInd.on_L==0);on_R_strength.on_R(on_R_commonInd.on_L==0)],'o-');
subplot(4,2,4)
plot([on_R_strength.on_L(on_R_commonInd.on_L==1);on_R_strength.on_R(on_R_commonInd.on_L==1)],'o-');

% subplot(4,1,4)
% Violin(off_R_strength.on_L,1);
% Violin(off_R_strength.on_R,2);
% Violin(off_R_strength.off_L,3);
% Violin(off_R_strength.off_R,4);

figure
subplot(4,1,1)
boxplot(on_L_strength.on_L,'Position',1);
hold on
boxplot(on_L_strength.on_R,'Position',2);
boxplot(on_L_strength.off_L,'Position',3);
boxplot(on_L_strength.off_R,'Position',4);
xlim([0 5])

subplot(4,1,2)
boxplot(on_R_strength.on_L,'Position',1);
hold on
boxplot(on_R_strength.on_R,'Position',2);
boxplot(on_R_strength.off_L,'Position',3);
boxplot(on_R_strength.off_R,'Position',4);
xlim([0 5])

subplot(4,1,3)
boxplot(off_L_strength.on_L,'Position',1);
hold on
boxplot(off_L_strength.on_R,'Position',2);
boxplot(off_L_strength.off_L,'Position',3);
boxplot(off_L_strength.off_R,'Position',4);
xlim([0 5])

subplot(4,1,4)
boxplot(off_R_strength.on_L,'Position',1);
hold on
boxplot(off_R_strength.on_R,'Position',2);
boxplot(off_R_strength.off_L,'Position',3);
boxplot(off_R_strength.off_R,'Position',4);
xlim([0 5])
end