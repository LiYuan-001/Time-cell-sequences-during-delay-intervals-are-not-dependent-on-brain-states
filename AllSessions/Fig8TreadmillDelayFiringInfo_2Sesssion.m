function Fig8TreadmillDelayFiringInfo_2Sesssion(inFile,AnalyzeSes)
close all

p.savePlot = 1;
p.writeToFile = 1;

% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\DelayInfoStat');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end
    
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    % load delayFire map
    delayFile = fullfile(sessInfo(i).mainDir,'Cell Property','Fig8DelayTimeMap_2Session.mat');
    load(delayFile);
    
    % initiate the data
    DelayInfoStat_2Session.rat = sessInfo(i).animal;
    DelayInfoStat_2Session.day = sessInfo(i).day;
    
    TList = Fig8DelayTimeMap_2Session.tList;
    clusterNum = length(Fig8DelayTimeMap_2Session.tList);
    DelayInfoStat_2Session.tList = Fig8DelayTimeMap_2Session.tList;
            
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    
    % get each phase names (no delay etc)
    sessDirs = {'on10','off10','on30','off30'};
    

    for j = 1:length(sessDirs)
        avgRate = nan(1,clusterNum);
        peakRate = nan(1,clusterNum);
        info_spike = nan(1,clusterNum);
        info_sec = nan(1,clusterNum);
        sparsity = nan(1,clusterNum);
        selectivity = nan(1,clusterNum);
        zCoherence = nan(1,clusterNum);
        % -----------------------------------------------------------------
        for k = 1:clusterNum
            
            spikeRate1_Smooth = Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate1_Smooth{k};
            spikeRate2_Smooth = Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate2_Smooth{k};
            
            spikeRate1_Combined_Smooth = Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate1_Combined_Smooth{k};
            spikeRate2_Combined_Smooth = Fig8DelayTimeMap_2Session.(sessDirs{j}).spikeRate2_Combined_Smooth{k};

            % calculate information/spike, info/sec, coherence, sparseness
            avgRate(k) = mean(spikeRate1_Combined_Smooth);
            peakRate(k) = max(spikeRate1_Combined_Smooth);
            
            posPDF = ones(1,length(spikeRate1_Combined_Smooth))*(1/length(spikeRate1_Combined_Smooth));
            [info_spike(k),info_sec(k),sparsity(k),selectivity(k)] = mapstat(spikeRate1_Combined_Smooth,posPDF);
            zCoherence(k) = fieldcohere(spikeRate1_Combined_Smooth); 
            if j == 1                
                h = figure(k);
                h.Position = [100 100 1600 800];
            else
                figure(k)
            end
            
            subplot(1,4,j)
            imagesc([0,size(spikeRate1_Smooth,2)],[9 8+size(spikeRate1_Smooth,1)],spikeRate1_Smooth/max(max(spikeRate1_Smooth)))
            hold on
            imagesc([0,size(spikeRate1_Combined_Smooth,2)],7,spikeRate1_Combined_Smooth/max(spikeRate1_Combined_Smooth))
            set(gca,'YDir','normal')
            colormap(jet)
            if j == 1
                TITLE1 = sprintf('%s%d%s%d%s%s%s%s','DelayZone2-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k},'-',sessDirs{j});
            else
                TITLE1 = sessDirs{j};
            end
            title({TITLE1},'Interpreter','None')
            axis off
            text1 = sprintf('%s%2.2f%s','Peak rate : ',peakRate(k),' Hz');
            text2 = sprintf('%s%2.2f%s%2.2f','Info/spike: ',info_spike(k),'  Info/sec: ',info_sec(k));
            text3 = sprintf('%s%2.2f%s%2.2f','Sparsity: ',sparsity(k),'  Selectivity: ',selectivity(k));
            text4 = sprintf('%s%2.2f','zCoherence: ',zCoherence(k));
            text([1,1,1,1],[5,4,3,2],{text1,text2,text3,text4},'FontSize',10,'Interpreter','None')
            ylim([0 20])
        end 
        DelayInfoStat_2Session.(sessDirs{j}).avgRate = avgRate;        
        DelayInfoStat_2Session.(sessDirs{j}).peakRate = peakRate;
        DelayInfoStat_2Session.(sessDirs{j}).info_spike = info_spike;        
        DelayInfoStat_2Session.(sessDirs{j}).info_sec = info_sec;
        DelayInfoStat_2Session.(sessDirs{j}).sparsity = sparsity;        
        DelayInfoStat_2Session.(sessDirs{j}).selectivity = selectivity;
        DelayInfoStat_2Session.(sessDirs{j}).zCoherence = zCoherence;
    end
    
    if p.savePlot == 1
        for k = 1:clusterNum
            figure(k)
            figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',TList{k},'-DelayInfoStat');
            print(figName,'-dpng','-r300');
        end
    end
    close all
    if p.writeToFile == 1
        save(fullfile(savedir2,'DelayInfoStat_2Session.mat'), 'DelayInfoStat_2Session');
    end
    clear DelayInfoStat
    fprintf('Finished analysis for session %d\n',i);
end

end

% Shannon information, sparseness, and selectivity  
function [information,information2,sparsity,selectivity] = mapstat(map,posPDF)

% Sparseness
n = size(map,1);
meanrate = nansum(nansum( map .* posPDF ));
meansquarerate = nansum(nansum( (map.^2) .* posPDF ));
if meansquarerate == 0
    sparsity = NaN;
else
    sparsity = meanrate^2 / meansquarerate;
end

% Selectivity
maxrate = max(max(map));
if meanrate == 0
   selectivity = NaN;
else
   selectivity = maxrate/meanrate;
end

% Shannon information
% "information density (bits/spike)"
[i1, i2] = find( (map>0) & (posPDF>0) );  % the limit of x*log(x) as x->0 is 0 
if length(i1)>0
    akksum = 0;
    for i = 1:length(i1)
        ii1 = i1(i);
        ii2 = i2(i);
        % information density (bits/spike)
        akksum = akksum + posPDF(ii1,ii2) * (map(ii1,ii2)/meanrate) * log2( map(ii1,ii2) / meanrate ); 
    end
    information = akksum;
else
    information = NaN;
end


% Shannon information
% "information rate (bits/sec)":
if length(find(posPDF>0))<1
    information2 = NaN;
else
    [i1, i2] = find( (map>0) & (posPDF>0) );  % the limit of x*log(x) as x->0 is 0 
    if ~isempty(i1)
        akksum = 0;
        for i = 1:length(i1)
            ii1 = i1(i);
            ii2 = i2(i);
            % information density (bits/spike)
            akksum = akksum + posPDF(ii1,ii2) * map(ii1,ii2) * log2( map(ii1,ii2) / meanrate );
        end
        information2 = akksum;
    else
        information2 = 0;
    end
end
end

function z = fieldcohere(map)
[n,m] = size(map);
tmp = zeros(n*m,2);
k=0;
for y = 1:n
    for x = 1:m
        k = k + 1;
        xstart = max([1,x-1]);
        ystart = max([1,y-1]);
        xend = min([m x+1]);
        yend = min([n y+1]);
        nn = sum(sum(isfinite(map(ystart:yend,xstart:xend)))) - isfinite(map(y,x));
        if (nn > 0)
            tmp(k,1) = map(y,x);
            tmp(k,2) = nansum([ nansum(nansum(map(ystart:yend,xstart:xend))) , -map(y,x) ]) / nn;
        else
            tmp(k,:) = [NaN,NaN];    
        end
    end
end
index = find( isfinite(tmp(:,1)) & isfinite(tmp(:,2)) );
if length(index) > 3
    cc = corrcoef(tmp(index,:));
    z = atanh(cc(2,1));
else
    z = NaN;
end
end