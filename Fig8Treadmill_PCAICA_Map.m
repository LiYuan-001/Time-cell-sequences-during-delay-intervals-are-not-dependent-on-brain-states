% calculate cell assemblies throught the sleep-Fig8-sleep session
% treat it as a unit and plot rate map
% Li Yuan, UCSD, May-02-2022
% 
function Fig8Treadmill_PCAICA_Map(inFile,AnalyzeSes)

close all
p.savePlot = 1;
p.writeToFile = 1;

p.binWidth = 5;
p.smoothing = 5;

p.savePlot = 1;
p.writeToFile = 1;


% % Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    sessDirs = sessInfo(i).sessDirs;
    sleepDirs = sessInfo(i).sleepDirs;
    
    if p.savePlot
        % Li Yuan, UCSD, 15-Nov-2019
        
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\Assembly spatial map');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load assembly file
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_WholeSes.mat');
    load(assemblyFile);
    
    % initiate the data
    CellAssembly_Map.rat = sessInfo(i).animal;
    CellAssembly_Map.day = sessInfo(i).day;
    CellAssembly_Map.timeBin = CellAssembly_WholeSes.binWidth;
    
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;    
    tList = SpikeProp.tList(rateLabel);
    
    patNum = CellAssembly_WholeSes.patNum ;
    AssmblPtrnCellIDs = CellAssembly_WholeSes.AssmblPtrnCellIDs;
    AssmblWght = CellAssembly_WholeSes.AssmblWght;
    AssmblStrength = CellAssembly_WholeSes.AssmblStrength;
    event_Time = CellAssembly_WholeSes.event_Time;
    event_strength = CellAssembly_WholeSes.event_strength;
    event_Num = CellAssembly_WholeSes.event_Num;
    
    %% fig 8 maze
    [~, pathData, pathDataIdeal, pathDataLin] = loadPathInfo_Treadmill(sessInfo(i));
    for j = 1:length(sessDirs)
        
        %             tinfo = trialInfo.(SessDirs{j});
        %             pinfo = parsingInfo.(SessDirs{j});
        pdata = pathData.(sessDirs{j});
        pIdeal = pathDataIdeal.(sessDirs{j});
        plin = pathDataLin.(sessDirs{j});
        
        % calculate original path rate map
        mapAxis1 = min(pdata.x)-p.binWidth/2:p.binWidth:max(pdata.x)+p.binWidth/2;
        mapAxis2 = min(pdata.y)-p.binWidth/2:p.binWidth:max(pdata.y)+p.binWidth/2;
        
        mapAxisIdeal1 = min(pIdeal.x)-p.binWidth/2:p.binWidth:max(pIdeal.x)+p.binWidth/2;
        mapAxisIdeal2 = min(pIdeal.y)-p.binWidth/2:p.binWidth:max(pIdeal.y)+p.binWidth/2;
        
        % Calulate what areas of the box that have been visited
        visited = visitedBins2(pdata.x,pdata.y,mapAxis1,mapAxis2);        
        visitedIdeal = visitedBins2(pIdeal.x,pIdeal.y,mapAxisIdeal1,mapAxisIdeal2);
        
        % for .mat files which is going to be write down
        CellAssembly_Map.(sessDirs{j}).mapAxis1 = mapAxis1;
        CellAssembly_Map.(sessDirs{j}).mapAxis2 = mapAxis2;
        CellAssembly_Map.(sessDirs{j}).mapAxisIdeal1 = mapAxisIdeal1;
        CellAssembly_Map.(sessDirs{j}).mapAxisIdeal2 = mapAxisIdeal2;
        CellAssembly_Map.(sessDirs{j}).visited = visited;
        CellAssembly_Map.(sessDirs{j}).visitedIdeal = visitedIdeal;
        
        % start calculate rate map & plot
        avgRate = zeros(patNum,1);
        peakRate = zeros(patNum,1);
        
        for k = 1:patNum            
            h=figure(k);
            if j ==1
                h.Position = [100 100 1600 900];
                subplot(4,length(sessDirs)+1,1)
                TEXT = tList(AssmblPtrnCellIDs{k});
                text(0,0,TEXT,'Interpreter','None');
                axis tight
                axis off
            end
            
            % get each spike time
            tSp = event_Time{k};
            [spkx,spky,newTsp,spkPosInd] = spk2posInd(tSp,pdata.x,pdata.y,pdata.t);
            
            NaNidx = isnan(spkx);
            spkx(NaNidx)=[];
            spky(NaNidx)=[];
            newTsp(NaNidx)=[];
            spkPosInd(NaNidx)=[];
            
            % plot path & spikes
            subplot(4,length(sessDirs)+1,j+1)
            plotpathspike(pdata, spkPosInd)
            axis off
            axis tight
            if j == 1
                TITLE1 = sprintf('%d%s%d%s%d%s%s',sessInfo(i).animal,'-Day-',sessInfo(i).day,'-Assembly-',k,'-',sessDirs{j});
                title(TITLE1,'Interpreter','None')
            else
                TITLE1 = sprintf('%s',sessDirs{j});
                title(TITLE1,'Interpreter','None')
            end
            set(gca,'YDir','reverse')
            
            % plot rate map
            subplot(4,length(sessDirs)+1,j+length(sessDirs)+2)
            
            [map, pospdf] = ratemap(spkx,spky,pdata.x,pdata.y,pdata.t,p.smoothing,mapAxis1,mapAxis2);
            %%% clear non visited part
            % Set unvisited parts of the box to NaN
            map(visited==0) = NaN;
            drawfield(map,mapAxis1,mapAxis2,'jet',max(max(map)),p.binWidth,p.smoothing);
            set(gcf,'color',[1 1 1]);
            axis off;
            set(gca,'YDir','reverse')
            
            avgRate(k) = length(newTsp)/(pdata.t(end)-pdata.t(1));
            peakRate(k) = max(max(map));
            TITLE_1 = sprintf('%s%2.2f%s','AvgRate: ',avgRate(k),'Hz');
            TITLE_2 = sprintf('%s%2.2f%s','PeakRate: ',peakRate(k),'Hz');
            title({TITLE_1;TITLE_2});
            
            % plot idealized path & spikes
            subplot(4,length(sessDirs)+1,j+2*(length(sessDirs)+1)+1)
            plotpathspike(pIdeal, spkPosInd)
            axis off;
            set(gca,'YDir','reverse')
            
            % plot map generated by idealized path
            subplot(4,length(sessDirs)+1,j+3*(length(sessDirs)+1)+1)
            
            [spkxIdeal,spkyIdeal,newTspIdeal,spkPosIndIdeal] = ...
                spk2posInd(tSp,pIdeal.x,pIdeal.y,pIdeal.t);
            NaNidx = isnan(spkxIdeal);
            spkxIdeal(NaNidx)=[];
            spkyIdeal(NaNidx)=[];
            newTspIdeal(NaNidx)=[];
            spkPosIndIdeal(NaNidx)=[];
            
            [mapIdeal, pospdfIdeal] = ratemap(spkxIdeal,spkyIdeal,pIdeal.x,pIdeal.y,pIdeal.t,p.smoothing,mapAxisIdeal1,mapAxisIdeal2);
            %%% clear non visited part
            % Set unvisited parts of the box to NaN
            mapIdeal(visitedIdeal==0) = NaN;
            drawfield(mapIdeal,mapAxisIdeal1,mapAxisIdeal2,'jet',max(max(mapIdeal)),p.binWidth,p.smoothing);
            set(gcf,'color',[1 1 1]);
            axis off;
            set(gca,'YDir','reverse')
            
            % calculate lin version map and plot
            
            CellAssembly_Map.(sessDirs{j}).map{k} = map;
            CellAssembly_Map.(sessDirs{j}).mapPdf{k} = pospdf;
            CellAssembly_Map.(sessDirs{j}).mapIdeal{k} = mapIdeal;
            CellAssembly_Map.(sessDirs{j}).mapPdfIdeal{k} = pospdfIdeal;
        end
        CellAssembly_Map.(sessDirs{j}).avgRate = avgRate;
        CellAssembly_Map.(sessDirs{j}).peakRate = peakRate;
    end
    
    %% write files down
    if p.savePlot == 1
        for k = 1:patNum
            figure(k)
            set(gcf, 'PaperPositionMode', 'auto')
            figName = sprintf('%s%s%d%s%d%s%d%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-Assembly-',k,'-Spatial map');
            print(figName,'-dpng','-r300');
        end
    end
    close all
    
    if p.writeToFile
        save(fullfile(savedir2,'CellAssembly_Map.mat'), 'CellAssembly_Map');
    end
    clear CellAssembly_Map
    
    fprintf('Finished analysis for session %d\n',i)
   
end