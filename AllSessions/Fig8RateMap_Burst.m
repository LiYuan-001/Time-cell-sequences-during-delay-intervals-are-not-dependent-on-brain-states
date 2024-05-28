function Fig8RateMap_Burst(inFile,AnalyzeSes)
%   
%  Li Yuan, UCSD, 15-Nov-2019
%  plot general map for fig8 maze nd open box session
% -------------------------------------------------------------------------
% set parameters for analysis
close all

p.binWidth = 5;
p.smoothing = 5;

% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'

p.savePlot = 1;
p.writeToFile = 1;

% -------------------------------------------------------------------------

% Read in input information
sessInfo = SessInfoImport(inFile);

p.avgRatesDir = cd;

	
for i = AnalyzeSes(1:end)
   
    if p.savePlot
        % Li Yuan, UCSD, 15-Nov-2019
        
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\Spatial_map_Burst');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end

    display(['About to do session ' sessInfo(i).mainDir]);
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    % load burst and single spike
    burstFile = fullfile(sessInfo(i).mainDir,'Cell Property','BurstActivity.mat');
    load(burstFile); 
    TList = BurstActivity.TList;
    
    if prod(~isempty(char(sessDirs{1})))
        
        [~, pathData, pathDataIdeal, pathDataLin] = loadPathInfo(sessInfo(i));

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
            rateMap_Burst.mapAxis1 = mapAxis1;
            rateMap_Burst.mapAxis2 = mapAxis2;
            rateMap_Burst.mapAxisIdeal1 = mapAxisIdeal1;
            rateMap_Burst.mapAxisIdeal2 = mapAxisIdeal2;
            rateMap_Burst.visited = visited;
            rateMap_Burst.visitedIdeal = visitedIdeal;            
            rateMap_Burst.session = sessDirs{j};
            
            % start calculate rate map & plot  
            for k = 1:length(TList)

                h=figure(k); 
                if j ==1
                    h.Position = [100 100 1600 2000];
                end

                % get each spike time 
                % burst time
                % single spike time
                nBursts = BurstActivity.(sessDirs{j}).nBurst(k);
                if nBursts > 0
                    BurstTs = BurstActivity.(sessDirs{j}).BurstProp{k}.BurststartTs;
                    [spkx,spky,newTsp,spkPosInd] = spk2posInd(BurstTs,pdata.x,pdata.y,pdata.t);
                    
                    NaNidx = isnan(spkx);
                    spkx(NaNidx)=[];
                    spky(NaNidx)=[];
                    newTsp(NaNidx)=[];
                    spkPosInd(NaNidx)=[];
                    
                    % plot path & spikes
                    subplot(4,length(sessDirs),j)
                    plotpathspike(pdata, spkPosInd)
                    axis off
                    axis tight
                    if j == 1
                        TITLE1 = sprintf('%d%s%d%s%s%s%s',sessInfo(i).animal,'-Day-',sessInfo(i).day,'-',TList{k}(1:end-2),'-',sessDirs{j});
                        TITLE2 = 'Burst map';
                        title({TITLE1;TITLE2},'Interpreter','None')
                    else
                        title(strcat(TList{k}(1:end-2),'-',sessDirs{j}),'Interpreter','None');
                    end
                    set(gca,'YDir','reverse')
                    
                    % plot rate map
                    subplot(4,length(sessDirs),j+length(sessDirs))
                    
                    [map_burst, pospdf_burst] = ratemap(spkx,spky,pdata.x,pdata.y,pdata.t,p.smoothing,mapAxis1,mapAxis2);
                    %%% clear non visited part
                    % Set unvisited parts of the box to NaN
                    map_burst(visited==0) = NaN;
                    drawfield(map_burst,mapAxis1,mapAxis2,'jet',max(max(map_burst)),p.binWidth,p.smoothing);
                    set(gcf,'color',[1 1 1]);
                    axis off;
                    set(gca,'YDir','reverse')
                    
                    [spkxIdeal,spkyIdeal,newTspIdeal,spkPosIndIdeal] = ...
                        spk2posInd(BurstTs,pIdeal.x,pIdeal.y,pIdeal.t);
                    NaNidx = isnan(spkxIdeal);
                    spkxIdeal(NaNidx)=[];
                    spkyIdeal(NaNidx)=[];
                    newTspIdeal(NaNidx)=[];
                    spkPosIndIdeal(NaNidx)=[];
                    
                    [mapIdeal_burst, pospdfIdeal_burst] = ratemap(spkxIdeal,spkyIdeal,pIdeal.x,pIdeal.y,pIdeal.t,p.smoothing,mapAxisIdeal1,mapAxisIdeal2);
                else
                    map_burst = [];
                    pospdf_burst = [];
                    mapIdeal_burst = [];
                    pospdfIdeal_burst = [];
                end
                
                
                nSinglespikes = BurstActivity.(sessDirs{j}).nSinglespikes(k);
                if nSinglespikes > 0
                    
                    SinglespikeTs = BurstActivity.(sessDirs{j}).SinglespikeProp{k}.SinglespikeTs;
                    [spkx,spky,newTsp,spkPosInd] = spk2posInd(SinglespikeTs,pdata.x,pdata.y,pdata.t);
                    
                    NaNidx = isnan(spkx);
                    spkx(NaNidx)=[];
                    spky(NaNidx)=[];
                    newTsp(NaNidx)=[];
                    spkPosInd(NaNidx)=[];
                    
                    % plot path & spikes
                    subplot(4,length(sessDirs),j+2*(length(sessDirs)))
                    plotpathspike(pdata, spkPosInd)
                    axis off
                    axis tight
                    if j == 1
                        TITLE1 = sprintf('%d%s%d%s%s%s%s',sessInfo(i).animal,'-Day-',sessInfo(i).day,'-',TList{k}(1:end-2),'-',sessDirs{j});
                        TITLE2 = 'Single spike map';
                        title({TITLE1;TITLE2},'Interpreter','None')
                    else
                        title(strcat(TList{k}(1:end-2),'-',sessDirs{j}),'Interpreter','None');
                    end
                    set(gca,'YDir','reverse')
                    
                    % plot rate map
                    subplot(4,length(sessDirs),j+3*(length(sessDirs)))
                    
                    [map_singleSpk, pospdf_singleSpk] = ratemap(spkx,spky,pdata.x,pdata.y,pdata.t,p.smoothing,mapAxis1,mapAxis2);
                    %%% clear non visited part
                    % Set unvisited parts of the box to NaN
                    map_singleSpk(visited==0) = NaN;
                    drawfield(map_singleSpk,mapAxis1,mapAxis2,'jet',max(max(map_singleSpk)),p.binWidth,p.smoothing);
                    set(gcf,'color',[1 1 1]);
                    axis off;
                    set(gca,'YDir','reverse')
                    
                    [spkxIdeal,spkyIdeal,newTspIdeal,spkPosIndIdeal] = ...
                        spk2posInd(SinglespikeTs,pIdeal.x,pIdeal.y,pIdeal.t);
                    NaNidx = isnan(spkxIdeal);
                    spkxIdeal(NaNidx)=[];
                    spkyIdeal(NaNidx)=[];
                    newTspIdeal(NaNidx)=[];
                    spkPosIndIdeal(NaNidx)=[];
                    
                    [mapIdeal_singleSpk, pospdfIdeal_singleSpk] = ratemap(spkxIdeal,spkyIdeal,pIdeal.x,pIdeal.y,pIdeal.t,p.smoothing,mapAxisIdeal1,mapAxisIdeal2);
                else
                    map_singleSpk = [];
                    pospdf_singleSpk = [];
                    mapIdeal_singleSpk = [];
                    pospdfIdeal_singleSpk = [];
                end

                
                rateMap_Burst.map{k} = map_burst;
                rateMap_Burst.mapPdf{k} = pospdf_burst;
                rateMap_Burst.mapIdeal{k} = mapIdeal_burst;
                rateMap_Burst.mapPdfIdeal{k} = pospdfIdeal_burst;
                
                rateMap_Burst.map{k} = map_singleSpk;
                rateMap_Burst.mapPdf{k} = pospdf_singleSpk;
                rateMap_Burst.mapIdeal{k} = mapIdeal_singleSpk;
                rateMap_Burst.mapPdfIdeal{k} = pospdfIdeal_singleSpk;
                
                rateMap_Burst.tFileNames{k} = TList{k}(1:end-2);
            end  
            
            if p.writeToFile
                save(fullfile(sessInfo(i).mainDir,sessDirs{j}, 'rateMap_Burst.mat'), 'rateMap_Burst');
            end
            clear rateMap_Burst

        end
   
        if p.savePlot == 1
            for k = 1:length(TList)
                figure(k)
                figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',TList{k}(1:end-2),'-Fig8Map_Burst');
                print(figName,'-dpng','-r300');
            end
        end
        clear rateMap
        close all
    end
    
    % start ploting 
     OpenDirs = sessInfo(i).openDirs;
     
    if prod(~isempty(char(OpenDirs{1})))
        posFile = fullfile(sessInfo(i).mainDir,'processedData','indataO.mat');
        load(posFile);
        for j = 1:length(OpenDirs)
            pos = indata(j);
            
            % calculate original path rate map
            % Observed diameter/sidelength of the box
            obsLength = max(max(pos.x)-min(pos.x),max(pos.y)-min(pos.y));
            bins = ceil(obsLength/p.binWidth);
            sLength = p.binWidth*bins;
            mapAxis = (-sLength/2+p.binWidth/2):p.binWidth:(sLength/2-p.binWidth/2);

            % Calulate what areas of the box that have been visited
            visited = visitedBins2(pos.x,pos.y,mapAxis,mapAxis); 
            
            % for .mat files which is going to be write down
            rateMap_Burst.mapAxis = mapAxis;
            rateMap_Burst.visited = visited;      
            rateMap_Burst.session = OpenDirs{j};
            
            % start calculate rate map & plot  
            for k = 1:length(TList)

                h=figure(k); 
                if j ==1
                    h.Position = [100 100 length(OpenDirs)*600 length(OpenDirs)*300];
                end

                % get each spike time 
                tSp = Spike_Session.(OpenDirs{j}){k};
                [spkx,spky,newTsp,spkPosInd] = spk2posInd(tSp,pos.x,pos.y,pos.t);
                NaNidx = isnan(spkx);
                spkx(NaNidx)=[];
                spky(NaNidx)=[];
                newTsp(NaNidx)=[];
                spkPosInd(NaNidx)=[];
                
                % plot path & spikes
                subplot(length(OpenDirs),2,2*j-1)
                plotpathspike(pos, spkPosInd)
                axis off
                axis equal
                title(strcat(TList{k}(1:end-2),'-',OpenDirs{j}),'Interpreter','None')
                set (gca,'YDir','reverse')
                
                % plot rate map
                subplot(length(OpenDirs),2,2*j)
                
                [map_burst pospdf_burst] = ratemap(spkx,spky,pos.x,pos.y,pos.t,p.smoothing,mapAxis,mapAxis);
                %%% clear non visited part
                % Set unvisited parts of the box to NaN
                map_burst(visited==0) = NaN;               
                drawfield(map_burst,mapAxis,mapAxis,'jet',max(max(map_burst)),p.binWidth,p.smoothing);
                set(gcf,'color',[1 1 1]);
                axis off;
                axis equal
                TITLE1 = sprintf('%s%2.1f%s','Avg Rate: ',nanmean(nanmean(map_burst)),' Hz');
                TITLE2 = sprintf('%s%2.1f%s','Peak Rate: ',max(max(map_burst)),' Hz');
                title({TITLE1;TITLE2})
                set (gca,'YDir','reverse')        
                
                if p.savePlot == 1
                    figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',TList{k}(1:end-2),'-OpenBox');
                    print(figName,'-dpng','-r300');
                end
                
                rateMap_Burst.map{k} = map_burst;
                rateMap_Burst.mapPdf{k} = pospdf_burst;
                rateMap_Burst.tFileNames{k} = TList{k}(1:end-2);
            end
            
            if p.writeToFile
                save(fullfile(sessInfo(i).mainDir,OpenDirs{j}, 'GeneralMap.mat'), 'rateMap');
            end
            clear rateMap
        end
    end 
    close all
    fprintf('Finished ratemap analysis for session %d\n',i);
end 

end
