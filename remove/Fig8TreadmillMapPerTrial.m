function Fig8TreadmillMapPerTrial(inFile,AnalyzeSes)
% 
%  Li Yuan, UCSD, 23-Jan-2020
%  plot trial per session map for fig8 maze
% -------------------------------------------------------------------------
%  set parameters for analysis

p.writeToFile = 1; % write ELRC file

p.binWidth = 5;
p.smoothing = 5;

% whether read in spike time from whole session or subsessions
p.spikeMode = 'subses'; % the other value is 'wholeses'

p.spatialRange.left = 0;
p.spatialRange.right = 360;
p.spatialRange.bottom = -2.5;
p.spatialRange.top = 2.5;

p.nBins.x = 72;
p.nBins.y = 1;

p.boundLabels = {'Return','Dl','St','Ch','Rw'};

p.plotEachRun = 1; % Whether plot each run or not

p.savePlot = 1;

p.PATH_RECT_POSITION = [-65 -85 130 185];

[~, linBnds] = fig8trialtemplate2();
p.boundLabels = {'Return','Dl','St','Ch','Rw'};
ptempl = parsingtemplate_Treadmill('fig8:rdscr','fig8rat');
p.ROI = 'delay'; 
roixy = ptempl(strcmp({ptempl.zone}, p.ROI));
linROI = plotTrialPerCell.roi2lin(roixy, fig8trialtemplate2);
linChoice = plotTrialPerCell.roi2lin(ptempl(strcmp({ptempl.zone}, 'choice')), fig8trialtemplate2);
xTickLin = midpt(sort([0 360 linROI linChoice]));
xTickLabelLin = p.boundLabels;

avgRates = [];

MidPanelW = .25;
MidPanelH = .15;
p.figpanels.mid.pos = [[.1 .45 MidPanelW MidPanelH]; [.1 .25 MidPanelW MidPanelH]; [.38 .45 MidPanelW MidPanelH]; [.38 .25 MidPanelW MidPanelH]];
p.figpanels.mid.pos(:, 2) = p.figpanels.mid.pos(:, 2) + .04;

% -------------------------------------------------------------------------

% Read in input information
sessInfo = SessInfoImport(inFile);
	
for i = AnalyzeSes(1:end)
   
    if p.savePlot
        % Li Yuan, UCSD, 15-Nov-2019
        
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\Spatial map per trial');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end

    display(['About to do session ' sessInfo(i).mainDir]);
    % get each phase names (no delay etc)
    SessDirs = sessInfo(i).sessDirs;
    
    % load spikes from each main session
    % get event timestamp
    % assign spikes into different events  
    [Spike_Session, TList] = loadSpikeTime_Treadmill(sessInfo(i),p);
    clusterNum = length(TList);
    if length(unique(TList))~=clusterNum
        error('TTList file has repeated clusters')
    end
    if prod(~isempty(char(SessDirs{1})))
        
        [trialInfo, pathData, pathDataIdeal, pathDataLin] = loadPathInfo_Treadmill(sessInfo(i));

        for j = 1:length(SessDirs)
            
            tinfo = trialInfo.(SessDirs{j});
%             pinfo = parsingInfo.(SessDirs{j});
            pdata = pathData.(SessDirs{j});
            pIdeal = pathDataIdeal.(SessDirs{j});
            plin = pathDataLin.(SessDirs{j});

                % get each spike time 
                tSp = Spike_Session.(SessDirs{j});
%                 [spkx,spky,newTsp,spkPosInd] = spk2posInd(tSp,pdata.x,pdata.y,pdata.t);
%                 
%                 NaNidx = isnan(spkx);
%                 spkx(NaNidx)=[];
%                 spky(NaNidx)=[];
%                 newTsp(NaNidx)=[];
%                 spkPosInd(NaNidx)=[];
                
                tridx = trial_startend_ind(tinfo);
                tridx = tridx(~tinfo.degen, :);
                trtype = plotTrialPerCell.categtable(tinfo);
                
                [M, xrange, ~, preRM] =plotTrialPerCell.computeLinMaps(plin, tSp, tridx, p);
                
                % for Stefan (requested on April 15, 2015)
                mux = @(x,y) x*2+~y+1; % simple multiplexer function of direction and success
                ratesByECLR.maps = M;
                ratesByECLR.tFileNames = TList;
                ratesByECLR.xRange = xrange;
                ratesByECLR.ECLR = mux(strcmpi('R', tinfo.direction), tinfo.success);
                
                ratesByECLR.spikemap = arrayfun(@(str) str.spikemap, preRM, 'un', 0);
                ratesByECLR.occupmap = arrayfun(@(str) str.occupmap, preRM, 'un', 0);
                ratesByECLR.valid = ~tinfo.degen;
                
                % calculate average linear map of left correct turn and
                % right correct turn
                
                if p.writeToFile
                    save(fullfile(sessInfo(i).mainDir,SessDirs{j}, 'ratesByECLR.mat'), 'ratesByECLR');
                end
             
%                 [avgRates, zonespikes] = plotTrialPerCell.computeAvgRates2(avgRates, tSp, tinfo, sessInfo, i, TList, SessDirs, group, eegx, opt, pdata, tridx, trtype);
                close all
                doPlot2; % plot and save figures, Li Yuan, 04-Feb-2020      
        end
   
%         for k = 1:length(TList)
%             figure(k)            
%             figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',TList{k}(1:end-2),'-Fig8MapPerTrial');
%             print(figName,'-dpng','-r300');
%         end
        
%         close all
    end  
    fprintf('Finished position analysis for session %d\n',i);
    close all
end 

end
