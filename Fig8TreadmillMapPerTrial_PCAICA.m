function Fig8TreadmillMapPerTrial_PCAICA(inFile,AnalyzeSes)
% 
%  Li Yuan, UCSD, 09-May-2022
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
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\Assembly spatial map per trial');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end

    display(['About to do session ' sessInfo(i).mainDir]);
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    % load average firing rate file
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    % load assembly file
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_Delay.mat');
    load(assemblyFile);
    
    clusterNum = length(SpikeProp.max_AvgRate);   
    rateCluster = SpikeProp.max_AvgRate;
    rateLabel = rateCluster<5 & rateCluster>0.1;    
    tList = SpikeProp.tList(rateLabel);
    
    patNum = CellAssembly_Delay.patNum ;
    AssmblPtrnCellIDs = CellAssembly_Delay.AssmblPtrnCellIDs;
    AssmblWght = CellAssembly_Delay.AssmblWght;
    AssmblStrength = CellAssembly_Delay.AssmblStrength;
    event_Time = CellAssembly_Delay.event_Time;
    event_strength = CellAssembly_Delay.event_strength;
    event_Num = CellAssembly_Delay.event_Num;
    
    if prod(~isempty(char(sessDirs{1})))
        
        [trialInfo, pathData, pathDataIdeal, pathDataLin] = loadPathInfo_Treadmill(sessInfo(i));

        for j = 1:length(sessDirs)
            
            tinfo = trialInfo.(sessDirs{j});
%             pinfo = parsingInfo.(SessDirs{j});
            pdata = pathData.(sessDirs{j});
            pIdeal = pathDataIdeal.(sessDirs{j});
            plin = pathDataLin.(sessDirs{j});

                % get each spike time 
                tSp = event_Time;
                newTsp = cell(patNum,1);
                for k = 1:patNum
                    [spkx,spky,newTsp{k},spkPosInd] = spk2posInd(tSp{k},pdata.x,pdata.y,pdata.t);
                    NaNidx = isnan(spkx);
                    newTsp{k} = newTsp{k}(~NaNidx);
                end
                
                tridx = trial_startend_ind(tinfo);
                tridx = tridx(~tinfo.degen, :);
                trtype = plotTrialPerCell.categtable(tinfo);
                
                [M, xrange, ~, preRM] =plotTrialPerCell.computeLinMaps(plin, newTsp, tridx, p);
                
                % for Stefan (requested on April 15, 2015)
                mux = @(x,y) x*2+~y+1; % simple multiplexer function of direction and success
                assemblyRatesByECLR.maps = M;
                assemblyRatesByECLR.xRange = xrange;
                assemblyRatesByECLR.ECLR = mux(strcmpi('R', tinfo.direction), tinfo.success);
                
                assemblyRatesByECLR.spikemap = arrayfun(@(str) str.spikemap, preRM, 'un', 0);
                assemblyRatesByECLR.occupmap = arrayfun(@(str) str.occupmap, preRM, 'un', 0);
                assemblyRatesByECLR.valid = ~tinfo.degen;
                
                % calculate average linear map of left correct turn and
                % right correct turn
                
                if p.writeToFile
                    save(fullfile(sessInfo(i).mainDir,sessDirs{j}, 'assemblyRatesByECLR.mat'), 'assemblyRatesByECLR');
                end
             
%                 [avgRates, zonespikes] = plotTrialPerCell.computeAvgRates2(avgRates, tSp, tinfo, sessInfo, i, TList, SessDirs, group, eegx, opt, pdata, tridx, trtype);
                close all
                doPlot2_Assembly; % plot and save figures, Li Yuan, 04-Feb-2020  
                clear assemblyRatesByECLR
        end

    end  
    fprintf('Finished position analysis for session %d\n',i);
    close all
end 

end