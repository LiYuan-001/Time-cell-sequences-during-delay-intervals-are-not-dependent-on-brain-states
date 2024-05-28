%  Li Yuan, UCSD, 3-May-2023
% -------------------------------------------------------------------------
function Fig8Treadmill_PCAICA_LFP(inFile,AnalyzeSes)
close all
% set parameters for analysis

p.savePlot = 1;
p.writeToFile = 1;
p.timeLen = 0.25; % unit sec
% parameters for import EEG from neuralynx / matlab utilities
FieldSelectionFlags = [1 0 0 0 1];
HeaderExtractionFlag = 1;
ExtractionMode = 1;
ExtractionModeVector = [];

p.simThre = 0.5;
onSession = {'on10_1','on10_2','on30_1','on30_2'};
offSession = {'off10_1','off10_2','off30_1','off30_2'};
session = {'on','off'};

% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    close all
    % selectivity for the whole map
    Assembly_LFPShape.session = sessInfo(i).mainDir;
    Assembly_LFPShape.rat = sessInfo(i).animal;
    Assembly_LFPShape.day = sessInfo(i).day;
    Assembly_LFPShape.timeBin = 25;
    
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\Assembly_LFP');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end
    
    if p.writeToFile
        savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
        if ~exist(savedir2, 'dir')
            mkdir(savedir2);
        end
    end
    
    display(['About to do session ' sessInfo(i).mainDir]);
    % get each phase names (no delay etc)
    sessDirs = sessInfo(i).sessDirs;
    
    % load assembly file
    assemblyFile = fullfile(sessInfo(i).mainDir,'Cell Property','CellAssembly_DelayOnOff_LR-25ms.mat');
    load(assemblyFile);
    
    binTime = CellAssembly_DelayLR.DelayOn.binTime;
    
    patNum_on = CellAssembly_DelayLR.DelayOn.patNum; 
    AssmblPtrnCellIDs_on = CellAssembly_DelayLR.DelayOn.AssmblPtrnCellIDs;
    AssmblWght_on = CellAssembly_DelayLR.DelayOn.AssmblWght;
    AssmblStrength_on = CellAssembly_DelayLR.DelayOn.AssmblStrength;
    event_Time_on = CellAssembly_DelayLR.DelayOn.event_Time;
    event_strength_on = CellAssembly_DelayLR.DelayOn.event_strength;
    event_Num_on = CellAssembly_DelayLR.DelayOn.event_Num;
    
    patNum_off = CellAssembly_DelayLR.DelayOff.patNum; 
    AssmblPtrnCellIDs_off = CellAssembly_DelayLR.DelayOff.AssmblPtrnCellIDs;
    AssmblWght_off = CellAssembly_DelayLR.DelayOff.AssmblWght;
    AssmblStrength_off = CellAssembly_DelayLR.DelayOff.AssmblStrength;
    event_Time_off = CellAssembly_DelayLR.DelayOff.event_Time;
    event_strength_off = CellAssembly_DelayLR.DelayOff.event_strength;
    event_Num_off = CellAssembly_DelayLR.DelayOff.event_Num;
    
    
    eegInd = sprintf('%s%d%s', 'CSC',sessInfo(i).EEGch2,'.ncs');
    cscFile = fullfile(sessInfo(i).mainDir,eegInd);
    % extract EEG by Nlx import
    [Timestamps,Samples,Header] = Nlx2MatCSC(cscFile, FieldSelectionFlags,HeaderExtractionFlag, ExtractionMode, ExtractionModeVector);
    p.t1 = Timestamps(1);
    
    % reshape EEG samples
    EEG_Raw=reshape(Samples,1,length(Samples(:)));
    % get Frequency and ADBvolts for this channel
    [p.Fs, p.ADBVolts,p.Invert] = ReadHeaderEEG2(Header);
    % get EEG to uV
    EEG_Raw = EEG_Raw*p.ADBVolts*10^6*p.Invert;
    % filter to 3-300 Hz;
    EEG = fftbandpass(EEG_Raw,p.Fs,4-1,4,300,300+1);
    % unit: sec
    EEGTs = timeStampResample(Samples,Timestamps,p.Fs)./10^6;   
    
    
    sim_on = zeros(1,patNum_on);
    sim_off = zeros(1,patNum_off);
    
    for m = 1:patNum_on
        sim_Temp = [];
        pattern_On = AssmblWght_on(:,m);
        for j = 1:patNum_off
            pattern_Off = AssmblWght_off(:,j);
            sim_Temp(j) = dot(pattern_On,pattern_Off)/(norm(pattern_On)*norm(pattern_Off));
        end
        sim_on(m) = max(sim_Temp);
    end
    
    for m = 1:patNum_off
        sim_Temp = [];
        pattern_Off = AssmblWght_off(:,m);
        for j = 1:patNum_on
            pattern_On = AssmblWght_on(:,j);
            sim_Temp(j) = dot(pattern_On,pattern_Off)/(norm(pattern_On)*norm(pattern_Off));
        end
        sim_off(m) = max(sim_Temp);
    end
    on_share = find(sim_on >= p.simThre);
    on_spec = find(sim_on < p.simThre);
    off_share = find(sim_off >= p.simThre);
    off_spec = find(sim_off < p.simThre);
    
    maxNum = max([length(on_share),length(on_spec),length(off_share),length(off_spec)]);
    
    h = figure(1);
    h.Position = [100,100,1200,900];
    subplot(maxNum+1,4,1)
    title('On share assembly')
    axis off
    subplot(maxNum+1,4,2)
    title('Off share assembly')
    axis off
    subplot(maxNum+1,4,3)
    title('On spec assembly')
    axis off
    subplot(maxNum+1,4,4)
    title('Off spec assembly')
    axis off
    
        if ~isempty(on_share)
            eventLFP_onShare = cell(1,length(on_share));
        end
        if ~isempty(off_share)
            eventLFP_offShare = cell(1,length(off_share));
        end
        if ~isempty(on_spec)
            eventLFP_onSpec = cell(1,length(on_spec));
        end
        if ~isempty(off_spec)
            eventLFP_offSpec = cell(1,length(off_spec));
        end
        
        for j = 1:length(onSession)
            sessName = onSession{j};
            delayFile = fullfile(sessInfo(i).mainDir,sessName,'Fig8DelayZonePos.mat');
            load(delayFile);
            
            % def1: delay starts at barrier
            % def2: delay starts at entrance
            delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
            delayTend1 = Fig8DelayZonePos.delayPos1.endT;
            %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
            %         delayTend2 = Fig8DelayZonePos.delayPos2.endT;
            
            trialNum = size(delayTstart1,2);
            if strcmp(sessName(end-3:end-2),'10')
                maxT = 10;
                delayTend1_2 = delayTstart1+maxT;
            elseif strcmp(sessName(end-3:end-2),'30')
                maxT = 30;
                delayTend1_2 = delayTstart1+maxT;
            else
                error('Delay time is wrong')
            end
            
            if ~isempty(on_share)
                for k = 1:length(on_share)
                    eventTime = event_Time_on{on_share(k)};
                    delay_event_Ind = find(any(eventTime>=delayTstart1(:) & eventTime<delayTend1_2(:)));
                    
                    for kk = 1:length(delay_event_Ind)
                        ts = eventTime(delay_event_Ind(kk))-p.timeLen;
                        [~,ts_startInd] = min(abs(ts-EEGTs));
                        ts = eventTime(delay_event_Ind(kk))+p.timeLen;
                        [~,ts_endInd] = min(abs(ts-EEGTs));
                        LFP_Temp = EEG(ts_startInd:ts_endInd);
                        eventLFP_onShare{k} = [eventLFP_onShare{k};LFP_Temp];                        
                    end
                end
            end
            
            if ~isempty(off_share)
                for k = 1:length(off_share)
                    eventTime = event_Time_off{off_share(k)};
                    delay_event_Ind = find(any(eventTime>=delayTstart1(:) & eventTime<delayTend1_2(:)));
                    
                    for kk = 1:length(delay_event_Ind)
                        ts = eventTime(delay_event_Ind(kk))-p.timeLen;
                        [~,ts_startInd] = min(abs(ts-EEGTs));
                        ts = eventTime(delay_event_Ind(kk))+p.timeLen;
                        [~,ts_endInd] = min(abs(ts-EEGTs));
                        LFP_Temp = EEG(ts_startInd:ts_endInd);
                        eventLFP_offShare{k} = [eventLFP_offShare{k};LFP_Temp];
                    end
                end
            end
            
             if ~isempty(on_spec)
                for k = 1:length(on_spec)
                    eventTime = event_Time_on{on_spec(k)};
                    delay_event_Ind = find(any(eventTime>=delayTstart1(:) & eventTime<delayTend1_2(:)));
                    
                    for kk = 1:length(delay_event_Ind)
                        ts = eventTime(delay_event_Ind(kk))-p.timeLen;
                        [~,ts_startInd] = min(abs(ts-EEGTs));
                        ts = eventTime(delay_event_Ind(kk))+p.timeLen;
                        [~,ts_endInd] = min(abs(ts-EEGTs));
                        LFP_Temp = EEG(ts_startInd:ts_endInd);
                        eventLFP_onSpec{k} = [eventLFP_onSpec{k};LFP_Temp];                        
                    end
                end
            end
            
            if ~isempty(off_spec)
                for k = 1:length(off_spec)
                    eventTime = event_Time_off{off_spec(k)};
                    delay_event_Ind = find(any(eventTime>=delayTstart1(:) & eventTime<delayTend1_2(:)));
                    
                    for kk = 1:length(delay_event_Ind)
                        ts = eventTime(delay_event_Ind(kk))-p.timeLen;
                        [~,ts_startInd] = min(abs(ts-EEGTs));
                        ts = eventTime(delay_event_Ind(kk))+p.timeLen;
                        [~,ts_endInd] = min(abs(ts-EEGTs));
                        LFP_Temp = EEG(ts_startInd:ts_endInd);
                        eventLFP_offSpec{k} = [eventLFP_offSpec{k};LFP_Temp];                        
                    end 
                end
                
            end   
        end
        
        % plot all on
        if ~isempty(on_share)
            for k = 1:length(on_share)
                figure(1)
                subplot(maxNum+1,8,k*8+1)
                stdshade(eventLFP_onShare{k},0.5,[1,0,0]);
            end
            Assembly_LFPShape.on.eventLFP_onShare = eventLFP_onShare;
        else
            Assembly_LFPShape.on.eventLFP_onShare = [];
        end
        
        if ~isempty(off_share)
            for k = 1:length(off_share)
                figure(1)
                subplot(maxNum+1,8,k*8+3)
                stdshade(eventLFP_offShare{k},0.5,[1,0,0]);
            end
            Assembly_LFPShape.on.eventLFP_offShare = eventLFP_offShare;
        else
            Assembly_LFPShape.on.eventLFP_offShare = [];
        end
        
        if ~isempty(on_spec)
            for k = 1:length(on_spec)
                figure(1)
                subplot(maxNum+1,8,k*8+5)
                stdshade(eventLFP_onSpec{k},0.5,[1,0,0]);
            end
            Assembly_LFPShape.on.eventLFP_onSpec = eventLFP_onSpec;
        else
            Assembly_LFPShape.on.eventLFP_onSpec = [];
        end
        
        if ~isempty(off_spec)
            for k = 1:length(off_spec)
                figure(1)
                subplot(maxNum+1,8,k*8+7)
                stdshade(eventLFP_offSpec{k},0.5,[1,0,0]);
            end
            Assembly_LFPShape.on.eventLFP_offSpec = eventLFP_offSpec;
        else
            Assembly_LFPShape.on.eventLFP_offSpec = [];
        end
            
        %%  
        if ~isempty(on_share)
            eventLFP_onShare = cell(1,length(on_share));
        end
        if ~isempty(off_share)
            eventLFP_offShare = cell(1,length(off_share));
        end
        if ~isempty(on_spec)
            eventLFP_onSpec = cell(1,length(on_spec));
        end
        if ~isempty(off_spec)
            eventLFP_offSpec = cell(1,length(off_spec));
        end
        
        for j = 1:length(offSession)
            sessName = offSession{j};
            delayFile = fullfile(sessInfo(i).mainDir,sessName,'Fig8DelayZonePos.mat');
            load(delayFile);
            
            % def1: delay starts at barrier
            % def2: delay starts at entrance
            delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
            delayTend1 = Fig8DelayZonePos.delayPos1.endT;
            %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
            %         delayTend2 = Fig8DelayZonePos.delayPos2.endT;
            
            trialNum = size(delayTstart1,2);
            if strcmp(sessName(end-3:end-2),'10')
                maxT = 10;
                delayTend1_2 = delayTstart1+maxT;
            elseif strcmp(sessName(end-3:end-2),'30')
                maxT = 30;
                delayTend1_2 = delayTstart1+maxT;
            else
                error('Delay time is wrong')
            end
            
            if ~isempty(on_share)
                for k = 1:length(on_share)
                    eventTime = event_Time_on{on_share(k)};
                    delay_event_Ind = find(any(eventTime>=delayTstart1(:) & eventTime<delayTend1_2(:)));
                    
                    for kk = 1:length(delay_event_Ind)
                        ts = eventTime(delay_event_Ind(kk))-p.timeLen;
                        [~,ts_startInd] = min(abs(ts-EEGTs));
                        ts = eventTime(delay_event_Ind(kk))+p.timeLen;
                        [~,ts_endInd] = min(abs(ts-EEGTs));
                        LFP_Temp = EEG(ts_startInd:ts_endInd);
                        eventLFP_onShare{k} = [eventLFP_onShare{k};LFP_Temp];                        
                    end
                end
            end
            
            if ~isempty(off_share)
                for k = 1:length(off_share)
                    eventTime = event_Time_off{off_share(k)};
                    delay_event_Ind = find(any(eventTime>=delayTstart1(:) & eventTime<delayTend1_2(:)));
                    
                    for kk = 1:length(delay_event_Ind)
                        ts = eventTime(delay_event_Ind(kk))-p.timeLen;
                        [~,ts_startInd] = min(abs(ts-EEGTs));
                        ts = eventTime(delay_event_Ind(kk))+p.timeLen;
                        [~,ts_endInd] = min(abs(ts-EEGTs));
                        LFP_Temp = EEG(ts_startInd:ts_endInd);
                        eventLFP_offShare{k} = [eventLFP_offShare{k};LFP_Temp];                        
                    end
                end
            end
            
             if ~isempty(on_spec)
                for k = 1:length(on_spec)
                    eventTime = event_Time_on{on_spec(k)};
                    delay_event_Ind = find(any(eventTime>=delayTstart1(:) & eventTime<delayTend1_2(:)));
                    
                    for kk = 1:length(delay_event_Ind)
                        ts = eventTime(delay_event_Ind(kk))-p.timeLen;
                        [~,ts_startInd] = min(abs(ts-EEGTs));
                        ts = eventTime(delay_event_Ind(kk))+p.timeLen;
                        [~,ts_endInd] = min(abs(ts-EEGTs));
                        LFP_Temp = EEG(ts_startInd:ts_endInd);
                        eventLFP_onSpec{k} = [eventLFP_onSpec{k};LFP_Temp];                        
                    end
                end
             end
            
            if ~isempty(off_spec)
                for k = 1:length(off_spec)
                    eventTime = event_Time_off{off_spec(k)};
                    delay_event_Ind = find(any(eventTime>=delayTstart1(:) & eventTime<delayTend1_2(:)));
                    
                    for kk = 1:length(delay_event_Ind)
                        ts = eventTime(delay_event_Ind(kk))-p.timeLen;
                        [~,ts_startInd] = min(abs(ts-EEGTs));
                        ts = eventTime(delay_event_Ind(kk))+p.timeLen;
                        [~,ts_endInd] = min(abs(ts-EEGTs));
                        LFP_Temp = EEG(ts_startInd:ts_endInd);
                        eventLFP_offSpec{k} = [eventLFP_offSpec{k};LFP_Temp];                        
                    end
                end
            end
        end

        % plot all off
        if ~isempty(on_share)
            for k = 1:length(on_share)
                figure(1)
                subplot(maxNum+1,8,k*8+2)
                stdshade(eventLFP_onShare{k},0.5,[0,0,0]);
            end
            Assembly_LFPShape.off.eventLFP_onShare = eventLFP_onShare;
        else
            Assembly_LFPShape.off.eventLFP_onShare = [];
        end
        
        if ~isempty(off_share)
            for k = 1:length(off_share)
                figure(1)
                subplot(maxNum+1,8,k*8+4)
                stdshade(eventLFP_offShare{k},0.5,[0,0,0]);
            end
            Assembly_LFPShape.off.eventLFP_offShare = eventLFP_offShare;
        else
            Assembly_LFPShape.off.eventLFP_offShare = [];
        end
        
        if ~isempty(on_spec)
            for k = 1:length(on_spec)
                figure(1)
                subplot(maxNum+1,8,k*8+6)
                stdshade(eventLFP_onSpec{k},0.5,[0,0,0]);
            end
            Assembly_LFPShape.off.eventLFP_onSpec = eventLFP_onSpec;
        else
            Assembly_LFPShape.off.eventLFP_onSpec = [];
        end
        
        if ~isempty(off_spec)
            for k = 1:length(off_spec)
                figure(1)
                subplot(maxNum+1,8,k*8+8)
                stdshade(eventLFP_offSpec{k},0.5,[0,0,0]);
            end
            Assembly_LFPShape.off.eventLFP_offSpec = eventLFP_offSpec;
        else
            Assembly_LFPShape.off.eventLFP_offSpec = [];
        end
            
            
    % save figure and save .mat file
    if p.savePlot
        figure(1)
        set(gcf, 'PaperPositionMode', 'auto')
        figName = sprintf('%s%s%d%s%d%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-AssemblyLFP');
        print(figName,'-dpng','-r300');
    end
    
    if p.writeToFile
        save(fullfile(savedir2,'Assembly_LFPShape.mat'), 'Assembly_LFPShape');
    end
    clear Assembly_LFPShape
    close all
    fprintf('Finished assembly LFP analysis %d\n',i);
end

end
