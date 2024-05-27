% Theta is defined by amp + length definition
% Li Yuan, Feb-22-2022, UCSD
function Fig8TreadmillEEGThetaLength(inFile,AnalyzeSes)

close all

p.savePlot = 0;
p.writeToFile = 1;

% parameters for import EEG from neuralynx / matlab utilities
% set parameter
p.delta = [1 4];
p.theta = [6 12];

p.thetaLength = 0.5; % unit: sec
p.thetaThres = [0.3,0.5,0.7,1]; % mean + x*std;


FieldSelectionFlags = [1 0 0 0 1];
HeaderExtractionFlag = 1;
ExtractionMode = 1;
ExtractionModeVector = [];


% Read in input information
sessInfo = SessInfoImport(inFile);  

for i = AnalyzeSes(1:end)
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\Theta Property');
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
    
    % initiate the data 
    ThetaProperty_Length.rat = sessInfo(i).animal;
    ThetaProperty_Length.day = sessInfo(i).day;
    ThetaProperty_Length.thetaLengthLimit = p.thetaLength;
    ThetaProperty_Length.thetaThres = p.thetaThres;
    
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    sessDirs = sessInfo(i).sessDirs;
     
    % load eeg file
    eegInd = sprintf('%s%d%s', 'CSC',sessInfo(i).EEGch,'.ncs');
    cscFile = fullfile(sessInfo(i).mainDir,eegInd);
    % extract EEG by Nlx import
    [Timestamps,Samples,Header] = Nlx2MatCSC(cscFile, FieldSelectionFlags,HeaderExtractionFlag, ExtractionMode, ExtractionModeVector);
    p.t1 = Timestamps(1);
    
    
    % reshape EEG samples
    EEG=reshape(Samples,length(Samples(:)),1);
    % get Frequency and ADBvolts for this channel
    [p.Fs, p.ADBVolts] = ReadHeaderEEG(Header);
    % get EEG to uV
    EEG = EEG*p.ADBVolts*10^6;
    EEGt = fftbandpass(EEG,p.Fs,p.theta(1)-1,p.theta(1),p.theta(2),p.theta(2)+1);
    [phaseT,AmpT] = thetaPhase2(EEGt);  
    
%     % theta phase accumulate from each cycle
%     % detect theta cycles
%     ThetaCycle = cycleLabel(phaseT);
%     % cycle in each region and spike count
%     cycleCount = length(unique(ThetaCycle.Cycle));
%     [~,cycleStartIdx] = unique(ThetaCycle.Cycle);
%     cycleEndIdx = cycleStartIdx-1;
%     cycleEndIdx = cycleEndIdx(2:end);
%     cycleEndIdx(cycleCount) = length(ThetaCycle.Cycle);
        
    % unit: sec
    EEGTs = timeStampResample(Samples,Timestamps,p.Fs)./10^6;          
    
    
    for n = 1:length(p.thetaThres)
        % calculate theta oscillation based on the threshold set
        [ThetaStartInd,ThetaEndInd,tsStart,tsEnd] = ThetaFinder3(EEGt,EEGTs,p.Fs,p.thetaThres(n),p.thetaLength);
        thresName = sprintf('%s%d','Thres',n);
        
        for j = 1:length(sessDirs)
            % load analyzed positions
            delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
            load(delayFile);
            % load path zone time
            pathZoneFile = fullfile(mainDir,sessDirs{j}, 'PathZone.mat');
            load(pathZoneFile);
            % load session x,y,t
            pathDataFile = fullfile(mainDir,sessDirs{j},'pathData.mat');
            pathData = load(pathDataFile);
            % load turning directions and correctness
            tInfoFile = fullfile(mainDir,sessDirs{j}, 'trialInfo.mat');
            tInfo = load(tInfoFile);
            
            % return area
            regTimeBorder_return = [PathZone.posStartT.Return,PathZone.posEndT.Return]; % return+base
            % stem area
            regTimeBorder_stem = [Fig8DelayZonePos.delayPos1.endT'+0.02,PathZone.posEndT.Center]; % return+base
            
            % def1: delay starts at barrier
            % def2: delay starts at entrance
            delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
            delayTend1 = Fig8DelayZonePos.delayPos1.endT;
            %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
            %         delayTend2 = Fig8DelayZonePos.delayPos2.endT;
            
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
            
            
            % calculate for delay, stem and return arm
            thetaStart_return = cell(trialNum,1);
            thetaEnd_return = cell(trialNum,1);
            thetaStart_stem = cell(trialNum,1);
            thetaEnd_stem = cell(trialNum,1);
            thetaStart_delay = cell(trialNum,1);
            thetaEnd_delay = cell(trialNum,1);
            
            % each trial
            for m = 1:trialNum
                
                %             % delay area
                %             [startTimeDifference,delaystartIdx] = min(abs(timeStamps-Fig8DelayZonePos.delayPos2.startT(k)));
                %             [endTimeDifference,delayendIdx] = min(abs(timeStamps-Fig8DelayZonePos.delayPos2.endT(k)));
                %             % accumulate delay Idx from each trial
                %             delayIdx = [delayIdx,delaystartIdx:delayendIdx];
                %
                %             if any([startTimeDifference,endTimeDifference] > (timeStamps(2)-timeStamps(1)))
                %                 error('EEG could not find matching position')
                %             end
                
                
                % return
                thetaStartTemp1 = [];
                thetaStartTemp2 = [];
                thetaEndTemp1 = [];
                thetaEndTemp2 = [];
                % return + base spikes under theta
                % condition 1: theta start after location start
                thetaInd = tsStart >= regTimeBorder_return(m,1) & tsStart <= regTimeBorder_return(m,2);
                if sum(thetaInd) > 0
                    thetaStartTemp1 = tsStart(thetaInd);
                    thetaEndTemp1 = tsEnd(thetaInd);
                    if thetaEndTemp1(end) > regTimeBorder_return(m,2)
                        thetaEndTemp1(end) = regTimeBorder_return(m,2);
                    end
                end
                % condition2: theta start before location start
                if any(tsStart < regTimeBorder_return(m,1) & tsEnd > regTimeBorder_return(m,1))
                    thetaInd = (tsStart < regTimeBorder_return(m,1) & tsEnd > regTimeBorder_return(m,1));
                    thetaStartTemp2 = regTimeBorder_return(m,1);
                    thetaEndTemp2 = tsEnd(thetaInd);
                    if thetaEndTemp2(end) > regTimeBorder_return(m,2)
                        thetaEndTemp2(end) = regTimeBorder_return(m,2);
                    end
                end
                
                thetaStart_return{m} = sort([thetaStartTemp1;thetaStartTemp2]);
                thetaEnd_return{m} = sort([thetaEndTemp1;thetaEndTemp2]);
                
                %--------------------------------------------------------------
                % delay area spike under theta
                thetaStartTemp1 = [];
                thetaStartTemp2 = [];
                thetaEndTemp1 = [];
                thetaEndTemp2 = [];
                % condition 1: theta start after location start
                thetaInd = tsStart >= delayTstart1(m) & tsStart <= delayTend1_2(m);
                if sum(thetaInd) > 0
                    thetaStartTemp1 = tsStart(thetaInd);
                    thetaEndTemp1 = tsEnd(thetaInd);
                    if thetaEndTemp1(end) > delayTend1_2(m)
                        thetaEndTemp1(end) = delayTend1_2(m);
                    end
                end
                % condition2: theta start before location start
                if any(tsStart < delayTstart1(m) & tsEnd > delayTstart1(m))
                    thetaInd = (tsStart < delayTstart1(m) & tsEnd > delayTstart1(m));
                    thetaStartTemp2 = delayTstart1(m);
                    thetaEndTemp2 = tsEnd(thetaInd);
                    if thetaEndTemp2(end) > delayTend1_2(m)
                        thetaEndTemp2(end) = delayTend1_2(m);
                    end
                end
                
                thetaStart_delay{m} = sort([thetaStartTemp1;thetaStartTemp2]);
                thetaEnd_delay{m} = sort([thetaEndTemp1;thetaEndTemp2]);
            end
            
            ThetaProperty_Length.(sessDirs{j}).thetaStart_return.(thresName) = thetaStart_return;
            ThetaProperty_Length.(sessDirs{j}).thetaEnd_return.(thresName) = thetaEnd_return;
            ThetaProperty_Length.(sessDirs{j}).thetaStart_delay.(thresName) = thetaStart_delay;
            ThetaProperty_Length.(sessDirs{j}).thetaEnd_delay.(thresName) = thetaEnd_delay;
            ThetaProperty_Length.(sessDirs{j}).ts_return.(thresName) = regTimeBorder_return;
            ThetaProperty_Length.(sessDirs{j}).ts_delay.(thresName) = [delayTstart1',delayTend1_2'];
        end
    end

    if p.writeToFile
        save(fullfile(savedir2,'ThetaProperty_Length.mat'), 'ThetaProperty_Length');
    end
    fprintf('Finished EEG power analysis for session %d\n',i);
    clear ThetaProperty_Length
end
end