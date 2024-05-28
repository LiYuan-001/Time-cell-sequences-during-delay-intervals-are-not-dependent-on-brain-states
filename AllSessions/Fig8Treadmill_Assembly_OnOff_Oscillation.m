% this code extract EEG period based on the evelope peak + threshold set
% Li Yuan, 27-Oct-2020
function Fig8Treadmill_Assembly_OnOff_Oscillation(inFile,AnalyzeSes)

% set parameters for analysis
p.writeToFile = 1;

% % Theta bandpass
% p.Fs1t = 5; p.Fp1t = 6; p.Fp2t = 11; p.Fs2t = 12;
% % Slow gamma bandpass
% p.Fs1s = 28; p.Fp1s = 30; p.Fp2s = 50; p.Fs2s = 52;
% % Mid gamma bandpass
% p.Fs1m = 53; p.Fp1m = 55; p.Fp2m = 90; p.Fs2m = 92;
% % Fast gamma bandpass
% p.Fs1f = 93; p.Fp1f = 95; p.Fp2f = 140; p.Fs2f = 142;
% Ripples bandpass
p.Fs1r = 148; p.Fp1r = 150; p.Fp2r = 250; p.Fs2r = 252;

% parameters for import EEG from neuralynx / matlab utilities
FieldSelectionFlags = [1 0 0 0 1];
HeaderExtractionFlag = 1;
ExtractionMode = 1;
ExtractionModeVector = [];

% Read in input information
sessInfo = SessInfoImport(inFile);  

for i = AnalyzeSes(1:end)
%     if p.savePlot
%         % directory for plot figures
%         % generate a folder for each rat eah day under the current folder
%         savedir = sprintf('%s%s%d%s%d%s',cd,'\',sessInfo(i).animal,'-day',sessInfo(i).day,'\EEGSpec');
%         if ~exist(savedir, 'dir')
%             mkdir(savedir);
%         end
%     end
    
    Fig8_Assembly_OnOff_SWR_Rate.session = sessInfo(i).mainDir;
    Fig8_Assembly_OnOff_SWR_Rate.rat = sessInfo(i).animal;
    Fig8_Assembly_OnOff_SWR_Rate.day = sessInfo(i).day;
    
    if p.writeToFile
        savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
        if ~exist(savedir2, 'dir')
            mkdir(savedir2);
        end
    end
    
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    sessDirs = sessInfo(i).sessDirs;
    sleepDirs = sessInfo(i).sleepDirs;
    
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
    
    % load eeg file
    eegInd = sprintf('%s%d%s', 'CSC',sessInfo(i).EEGch2,'.ncs');
    cscFile = fullfile(mainDir,eegInd);
    % extract EEG by Nlx import
    [Timestamps,Samples,Header] = Nlx2MatCSC(cscFile, FieldSelectionFlags,HeaderExtractionFlag, ExtractionMode, ExtractionModeVector);
    p.t1 = Timestamps(1);
    
    % reshape EEG samples
    EEG=reshape(Samples,length(Samples(:)),1);
    % get Frequency and ADBvolts for this channel
    [p.Fs, p.ADBVolts] = ReadHeaderEEG(Header);
    % reasample timestamps;
    % unit: sec
    EEGTs = timeStampResample(Samples,Timestamps,p.Fs)./10^6;   
    % get EEG to uV
    EEG = EEG*p.ADBVolts*10^6;
        
    % Band pass filter
%     eeg.T = fftbandpass(EEG,p.Fs, p.Fs1t, p.Fp1t, p.Fp2t, p.Fs2t); % Theta
%     eeg.S = fftbandpass(EEG,p.Fs, p.Fs1s, p.Fp1s, p.Fp2s, p.Fs2s); % Slow-gamma
%     eeg.M = fftbandpass(EEG,p.Fs, p.Fs1m, p.Fp1m, p.Fp2m, p.Fs2m); % Mid-gamma
%     eeg.F = fftbandpass(EEG,p.Fs, p.Fs1f, p.Fp1f, p.Fp2f, p.Fs2f); % Fast-gamma
    eeg.R = fftbandpass(EEG,p.Fs, p.Fs1r, p.Fp1r, p.Fp2r, p.Fs2r); % ripples
    
    % Calculate eeg phase based on TROUGHs with Hilbert transform
%     disp('Phase & Amplitute detection with hilbert transform')
%     [phase.T,amp.T] = assignPhase2(eeg.T);
%     [phase.S,amp.S] = assignPhase2(eeg.S);
%     [phase.M,amp.M] = assignPhase2(eeg.M);
%     [phase.F,amp.F] = assignPhase2(eeg.F);
%     [phase.R,amp.R] = assignPhase2(eeg.R);
    
    % Calculate gamma period
    % function [evtOnIdx,evtOffIdx,eegLabel] = eventFinder(LFP_Raw,LFP_Fil,Fs,ThresSD,TJoingtLim,TLthLim,PLOT)
%     [evtOnIdx.T,evtOffIdx.T,eegLabel.T] = eventFinder(EEG,eeg.T,p.Fs,0.3,0.5,0.5,0);
%     [evtOnIdx.S,evtOffIdx.S,eegLabel.S] = eventFinder(EEG,eeg.S,p.Fs,1,0.1,0.1,0);
%     [evtOnIdx.M,evtOffIdx.M,eegLabel.M] = eventFinder(EEG,eeg.M,p.Fs,1,0.1,0.1,0);
%     [evtOnIdx.F,evtOffIdx.F,eegLabel.F] = eventFinder(EEG,eeg.F,p.Fs,1,0.05,0.05,0);
%     [evtOnIdx.R,evtOffIdx.R,eegLabel.R] = eventFinder_SWR(EEG,eeg.R,p.Fs,3,0.05,0.020,0);
%     swr_Ts = EEGTs(eegLabel.R==1);
%     swr_Ts_Start = EEGTs(evtOnIdx.R);
%     
    % Hilbert transform
    Z = hilbert(eeg.R);    
    % Wave amplitude
    amp_SWR = abs(Z);
    % POWER
    pow_SWR = amp_SWR.^2;
    
    event_swr_pow_on = cell(patNum_on,1);
    % calculate swr strength during the assembly event
    for k = 1:patNum_on
        event_Ts = event_Time_on{k};
        for mm = 1:length(event_Ts)
            start_ts_temp = event_Ts(mm)-0.1;
            end_ts_temp = event_Ts(mm)+0.2;
            pow_Temp = nanmean(pow_SWR(EEGTs>=start_ts_temp & EEGTs<=end_ts_temp));
            event_swr_pow_on{k} = [event_swr_pow_on{k},pow_Temp];
        end
    end
    
    event_swr_pow_off = cell(patNum_off,1);
    % calculate swr strength during the assembly event
    for k = 1:patNum_off
        event_Ts = event_Time_off{k};
        for mm = 1:length(event_Ts)
            start_ts_temp = event_Ts(mm)-0.1;
            end_ts_temp = event_Ts(mm)+0.2;
            pow_Temp = nanmean(pow_SWR(EEGTs>=start_ts_temp & EEGTs<=end_ts_temp));
            event_swr_pow_off{k} = [event_swr_pow_off{k},pow_Temp];
        end
    end
    
    
    for j = 1:length(sessDirs)
     
        % load maps for each cluster
        pathZoneFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'PathZone.mat');
        load(pathZoneFile);
        delayFile = fullfile(sessInfo(i).mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        
        % def1: delay starts at barrier
        % def2: delay starts at entrance
        delayTstart1 = Fig8DelayZonePos.delayPos1.startT; % change unit to ms
        delayTend1 = Fig8DelayZonePos.delayPos1.endT;
        %         delayTstart2 = Fig8DelayZonePos.delayPos2.startT;
        %         delayTend2 = Fig8DelayZonePos.delayPos2.endT;
        
        trialNum = size(delayTstart1,2);
        if strcmp(sessDirs{j}(end-3:end-2),'10')
            maxT = 10;
            delayTend1_2 = delayTstart1+maxT;
        elseif strcmp(sessDirs{j}(end-3:end-2),'30')
            maxT = 30;
            delayTend1_2 = delayTstart1+maxT;
        else
            error('Delay time is wrong')
        end
        
        %% assembly rate and matching swr power
%         pat_Delay_Event = zeros(patNum_on,length(bin_Time_SWR));
%         pat_Delay_Strength = zeros(patNum_on,length(bin_Time_SWR));
        for k = 1:patNum_on
            % get event time and strength from the assembly file
            tSp = event_Time_on{k};
            strength = event_strength_on{k};
            spike_All = nan(length(delayTstart1),6);
            spik_swr_pow = nan(length(delayTstart1),6);
            strength_All = nan(length(delayTstart1),6);
            time_All = nan(length(delayTstart1),6);
            
            binCount = 0;
            for m = 1:length(delayTstart1)
%                 bin_Time_Idx = bin_Time>=Fig8DelayZonePos.delayPos1.startT(m) & bin_Time<delayTend1_2(m);
%                 bin_Time_Temp = bin_Time(bin_Time_Idx);
%                 for mm = 1:length(bin_Time_Temp)-1
%                     if sum(tSp>=bin_Time_Temp(mm) & tSp<bin_Time_Temp(mm+1)) > 0
%                         idx_temp = (tSp>=bin_Time_Temp(mm) & tSp<bin_Time_Temp(mm+1));
%                         pat_Delay_Event(k,mm+binCount) = 1;
%                         pat_Delay_Strength(k,mm+binCount) = strength(idx_temp);
%                     end
%                     
%                 end
%                 binCount = binCount + length(bin_Time_Temp);
                
                % return
                timeTemp = PathZone.posEndT.Return(m)-PathZone.posStartT.Return(m);
                time_All(m,1) = timeTemp;
                % base
                timeTemp = PathZone.posEndT.Base(m)-PathZone.posStartT.Base(m);
                time_All(m,2) = timeTemp;
                % delay
                timeTemp = delayTend1_2(m)-Fig8DelayZonePos.delayPos1.startT(m);
                time_All(m,3) = timeTemp;
                % stem
                timeTemp = PathZone.posEndT.Center(m)-(Fig8DelayZonePos.delayPos1.endT(m)+1/30);
                time_All(m,4) = timeTemp;
                % choice
                timeTemp = PathZone.posEndT.Choice(m)-PathZone.posStartT.Choice(m);
                time_All(m,5) = timeTemp;
                % reward
                timeTemp = PathZone.posEndT.Reward(m)-PathZone.posStartT.Reward(m);
                time_All(m,6) = timeTemp;
            
                % get spike count
                % return
                spikeInd = (tSp>=PathZone.posStartT.Return(m) & tSp<PathZone.posEndT.Return(m));
                spike_All(m,1) = sum(spikeInd);
                spik_swr_pow(m,1) = sum(event_swr_pow_on{k}(spikeInd));
                strength_All(m,1) = sum(strength(spikeInd));
                
                % base
                spikeInd = (tSp>=PathZone.posStartT.Base(m) & tSp<PathZone.posEndT.Base(m));
                spike_All(m,2) = sum(spikeInd);
                spik_swr_pow(m,2) = sum(event_swr_pow_on{k}(spikeInd));
                strength_All(m,2) = sum(strength(spikeInd));
                
                % delay
                spikeInd = (tSp>=Fig8DelayZonePos.delayPos1.startT(m) & tSp<delayTend1_2(m));
                spike_All(m,3) = sum(spikeInd);
                spik_swr_pow(m,3) = sum(event_swr_pow_on{k}(spikeInd));
                strength_All(m,3) = sum(strength(spikeInd));
                
                % stem
                spikeInd = (tSp>=(Fig8DelayZonePos.delayPos1.endT(m)+1/30) & tSp<PathZone.posEndT.Center(m));
                spike_All(m,4) = sum(spikeInd);
                spik_swr_pow(m,4) = sum(event_swr_pow_on{k}(spikeInd));
                strength_All(m,4) = sum(strength(spikeInd));
                
                % choice
                spikeInd = (tSp>=PathZone.posStartT.Choice(m) & tSp<PathZone.posEndT.Choice(m));
                spike_All(m,5) = sum(spikeInd);
                spik_swr_pow(m,5) = sum(event_swr_pow_on{k}(spikeInd));
                strength_All(m,5) = sum(strength(spikeInd));
                
                % reward
                spikeInd = (tSp>=PathZone.posStartT.Reward(m) & tSp<PathZone.posEndT.Reward(m));
                spike_All(m,6) = sum(spikeInd);
                spik_swr_pow(m,6) = sum(event_swr_pow_on{k}(spikeInd));
                strength_All(m,6) = sum(strength(spikeInd));       
            end 
            
            % return
            rateReturn = nansum(nansum(spike_All(:,1)))./sum(sum(time_All(:,1)));
            % base
            rateBase = nansum(nansum(spike_All(:,2)))./sum(sum(time_All(:,2)));
            % delay zone
            rateDelay = nansum(nansum(spike_All(:,3)))./sum(sum(time_All(:,3)));
            % stem zone
            rateStem = nansum(nansum(spike_All(:,4)))./sum(sum(time_All(:,4)));
            % choice zone
            rateChoice = nansum(nansum(spike_All(:,5)))./sum(sum(time_All(:,5)));
            % reward zone
            rateReward = nansum(nansum(spike_All(:,6)))./sum(sum(time_All(:,6)));
            
            Fig8_Assembly_OnOff_SWR_Rate.on.(sessDirs{j}).Assembly.rateReturn(k,:) = [rateReturn,sum(sum(time_All(:,1))),nansum(nansum(spike_All(:,1))),sum(spik_swr_pow(:,1)),sum(strength_All(:,1))];
            Fig8_Assembly_OnOff_SWR_Rate.on.(sessDirs{j}).Assembly.rateBase(k,:) = [rateBase,sum(sum(time_All(:,2))),nansum(nansum(spike_All(:,2))),sum(spik_swr_pow(:,2)),sum(strength_All(:,2))];
            Fig8_Assembly_OnOff_SWR_Rate.on.(sessDirs{j}).Assembly.rateDelay(k,:) = [rateDelay,sum(sum(time_All(:,3))),nansum(nansum(spike_All(:,3))),sum(spik_swr_pow(:,3)),sum(strength_All(:,3))];
            Fig8_Assembly_OnOff_SWR_Rate.on.(sessDirs{j}).Assembly.rateStem(k,:) = [rateStem,sum(sum(time_All(:,4))),nansum(nansum(spike_All(:,4))),sum(spik_swr_pow(:,4)),sum(strength_All(:,4))];
            Fig8_Assembly_OnOff_SWR_Rate.on.(sessDirs{j}).Assembly.rateChoice(k,:) = [rateChoice,sum(sum(time_All(:,5))),nansum(nansum(spike_All(:,5))),sum(spik_swr_pow(:,5)),sum(strength_All(:,5))];
            Fig8_Assembly_OnOff_SWR_Rate.on.(sessDirs{j}).Assembly.rateReward(k,:) = [rateReward,sum(sum(time_All(:,6))),nansum(nansum(spike_All(:,6))),sum(spik_swr_pow(:,6)),sum(strength_All(:,6))];
        
        end
        
%         pat_Delay_Event = zeros(patNum_off,length(bin_Time_SWR));
%         pat_Delay_Strength = zeros(patNum_off,length(bin_Time_SWR));
        for k = 1:patNum_off
            % get event time and strength from the assembly file
            tSp = event_Time_off{k};
            strength = event_strength_off{k};
            spike_All = nan(length(delayTstart1),6);
            spik_swr_pow = nan(length(delayTstart1),6);
            strength_All = nan(length(delayTstart1),6);
            time_All = nan(length(delayTstart1),6);
            
            binCount = 0;
            for m = 1:length(delayTstart1)
%                 bin_Time_Idx = bin_Time>=Fig8DelayZonePos.delayPos1.startT(m) & bin_Time<delayTend1_2(m);
%                 bin_Time_Temp = bin_Time(bin_Time_Idx);
%                 for mm = 1:length(bin_Time_Temp)-1
%                     if sum(tSp>=bin_Time_Temp(mm) & tSp<bin_Time_Temp(mm+1)) > 0
%                         idx_temp = (tSp>=bin_Time_Temp(mm) & tSp<bin_Time_Temp(mm+1));
%                         pat_Delay_Event(k,mm+binCount) = 1;
%                         pat_Delay_Strength(k,mm+binCount) = strength(idx_temp);
%                     end
%                     
%                 end
%                 binCount = binCount + length(bin_Time_Temp);
                
                % return
                timeTemp = PathZone.posEndT.Return(m)-PathZone.posStartT.Return(m);
                time_All(m,1) = timeTemp;
                % base
                timeTemp = PathZone.posEndT.Base(m)-PathZone.posStartT.Base(m);
                time_All(m,2) = timeTemp;
                % delay
                timeTemp = delayTend1_2(m)-Fig8DelayZonePos.delayPos1.startT(m);
                time_All(m,3) = timeTemp;
                % stem
                timeTemp = PathZone.posEndT.Center(m)-(Fig8DelayZonePos.delayPos1.endT(m)+1/30);
                time_All(m,4) = timeTemp;
                % choice
                timeTemp = PathZone.posEndT.Choice(m)-PathZone.posStartT.Choice(m);
                time_All(m,5) = timeTemp;
                % reward
                timeTemp = PathZone.posEndT.Reward(m)-PathZone.posStartT.Reward(m);
                time_All(m,6) = timeTemp;
            
                % get spike count
                % return
                spikeInd = (tSp>=PathZone.posStartT.Return(m) & tSp<PathZone.posEndT.Return(m));
                spike_All(m,1) = sum(spikeInd);
                spik_swr_pow(m,1) = sum(event_swr_pow_off{k}(spikeInd));
                strength_All(m,1) = sum(strength(spikeInd));
                
                % base
                spikeInd = (tSp>=PathZone.posStartT.Base(m) & tSp<PathZone.posEndT.Base(m));
                spike_All(m,2) = sum(spikeInd);
                spik_swr_pow(m,2) = sum(event_swr_pow_off{k}(spikeInd));
                strength_All(m,2) = sum(strength(spikeInd));
                
                % delay
                spikeInd = (tSp>=Fig8DelayZonePos.delayPos1.startT(m) & tSp<delayTend1_2(m));
                spike_All(m,3) = sum(spikeInd);
                spik_swr_pow(m,3) = sum(event_swr_pow_off{k}(spikeInd));
                strength_All(m,3) = sum(strength(spikeInd));
                
                % stem
                spikeInd = (tSp>=(Fig8DelayZonePos.delayPos1.endT(m)+1/30) & tSp<PathZone.posEndT.Center(m));
                spike_All(m,4) = sum(spikeInd);
                spik_swr_pow(m,4) = sum(event_swr_pow_off{k}(spikeInd));
                strength_All(m,4) = sum(strength(spikeInd));
                
                % choice
                spikeInd = (tSp>=PathZone.posStartT.Choice(m) & tSp<PathZone.posEndT.Choice(m));
                spike_All(m,5) = sum(spikeInd);
                spik_swr_pow(m,5) = sum(event_swr_pow_off{k}(spikeInd));
                strength_All(m,5) = sum(strength(spikeInd));
                
                % reward
                spikeInd = (tSp>=PathZone.posStartT.Reward(m) & tSp<PathZone.posEndT.Reward(m));
                spike_All(m,6) = sum(spikeInd);
                spik_swr_pow(m,6) = sum(event_swr_pow_off{k}(spikeInd));
                strength_All(m,6) = sum(strength(spikeInd));       
            end 
            
            % return
            rateReturn = nansum(nansum(spike_All(:,1)))./sum(sum(time_All(:,1)));
            % base
            rateBase = nansum(nansum(spike_All(:,2)))./sum(sum(time_All(:,2)));
            % delay zone
            rateDelay = nansum(nansum(spike_All(:,3)))./sum(sum(time_All(:,3)));
            % stem zone
            rateStem = nansum(nansum(spike_All(:,4)))./sum(sum(time_All(:,4)));
            % choice zone
            rateChoice = nansum(nansum(spike_All(:,5)))./sum(sum(time_All(:,5)));
            % reward zone
            rateReward = nansum(nansum(spike_All(:,6)))./sum(sum(time_All(:,6)));
            
            Fig8_Assembly_OnOff_SWR_Rate.off.(sessDirs{j}).Assembly.rateReturn(k,:) = [rateReturn,sum(sum(time_All(:,1))),nansum(nansum(spike_All(:,1))),sum(spik_swr_pow(:,1)),sum(strength_All(:,1))];
            Fig8_Assembly_OnOff_SWR_Rate.off.(sessDirs{j}).Assembly.rateBase(k,:) = [rateBase,sum(sum(time_All(:,2))),nansum(nansum(spike_All(:,2))),sum(spik_swr_pow(:,2)),sum(strength_All(:,2))];
            Fig8_Assembly_OnOff_SWR_Rate.off.(sessDirs{j}).Assembly.rateDelay(k,:) = [rateDelay,sum(sum(time_All(:,3))),nansum(nansum(spike_All(:,3))),sum(spik_swr_pow(:,3)),sum(strength_All(:,3))];
            Fig8_Assembly_OnOff_SWR_Rate.off.(sessDirs{j}).Assembly.rateStem(k,:) = [rateStem,sum(sum(time_All(:,4))),nansum(nansum(spike_All(:,4))),sum(spik_swr_pow(:,4)),sum(strength_All(:,4))];
            Fig8_Assembly_OnOff_SWR_Rate.off.(sessDirs{j}).Assembly.rateChoice(k,:) = [rateChoice,sum(sum(time_All(:,5))),nansum(nansum(spike_All(:,5))),sum(spik_swr_pow(:,5)),sum(strength_All(:,5))];
            Fig8_Assembly_OnOff_SWR_Rate.off.(sessDirs{j}).Assembly.rateReward(k,:) = [rateReward,sum(sum(time_All(:,6))),nansum(nansum(spike_All(:,6))),sum(spik_swr_pow(:,6)),sum(strength_All(:,6))];
        
        end
    end
    
    if p.writeToFile
        save(fullfile(savedir2,'Fig8_Assembly_OnOff_SWR_Rate-25ms.mat'), 'Fig8_Assembly_OnOff_SWR_Rate');
    end
    
    fprintf('Finished spike shape analysis for session %d\n',i);
    clear Fig8_Assembly_OnOff_SWR_Rate
end
end

function [Frequency, ADBitVolts] = ReadHeaderEEG(Header)
% extract ADBvolt for that EEG
j = 1;
while isempty(strfind(Header{j}, '-ADBitVolts'))
    j=j+1;
end

if j <= length(Header)
    bvInd = strfind(Header{j}, '-ADBitVolts');
    ADBitVolts = str2double(Header{j}(11+bvInd:end));
else
    sprintf('No ADBvolts detected')
    ADBitVolts = 6.1e-8;
end

% Extrac sampling frequency for EEG
k = 1;
while isempty(strfind(Header{k}, '-SamplingFrequency'))
    k=k+1;
end

if k <= length(Header)
    bvInd = strfind(Header{k}, '-SamplingFrequency');
    Frequency = str2double(Header{k}(18+bvInd:end));
else
    sprintf('No Frequency detected')
    Frequency = 2000;
end
end

function xf = fftbandpass(x,Fs,Fs1,Fp1,Fp2,Fs2)
% function XF = fftbandpass(X,FS,FS1,FP1,FP2,FS2)
%
% Bandpass filter for the signal X (time x trials). An acausal fft 
% algorithm is applied (i.e. no phase shift). The filter functions is         
% constructed from a Hamming window. 
%
% Fs : sampling frequency
%
% The passbands (Fp1 Fp2) and stop bands (Fs1 Fs2) are defined as
%                 -----------                      
%                /           \
%               /             \
%              /               \
%             /                 \
%   ----------                   ----------------- 
%           Fs1  Fp1       Fp2  Fs2              
%
% If no output arguments are assigned the filter function H(f) and
% impulse response are plotted. 
%
% NOTE: for long data traces the filter is very slow.
%
%------------------------------------------------------------------------
% Ole Jensen, Brain Resarch Unit, Low Temperature Laboratory,
% Helsinki University of Technology, 02015 HUT, Finland,
% Report bugs to ojensen@neuro.hut.fi
%------------------------------------------------------------------------

%    Copyright (C) 2000 by Ole Jensen 
%    This program is free software; you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation; either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You can find a copy of the GNU General Public License
%    along with this package (4DToolbox); if not, write to the Free Software
%    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

if size(x,1) == 1
    x = x';
end
% Make x even
Norig = size(x,1); 
if rem(Norig,2)
    x = [x' zeros(size(x,2),1)]';                
end

% Normalize frequencies  
Ns1 = Fs1/(Fs/2);
Ns2 = Fs2/(Fs/2);
Np1 = Fp1/(Fs/2);
Np2 = Fp2/(Fs/2);

% Construct the filter function H(f)
N = size(x,1);
Nh = N/2;

B = fir2(N-1,[0 Ns1 Np1 Np2 Ns2 1],[0 0 1 1 0 0]); 
H = abs(fft(B));  % Make zero-phase filter function
IPR = real(ifft(H));
if nargout == 0 
    subplot(2,1,1)
    f = Fs*(0:Nh-1)/(N);
    plot(f,H(1:Nh));
    xlim([0 2*Fs2])
    ylim([0 1]); 
    title('Filter function H(f)')
    xlabel('Frequency (Hz)')
    subplot(2,1,2)
    plot((1:Nh)/Fs,IPR(1:Nh))
    xlim([0 2/Fp1])
    xlabel('Time (sec)')
    ylim([min(IPR) max(IPR)])
    title('Impulse response')
end


if size(x,2) > 1
    for k=1:size(x,2)
        xf(:,k) = real(ifft(fft(x(:,k)) .* H'));
    end
    xf = xf(1:Norig,:);
else
    xf = real(ifft(fft(x') .* H));
    xf = xf(1:Norig);
end
end

function [Phase, Amp] = assignPhase2(Eeg)

% Hilbert transform
Z = hilbert(Eeg);

% Wave amplitude
Amp = abs(Z);

% EEG phase in rad
Phase = angle(Z);

% Rad to Degree (-180 to +180)
Phase = Phase / pi *180;

% Degree (0 to +360)
Phase = Phase + 180;
end

% This code detects oscillation events based on evelop graph of filtered EEG
% if one cycle of the envelop is longer than defined Time and
% peak if over threshold, this cycle is treated as an event episode
% Input: LFP, filtered LFP in the interested frequency range
% Fs: Sampling frequency
% ThresSD: threshold of signal envelop being included (mean+ThresSD*SD) 
% TJoingtLim: maximum time between two events to be concatenated as one.(unit: sec)
% TLthLim: Minimum time length for an event episode.(unit: sec)
% PLOT: plot raw signal and chosen events episodes for visual check
% Output: evtOnIdx: Event onset inidx in filtered signal
% evtOffIdx: Event offset inidx in filtered signal
% eegLabel: filtered signal index in the events are labeled as 1, else 0
%
% running example: [evtOnIdx,evtOffIdx,eegLabel] = eventFinder(EEG,EEGt,1,0.5,0.2,1,1);
% Li, 07-Feb-2020
function [evtOnIdx,evtOffIdx,eegLabel] = eventFinder_SWR(LFP_Raw,LFP_Fil,Fs,ThresSD,TJoingtLim,TLthLim,PLOT)
if nargin == 2
    ThresSD = 1;
    TJoingtLim = 0.2;
    PLOT = 0;
elseif nargin == 2
    TJoingtLim = 0.2;
    PLOT = 0;
elseif nargin == 4
    PLOT = 0;
end

% Direction of vector
if size(LFP_Fil,1)>1
    LFP_Fil=LFP_Fil';
end

eegLabel = zeros(length(LFP_Fil),1);

% % get the upper boundry of filtered LFP signal envlop
% % Hilbert transform
% Z = hilbert(LFP_Fil);
% % Wave amplitude
% EnvUpper = abs(Z);
[swr_peaks,swr_peakIdx] = findpeaks(LFP_Fil);
% get mean and SD of signal envlop
swr_mean = nanmean(swr_peaks);
swr_SD = std(swr_peaks,'omitnan');

% calculate threshold value for defining the signal episode
swr_thres = swr_mean+ThresSD*swr_SD;

% get peak and trough index of envelop (where later start and end of an episode
% is located on)
% envelop peak over threshold is labeled as 1
% envelop peak below threshold is labeled as 0
% envelop trough is labeled as -1

% the peak of envelop
PEAKidx = find(swr_peaks >= swr_thres); % PEAK which over thre1 idx on peakidx
peakLabel = zeros(1,length(swr_peaks));
peakLabel(PEAKidx) = 1;
% to calculate envelop trough index
maxEnv = max(LFP_Fil);
InverseEnv = maxEnv-LFP_Fil;
[swrTrough,swr_TroughIdx] = findpeaks(InverseEnv);
troughLabel = zeros(1,length(swr_TroughIdx))-1;

% combine peak and trough index into one vector
peakTroughIdx = [swr_peakIdx,swr_TroughIdx];
peakTroughLabel = [peakLabel,troughLabel];
[peakTroughIdx,sortIdx] = sort(peakTroughIdx);
peakTroughLabel = peakTroughLabel(sortIdx);

% start detecting event
StartIdx = [];
EndIdx = [];
% find start of the over threshold peak and end of the overthreshold
% peak
for i = 1:length(peakTroughLabel)-2
    % event start is from an envelop trough
    if peakTroughLabel(i) == 0&& peakTroughLabel(i+1)==-1 && peakTroughLabel(i+2)==1
        StartIdx = [StartIdx,i+1];
        % event end is an envelop trough
    elseif peakTroughLabel(i) == 1 && peakTroughLabel(i+1)==-1 && peakTroughLabel(i+2)==0
        EndIdx = [EndIdx,i+1];
    end
end

% exclude special cases at boundries
if EndIdx(1) < StartIdx(1)
    StartIdx = [1,StartIdx];
end
if StartIdx(end) > EndIdx(end)
    EndIdx = [EndIdx,length(peakTroughLabel)];
end

eventStartIdx = peakTroughIdx(StartIdx);
eventEndIdx = peakTroughIdx(EndIdx);

Label1 = zeros(1,length(LFP_Fil));
Label1(eventStartIdx) = 1;


% if there are next start is within TimeJointLim sec, join two event as one
% event
nEvent = 0;
breakLength = TJoingtLim*Fs;
i = 1;
while i <=length(eventStartIdx)
    Start = eventStartIdx(i);
    Ending = eventEndIdx(i);
    nEvent = nEvent+1;
     
    if Ending+breakLength <= length(LFP_Fil)
        step = breakLength;
    else
        step = length(LFP_Fil)-Ending;
    end
    
    % check whithin the timelimit, if there is new start
    while any(Label1(Ending+1:Ending+step)) == 1        
        if i+1 <= length(eventStartIdx) 
            i = i+1;
        end
        Ending = eventEndIdx(i);
        
        % redeine step again incase i is close to end of sequence
        if Ending+breakLength <=length(LFP_Fil)
            step = breakLength;
        else
            step = length(LFP_Fil)-Ending;
        end        
    end       
    i=i+1;
    evtOnIdx(nEvent) = Start;
    evtOffIdx(nEvent) = Ending;
    
end

% if from an event start to end is < timelength set, exclude this event
k=zeros(1,length(evtOnIdx));
for j = 1:length(evtOnIdx)
    if evtOffIdx(j)-evtOnIdx(j) >= TLthLim*Fs
        k(j) = 1;
    end
end
evtOnIdx = evtOnIdx(k==1);
evtOffIdx = evtOffIdx(k==1);
nEvent = sum(k==1);

% label eeg index inside event as 1
for m = 1:nEvent
    eegLabel(evtOnIdx(m):evtOffIdx(m)) = 1;
end

% plot figures to see ripples
if PLOT
     % plot original EEG to see singal quality
    figure
    plot(LFP_Raw+1000)   
    hold on
    plot(LFP_Fil,'b')
    plot(swr_peakIdx(PEAKidx),LFP_Fil(swr_peakIdx((PEAKidx))),'r^')
%     plot(swr_TroughIdx,LFP_Fil(swr_TroughIdx),'g.')
    
    for i = 1:nEvent
        plot(evtOnIdx(i):evtOffIdx(i),LFP_Fil(evtOnIdx(i):evtOffIdx(i))-500)
    end
end
ylim('auto')
hold off
end
