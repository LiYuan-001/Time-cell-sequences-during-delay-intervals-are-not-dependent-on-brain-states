% Theta-Delta ratio was calculated based on band power (wavelet) in each
% session
% Theta-Delta ratio can be used to determine whether there is theta
% oscillation
% Li Yuan, 08-Mar-2021, UCSD
function Fig8TreadmillEEGThetaDeltaRatio(inFile,AnalyzeSes)

% set parameter
p.delta = [1 4];
p.theta = [6 12];


p.space = 'lin'; % 'log' or 'lin'  spacing of f's
p.frange = [1 100];
p.nCycles = 10; % cycles for morelet wavelet
p.dt = 0.005; 

p.smoothTime = 0.15;

p.savePlot = 1;
p.writeToFile = 1;
% Read in input information
sessInfo = SessInfoImport(inFile);  

% parameters for import EEG from neuralynx / matlab utilities
FieldSelectionFlags = [1 0 0 0 1];
HeaderExtractionFlag = 1;
ExtractionMode = 1;
ExtractionModeVector = [];
    
for i = AnalyzeSes(1:end)
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\EEGSpec\ThetaDeltaRatio');
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
    ThetaDeltaRatio.rat = sessInfo(i).animal;
    ThetaDeltaRatio.day = sessInfo(i).day;
    ThetaDeltaRatio.smoothTime = p.smoothTime;
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    sessDirs = sessInfo(i).sessDirs;
    
    % load eeg file
    eegInd = sprintf('%s%d%s', 'CSC',sessInfo(i).EEGch,'.ncs');
    cscFile = fullfile(mainDir,eegInd);
    % extract EEG by Nlx import
    [Timestamps,Samples,Header] = Nlx2MatCSC(cscFile, FieldSelectionFlags,HeaderExtractionFlag, ExtractionMode, ExtractionModeVector);
    p.t1 = Timestamps(1);
    
    % reshape EEG samples
    EEG = reshape(Samples,length(Samples(:)),1);
    % get Frequency and ADBvolts for this channel
    [p.Fs, p.ADBVolts] = ReadHeaderEEG(Header);
    % reasample timestamps;
    % unit: sec
    EEGTs = timeStampResample(Samples,Timestamps,p.Fs)./10^6;   
    % get EEG to uV
    EEG = EEG*p.ADBVolts*10^6;
    
    for j = 1:length(sessDirs)
        close all
        % load analyzed positions
        posFile = fullfile(mainDir,sessDirs{j}, 'locInfo.mat');
        locInfo = load(posFile);
        pathFile = fullfile(mainDir,sessDirs{j}, 'pathData.mat');
        pathData = load(pathFile);
        pathZoneFile = fullfile(mainDir,sessDirs{j}, 'PathZone.mat');
        load(pathZoneFile);
        delayFile = fullfile(mainDir,sessDirs{j}, 'Fig8DelayZonePos.mat');
        load(delayFile);
        
        
        [minstart,EEGStartIdx] =  min(abs(EEGTs-pathData.t(1)));
        [minend,EEGEndIdx] = min(abs(EEGTs-pathData.t(end)));
        EEG_Session = EEG(EEGStartIdx:EEGEndIdx);
        % double check timestamp
        if any([minstart,minend] > (EEGTs(2)-EEGTs(1)))
            error('time is wrong')
        end
        
        lfp.data = EEG_Session;
        lfp.timestamps = EEGTs(EEGStartIdx:EEGEndIdx);
        lfp.samplingRate = p.Fs;
        
        spec = PowerSpectrum_Wavelet(lfp,p.nCycles,p.dt,p.frange,'nfreqs',100);
        
        % double check timestamp
        if any([min(spec.timestamps - pathData.t(1)),min(spec.timestamps - pathData.t(end))] > (spec.timestamps(2)-spec.timestamps(1)))
            error('time is wrong')
        end
        % N time points power
        power_dt = bandpower(abs(spec.data)',spec.freqs,p.delta,'psd');
        power_th = bandpower(abs(spec.data)',spec.freqs,p.theta,'psd');
        % tdr: theta-delta ratio
        tdr = power_th./power_dt;
        % smooth tdr with 0.1 sec before and after the current timepoint
        smthWindow = p.smoothTime/p.dt;
        tdr_Smooth = tdr;
        for m = smthWindow+1:length(tdr)-smthWindow
            tdr_Smooth(m) = nanmean(tdr(m-smthWindow:m+smthWindow));
        end
        % Periodogram of whole session
        h1=figure(1);
        h1.Position = [100 100 1200 900];

        % Spectrogram of whole session
        h(1) = subplot(2,1,1);
        imagesc(spec.timestamps,log2(spec.freqs),spec.amp')
        SpecColorRange( spec.amp );
        colormap jet
        LogScale('y',2)
        ylabel('f (Hz)')
        xlabel('time (sec)')
        axis xy
        TITLE1 = sprintf('%s%d%s%d%s','LFP Power -',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',sessDirs{j});
        TITLE2 = 'Spectrogram';
        title({TITLE1,TITLE2},'Interpreter','None')
         
        h(2) = subplot(2,1,2);
        plot(spec.timestamps,tdr_Smooth)
        ylabel('Theta / Delta power')
        xlabel('time (sec)')
        title('Theta-Delta Ratio')
% 
%         h(3) = subplot(3,1,3);
%         plot(lfp.timestamps,lfp.data);
        
        linkaxes(h,'x');
         
        if p.savePlot
            figure(1)
            figName = sprintf('%s%s%d%s%d%s%s%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',sessDirs{j},'-ThetaDeltaRatio');
            print(figName,'-dpng','-r300');
        end
        
        ThetaDeltaRatio.(sessDirs{j}).timeStamp = spec.timestamps;
        ThetaDeltaRatio.(sessDirs{j}).thetaDeltaRatio = tdr_Smooth;
        ThetaDeltaRatio.(sessDirs{j}).thetaPow = power_th;
        ThetaDeltaRatio.(sessDirs{j}).deltaPow = power_dt;

    end
    if p.writeToFile
        save(fullfile(savedir2,'ThetaDeltaRatio.mat'), 'ThetaDeltaRatio');
    end
    fprintf('Finished theta delta ratio analysis for session %d\n',i);
    clear ThetaDeltaRatio
end
end