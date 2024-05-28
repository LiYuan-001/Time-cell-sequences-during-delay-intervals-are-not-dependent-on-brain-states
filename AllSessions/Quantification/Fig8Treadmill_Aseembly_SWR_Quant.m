function Fig8Treadmill_Aseembly_SWR_Quant(inFile,AnalyzeSes)

close all

p.savePlot = 0;
p.writeToFile = 0;

allAssemblyNum = 0;

% % Read in input information
sessInfo = SessInfoImport(inFile);
% get each phase names (no delay etc)
sessDirs = {'on10','off10','on30','off30'};
    
if p.savePlot
    % directory for plot figures
    % generate a folder for each rat eah day under the current folder
    savedir = sprintf('%s%s',cd,'\Figures\Quantification Figures\Time cell num-dist quant');
    if ~exist(savedir, 'dir')
        mkdir(savedir);
    end
end

armNames = {'rateReturn','rateDelay','rateStem','rateChoice','rateReward'};
session = {'on','off'};
timeNames = {'timeReturn','timeDelay','timeStem','timeChoice','timeReward'};
                
for j = 1:length(session)
    for k = 1:length(armNames)
        % arm compare value
        swr_Rate.(session{j}).(armNames{k}) = [];
        assembly_Rate.(session{j}).(armNames{k}) = [];
        assembly_strength.(session{j}).(armNames{k}) = [];
        assembly_swr_pow.(session{j}).(armNames{k}) = [];
    end
end


for i = AnalyzeSes(1:end)
%     % load average firing rate file
%     SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
%     load(SpikeShapeFile);
    % load arm rate file
    armRateFile = fullfile(sessInfo(i).mainDir,'Cell Property', 'Fig8_Assembly_SWR_Rate-25ms.mat');
    load(armRateFile);
    
    assemblyNum = size(Fig8_Assembly_SWR_Rate.on10_1.Assembly.rateReturn,1);
    allAssemblyNum = allAssemblyNum + assemblyNum;
    
    
    for k = 1:length(armNames)
        
        %% swr
        spikeNumTemp = (Fig8_Assembly_SWR_Rate.on10_1.SWR.(armNames{k})(3) + Fig8_Assembly_SWR_Rate.on30_1.SWR.(armNames{k})(3)+...
            Fig8_Assembly_SWR_Rate.on10_2.SWR.(armNames{k})(3) + Fig8_Assembly_SWR_Rate.on30_2.SWR.(armNames{k})(3));
        timeTemp =(Fig8_Assembly_SWR_Rate.on10_1.SWR.(armNames{k})(2) + Fig8_Assembly_SWR_Rate.on30_1.SWR.(armNames{k})(2)+...
            Fig8_Assembly_SWR_Rate.on10_2.SWR.(armNames{k})(2) + Fig8_Assembly_SWR_Rate.on30_2.SWR.(armNames{k})(2));
        swr_Rate.on.(armNames{k}) = [swr_Rate.on.(armNames{k});spikeNumTemp./timeTemp];
        
        spikeNumTemp = (Fig8_Assembly_SWR_Rate.off10_1.SWR.(armNames{k})(3) + Fig8_Assembly_SWR_Rate.off30_1.SWR.(armNames{k})(3)+...
            Fig8_Assembly_SWR_Rate.off10_2.SWR.(armNames{k})(3) + Fig8_Assembly_SWR_Rate.off30_2.SWR.(armNames{k})(3));
        timeTemp =(Fig8_Assembly_SWR_Rate.off10_1.SWR.(armNames{k})(2) + Fig8_Assembly_SWR_Rate.off30_1.SWR.(armNames{k})(2)+...
            Fig8_Assembly_SWR_Rate.off10_2.SWR.(armNames{k})(2) + Fig8_Assembly_SWR_Rate.off30_2.SWR.(armNames{k})(2));
        swr_Rate.off.(armNames{k}) = [swr_Rate.off.(armNames{k});spikeNumTemp./timeTemp];
        
        %% assembly
        spikeNumTemp = (Fig8_Assembly_SWR_Rate.on10_1.Assembly.(armNames{k})(:,3) + Fig8_Assembly_SWR_Rate.on30_1.Assembly.(armNames{k})(:,3)+...
            Fig8_Assembly_SWR_Rate.on10_2.Assembly.(armNames{k})(:,3) + Fig8_Assembly_SWR_Rate.on30_2.Assembly.(armNames{k})(:,3));
        timeTemp =(Fig8_Assembly_SWR_Rate.on10_1.Assembly.(armNames{k})(:,2) + Fig8_Assembly_SWR_Rate.on30_1.Assembly.(armNames{k})(:,2)+...
            Fig8_Assembly_SWR_Rate.on10_2.Assembly.(armNames{k})(:,2) + Fig8_Assembly_SWR_Rate.on30_2.Assembly.(armNames{k})(:,2));
        assembly_Rate.on.(armNames{k}) = [assembly_Rate.on.(armNames{k});spikeNumTemp./timeTemp];
        
        strengthTemp = (Fig8_Assembly_SWR_Rate.on10_1.Assembly.(armNames{k})(:,5) + Fig8_Assembly_SWR_Rate.on30_1.Assembly.(armNames{k})(:,5)+...
            Fig8_Assembly_SWR_Rate.on10_2.Assembly.(armNames{k})(:,5) + Fig8_Assembly_SWR_Rate.on30_2.Assembly.(armNames{k})(:,5));
        assembly_strength.on.(armNames{k}) = [assembly_strength.on.(armNames{k});strengthTemp./spikeNumTemp];
        
        swr_pow_Temp = (Fig8_Assembly_SWR_Rate.on10_1.Assembly.(armNames{k})(:,4) + Fig8_Assembly_SWR_Rate.on30_1.Assembly.(armNames{k})(:,4)+...
            Fig8_Assembly_SWR_Rate.on10_2.Assembly.(armNames{k})(:,4) + Fig8_Assembly_SWR_Rate.on30_2.Assembly.(armNames{k})(:,4));
        assembly_swr_pow.on.(armNames{k}) = [assembly_swr_pow.on.(armNames{k});swr_pow_Temp./spikeNumTemp];
        
        
        spikeNumTemp = (Fig8_Assembly_SWR_Rate.off10_1.Assembly.(armNames{k})(:,3) + Fig8_Assembly_SWR_Rate.off30_1.Assembly.(armNames{k})(:,3)+...
            Fig8_Assembly_SWR_Rate.off10_2.Assembly.(armNames{k})(:,3) + Fig8_Assembly_SWR_Rate.off30_2.Assembly.(armNames{k})(:,3));
        timeTemp =(Fig8_Assembly_SWR_Rate.off10_1.Assembly.(armNames{k})(:,2) + Fig8_Assembly_SWR_Rate.off30_1.Assembly.(armNames{k})(:,2)+...
            Fig8_Assembly_SWR_Rate.off10_2.Assembly.(armNames{k})(:,2) + Fig8_Assembly_SWR_Rate.off30_2.Assembly.(armNames{k})(:,2));
        assembly_Rate.off.(armNames{k}) = [assembly_Rate.off.(armNames{k});spikeNumTemp./timeTemp];
        
        strengthTemp = (Fig8_Assembly_SWR_Rate.off10_1.Assembly.(armNames{k})(:,5) + Fig8_Assembly_SWR_Rate.off30_1.Assembly.(armNames{k})(:,5)+...
            Fig8_Assembly_SWR_Rate.off10_2.Assembly.(armNames{k})(:,5) + Fig8_Assembly_SWR_Rate.off30_2.Assembly.(armNames{k})(:,5));
        assembly_strength.off.(armNames{k}) = [assembly_strength.off.(armNames{k});strengthTemp./spikeNumTemp];
        
        swr_pow_Temp = (Fig8_Assembly_SWR_Rate.off10_1.Assembly.(armNames{k})(:,4) + Fig8_Assembly_SWR_Rate.off30_1.Assembly.(armNames{k})(:,4)+...
            Fig8_Assembly_SWR_Rate.off10_2.Assembly.(armNames{k})(:,4) + Fig8_Assembly_SWR_Rate.off30_2.Assembly.(armNames{k})(:,4));
        assembly_swr_pow.off.(armNames{k}) = [assembly_swr_pow.off.(armNames{k});swr_pow_Temp./spikeNumTemp];
        
    end    
end
% pyrIdx = avgRate>0.1 & avgRate<5;
% 
% %% delay definition 1: delay aligned by barrier
% % as long as it is over average rate threshold in one session
% delay_onIdx = sum(rate_Delay>p.avgRateThres,1) & pyrIdx;
% 
% % delay active cell rate difference in treadmill on sessions

for k = 1:length(armNames)
    
    figure(1)
    xplot = zeros(2,length(swr_Rate.on.(armNames{k})));
    xplot(1,:) = 3*k-1.7 + xplot(1,:);
    xplot(2,:) = 3*k-1.3 + xplot(2,:);
    plot(xplot,[swr_Rate.on.(armNames{k})';swr_Rate.off.(armNames{k})'],'o-','Color',[0.7,0.7,0.7])
    hold on
    
    Violin(swr_Rate.on.(armNames{k}),3*k-2,'ShowData',false,'ViolinColor',[1,0,0])
    Violin(swr_Rate.off.(armNames{k}),3*k-1,'ShowData',false,'ViolinColor',[0,0,0])
    set(gca,'XTick',[1:3:15],'XTickLabel',armNames);
    xlim([0 15])
    title('SWR Rate')
    
    figure(2)
    Violin(swr_Rate.on.(armNames{k})-swr_Rate.off.(armNames{k}),k,'ShowData',false,'ViolinColor',[0,0,1])
    set(gca,'XTick',[1:5],'XTickLabel',armNames)
    xlim([0 6])
    title('SWR Rate')

    figure(3)
    xplot = zeros(2,length(assembly_Rate.on.(armNames{k})));
    xplot(1,:) = 3*k-1.7 + xplot(1,:);
    xplot(2,:) = 3*k-1.3 + xplot(2,:);
    plot(xplot,[assembly_Rate.on.(armNames{k})';assembly_Rate.off.(armNames{k})'],'o-','Color',[0.7,0.7,0.7])
    hold on
    
    boxplot(assembly_Rate.on.(armNames{k}),'Position',3*k-2,'Color',[1,0,0])
    hold on
    boxplot(assembly_Rate.off.(armNames{k}),'Position',3*k-1,'Color',[0,0,0])
    set(gca,'XTick',[1:3:15],'XTickLabel',armNames);
    xlim([0 15])
    title('Assembly Rate')
    
    figure(4)
    Violin(assembly_Rate.on.(armNames{k})-assembly_Rate.off.(armNames{k}),k,'ShowData',false,'ViolinColor',[0,0,1])
    set(gca,'XTick',[1:5],'XTickLabel',armNames)
    xlim([0 6])
    title('Assembly Rate')
    
    figure(5)
    xplot = zeros(2,length(assembly_strength.on.(armNames{k})));
    xplot(1,:) = 3*k-1.7 + xplot(1,:);
    xplot(2,:) = 3*k-1.3 + xplot(2,:);
    plot(xplot,[assembly_strength.on.(armNames{k})';assembly_strength.off.(armNames{k})'],'o-','Color',[0.7,0.7,0.7])
    hold on
    
    boxplot(assembly_strength.on.(armNames{k}),'Position',3*k-2,'Color',[1,0,0])
    hold on
    boxplot(assembly_strength.off.(armNames{k}),'Position',3*k-1,'Color',[0,0,0])
    
    set(gca,'XTick',[1:3:15],'XTickLabel',armNames);
    xlim([0 15])
    title('Assembly strength')
    
    figure(6)
    Violin(assembly_strength.on.(armNames{k})-assembly_strength.off.(armNames{k}),k,'ShowData',false,'ViolinColor',[0,0,1])
    set(gca,'XTick',[1:5],'XTickLabel',armNames)
    xlim([0 6])
    title('Assembly strength')
    
    figure(7)
    xplot = zeros(2,length(assembly_swr_pow.on.(armNames{k})));
    xplot(1,:) = 3*k-1.7 + xplot(1,:);
    xplot(2,:) = 3*k-1.3 + xplot(2,:);
    plot(xplot,10*log10([assembly_swr_pow.on.(armNames{k})';assembly_swr_pow.off.(armNames{k})']),'o-','Color',[0.7,0.7,0.7])
    hold on
    
    Violin(10*log10(assembly_swr_pow.on.(armNames{k})),3*k-2,'ShowData',false,'ViolinColor',[1,0,0])
    Violin(10*log10(assembly_swr_pow.off.(armNames{k})),3*k-1,'ShowData',false,'ViolinColor',[0,0,0])
    set(gca,'XTick',[1:3:15],'XTickLabel',armNames);
    xlim([0 15])
    ylabel('Power(db)')
    title('Assembly swr power')
    
    figure(8)
    Violin(10*log10(assembly_swr_pow.on.(armNames{k}))-10*log10(assembly_swr_pow.off.(armNames{k})),k,'ShowData',false,'ViolinColor',[0,0,1])
    set(gca,'XTick',[1:5],'XTickLabel',armNames)
    xlim([0 6])
    ylabel('Power(db)')
    title('Assembly swr power')
%     neuronRate_mean.on.(armNames{k}) = mean(swr_Rate.on.(armNames{k})(pyrIdx));
%     neuronRate_mean.off.(armNames{k}) = mean(swr_Rate.off.(armNames{k})(pyrIdx));
%     neuronRate_err.on.(armNames{k}) = std(swr_Rate.on.(armNames{k})(pyrIdx))./sqrt(sum(pyrIdx));
%     neuronRate_err.off.(armNames{k}) = std(swr_Rate.off.(armNames{k})(pyrIdx))./sqrt(sum(pyrIdx));
%     
%     burstRate_mean.on.(armNames{k}) = mean(assembly_Rate.on.(armNames{k})(pyrIdx));
%     burstRate_mean.off.(armNames{k}) = mean(assembly_Rate.off.(armNames{k})(pyrIdx));
%     burstRate_err.on.(armNames{k}) = std(assembly_Rate.on.(armNames{k})(pyrIdx))./sqrt(sum(pyrIdx));
%     burstRate_err.off.(armNames{k}) = std(assembly_Rate.off.(armNames{k})(pyrIdx))./sqrt(sum(pyrIdx));
%     
%     singleSpkRate_mean.on.(armNames{k}) = mean(assembly_strength.on.(armNames{k})(pyrIdx));
%     singleSpkRate_mean.off.(armNames{k}) = mean(assembly_strength.off.(armNames{k})(pyrIdx));
%     singleSpkRate_err.on.(armNames{k}) = std(assembly_strength.on.(armNames{k})(pyrIdx))./sqrt(sum(pyrIdx));
%     singleSpkRate_err.off.(armNames{k}) = std(assembly_strength.off.(armNames{k})(pyrIdx))./sqrt(sum(pyrIdx));
% %     
%     AssemblyRate_mean.on.(armNames{k}) = mean(AssemblyRate.on.(armNames{k}));
%     AssemblyRate_mean.off.(armNames{k}) = mean(AssemblyRate.off.(armNames{k}));
%     AssemblyRate_err.on.(armNames{k}) = std(AssemblyRate.on.(armNames{k}))./sqrt(length(AssemblyRate.on.(armNames{k})));
%     AssemblyRate_err.off.(armNames{k}) = std(AssemblyRate.off.(armNames{k}))./sqrt(length(AssemblyRate.off.(armNames{k})));
    
%     figure(1)
%     Violin(log2(neuronRate.on.(armNames{k})(pyrIdx)),2*k-1,'ViolinColor',[1,0,0],'ShowData',false)
%     hold on
%     Violin(log2(neuronRate.off.(armNames{k})(pyrIdx)),2*k,'ViolinColor',[0,0,0],'ShowData',false)
%     set(gca,'XTick',[1:2:10],'XTickLabel',armNames);
%     xlim([0 10])
%     title('Pyr Rate')
%     
    
%     figure(1)
%     errorbar(k,neuronRate_mean.on.(armNames{k}),neuronRate_err.on.(armNames{k}),'rd');
%     hold on
%     errorbar(k+0.2,neuronRate_mean.off.(armNames{k}),neuronRate_err.off.(armNames{k}),'kd');
%     title('Pyr Rate')
%     set(gca,'XTick',[1:5],'XTickLabel',armNames);
%     xlim([0 6])
%     
%     
%     figure(2)
%     errorbar(k,burstRate_mean.on.(armNames{k}),burstRate_err.on.(armNames{k}),'rd');
%     hold on
%     errorbar(k+0.2,burstRate_mean.off.(armNames{k}),burstRate_err.off.(armNames{k}),'kd');
%     title('Burst rate')
%     set(gca,'XTick',[1:5],'XTickLabel',armNames);
%     xlim([0 6])
%     
%     
%     figure(3)
%     errorbar(k,singleSpkRate_mean.on.(armNames{k}),singleSpkRate_err.on.(armNames{k}),'rd');
%     hold on
%     errorbar(k+0.2,singleSpkRate_mean.off.(armNames{k}),singleSpkRate_err.off.(armNames{k}),'kd');
%     title('Single spike rate')
%     set(gca,'XTick',[1:5],'XTickLabel',armNames);
%     xlim([0 6])
% %     
%     figure(4)
%     errorbar(k,AssemblyRate_mean.on.(armNames{k}),AssemblyRate_err.on.(armNames{k}),'rd');
%     hold on
%     errorbar(k+0.2,AssemblyRate_mean.off.(armNames{k}),AssemblyRate_err.off.(armNames{k}),'kd');
%     title('Assembly rate')
end

% set(gca,'XTick',[1:5],'XTickLabel',armNames);
% xlim([0 6])
% delay active cells rate difference in treadmill off sessions

    
end