% time field in this code is firstly defined by whether the within session
% shuffled stability passed 5 percentile
% stability is definsed as median value of all trial by trial correlation
% in one session
% then extract field by peak rate
function Fig8TreadmillDelayField_TrialbyTrialShuffle_PartDelay(inFile,AnalyzeSes)
close all

p.savePlot = 0;
p.writeToFile = 1;

p.peak = 1; % unit Hz
% p.sdThreshold = 2;
p.fieldTreshold = 0.2;
% p.timeGap = 1; % if field is seperated by 1 bin, concatenate time field
SecName = {'first','second','third'};

% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    
    if p.savePlot
        % directory for plot figures
        % generate a folder for each rat eah day under the current folder
        savedir = sprintf('%s%s%d%s%d%s',cd,'\Figures\',sessInfo(i).animal,'-day',sessInfo(i).day,'\DelayTimeField');
        if ~exist(savedir, 'dir')
            mkdir(savedir);
        end
    end
    
    savedir2 = sprintf('%s%s',sessInfo(i).mainDir,'\Cell Property');
    
    % load delayFire map
    delayFile = fullfile(savedir2, 'Fig8DelayTimeMap.mat');
    load(delayFile);
    % load correlation values
    stabilityFile = fullfile(sessInfo(i).mainDir,'Cell Property','DelayFireStability_TrialbyTrial_Shuffle.mat');
    load(stabilityFile);    
    stability30to10File = fullfile(sessInfo(i).mainDir,'Cell Property','DelayFireStability30to10_Shuffle.mat');
    load(stability30to10File);    
    
    % initiate the data
    DelayField_PartDelay.rat = sessInfo(i).animal;
    DelayField_PartDelay.day = sessInfo(i).day;
    DelayField_PartDelay.timeBin = DelayFire.timeBin;
    DelayField_PartDelay.gaussSigma = DelayFire.gaussSigma;

    TList = DelayFire.tList;
    clusterNum = length(DelayFire.tList);
    DelayField_PartDelay.tList = DelayFire.tList;
            
    % get processed data from each subfolder
    mainDir = sessInfo(i).mainDir;
    
    % get each phase names (no delay etc)
    SessDirs = sessInfo(i).sessDirs;    

    for j = 1:length(SessDirs)
        
        if contains(SessDirs{j},'30')
            for section = 1:3
                
                % def1, delay start at barrier
                shuffleSig_1 = DelayFireStability30to10_Shuffle.(SessDirs{j}).(SecName{section}).delayCorr_Def1_Trial_Shuffle_Label95;
                timeField_1 = nan(1,clusterNum);
                fieldLabel_1 = cell(1,clusterNum);
                % def2, delay start at entrance
                shuffleSig_2 = DelayFireStability30to10_Shuffle.(SessDirs{j}).(SecName{section}).delayCorr_Def2_Trial_Shuffle_Label95;
                timeField_2 = nan(1,clusterNum);
                fieldLabel_2 = cell(1,clusterNum);
        
                for k = 1:clusterNum

                    spikeRate1_Smooth = DelayFire.(SessDirs{j}).spikeRate1_Smooth{k};
                    spikeRate2_Smooth = DelayFire.(SessDirs{j}).spikeRate2_Smooth{k};                    
                    spikeRate1_Combined_Smooth = DelayFire.(SessDirs{j}).spikeRate1_Combined_Smooth{k};                    
                    spikeRate2_Combined_Smooth = DelayFire.(SessDirs{j}).spikeRate2_Combined_Smooth{k};
                    
                    binNum = floor(size(spikeRate2_Smooth,2)/3);   
                    
                    spikeRate1_Smooth_Temp = spikeRate1_Smooth(:,(section-1)*binNum+1:section*binNum);
                    spikeRate1_Combined_Smooth_Temp = spikeRate1_Combined_Smooth(:,(section-1)*binNum+1:section*binNum);               
                    peakRate_1 = max(spikeRate1_Combined_Smooth_Temp);
                    
                    spikeRate2_Smooth_Temp = spikeRate2_Smooth(:,(section-1)*binNum+1:section*binNum);
                    spikeRate2_Combined_Smooth_Temp = spikeRate2_Combined_Smooth(:,(section-1)*binNum+1:section*binNum);               
                    peakRate_2 = max(spikeRate2_Combined_Smooth_Temp);
            
                    % def 1
                    if shuffleSig_1(k) == 1
                        % calculate information/spike, info/sec, coherence, sparseness
                        [nFields_1,FieldBinX_1,fieldLabelTemp_1] = timeFieldSearch(spikeRate1_Combined_Smooth_Temp,p);
                    else
                        nFields_1 = 0;
                        fieldLabelTemp_1 = zeros(length(spikeRate1_Combined_Smooth_Temp),1);
                    end
                    timeField_1(k) = nFields_1;
                    fieldLabel_1{k} = fieldLabelTemp_1;
                    % def 2
                    if shuffleSig_2(k) == 1
                        % calculate information/spike, info/sec, coherence, sparseness
                        [nFields_2,FieldBinX_2,fieldLabelTemp_2] = timeFieldSearch(spikeRate2_Combined_Smooth_Temp,p);
                    else
                        nFields_2 = 0;
                        fieldLabelTemp_2 = zeros(length(spikeRate2_Combined_Smooth_Temp),1);
                    end
                    timeField_2(k) = nFields_2;
                    fieldLabel_2{k} = fieldLabelTemp_2;
            
            
            
%                     if j == 1
%                         h = figure(k);
%                         h.Position = [100 100 1600 800];
%                     else
%                         figure(k)
%                     end
%             
%                     subplot(2,4,j)
%                     imagesc([0,size(spikeRate2_Smooth_Temp,2)],[9 8+size(spikeRate2_Smooth_Temp,1)],spikeRate2_Smooth_Temp/max(max(spikeRate2_Smooth_Temp)))
%                     hold on
%                     imagesc([0,size(spikeRate2_Combined_Smooth_Temp,2)],7,spikeRate2_Combined_Smooth_Temp/max(spikeRate2_Combined_Smooth_Temp))
%                     set(gca,'YDir','normal')
%                     colormap(jet)
%                     if nFields ==1
%                         plot(FieldBinX-1.5,ones(1,length(FieldBinX))*6.3,'r','LineWidth',3)
%                     end
%                     
%                     if j == 1
%                         TITLE1 = sprintf('%s%d%s%d%s%s%s%s','DelayZone2-',sessInfo(i).animal,' Day-',sessInfo(i).day,'-',TList{k},'-',SessDirs{j});
%                     else
%                         if  contains(SessDirs{j},'30')
%                             TITLE1 = strcat(SessDirs{j},' First 10 sec');
%                         else
%                             TITLE1 = SessDirs{j};
%                         end
%                         
%                     end
%                     title({TITLE1},'Interpreter','None')
%                     axis off
%                     text2 = sprintf('%s%d','Shuffled significance: ',shuffleSig(k));
%                     text1 = sprintf('%s%2.2f%s','Peak rate : ',peakRate,' Hz');
%                     text3 = sprintf('%s%d','TimeField: ',nFields);
%                     text([1,1,1],[5,4,3],{text1,text2,text3},'FontSize',10,'Interpreter','None')
%                     ylim([0 20])
                end
                DelayField_PartDelay.(SessDirs{j}).(SecName{section}).timeField_1 = timeField_1;
                DelayField_PartDelay.(SessDirs{j}).(SecName{section}).fieldLabel_1 = fieldLabel_1;
                DelayField_PartDelay.(SessDirs{j}).(SecName{section}).timeField_2 = timeField_2;
                DelayField_PartDelay.(SessDirs{j}).(SecName{section}).fieldLabel_2 = fieldLabel_2;
        
            end
        end
    end
    
    if p.savePlot == 1
        for k = 1:clusterNum
            figure(k)
            figName = sprintf('%s%s%d%s%d%s%s%s',savedir,'\Rat-',sessInfo(i).animal,'-Day',sessInfo(i).day,'-',TList{k},'-DelayTimeField_Shuffle_PartDelay');
            print(figName,'-dpng','-r300');
        end
    end
    close all
    if p.writeToFile == 1
        save(fullfile(savedir2,'DelayTimeField_PartDelay.mat'), 'DelayField_PartDelay');
    end
    clear DelayField_Shuffle_PartDelay
    fprintf('Finished analysis for session %d\n',i);
end

end
