% Spike_Session is a strcut contains all sub sesssion as fields, and spike
% shapes are a 3D matrix [32samples,4ch,spikeNumber] written in each
% subfield
% Spike_Whole is a 3D matrix [32samples,4ch,spikeNumber] which contains all
% spikes from same cluster across all sub sessions
function [Spike_Session,Spike_Whole,Fs_Header,ADBitVolts] = loadSpike(sessInfo,tSp_Session,TList)

FieldSelectionFlags = [1,0,1,0,1];
HeaderExtractionFlag = 1;
ExtractMode = 5;

subDirs = fieldnames(tSp_Session);
SpikeShape = cell(length(TList),length(subDirs));
Spike_Whole = cell(length(TList),1);
ADBitVolts = zeros(12,1);

for j = 1:length(subDirs)
    for k = 1:length(TList)
        ExtractionModeVector = tSp_Session.(subDirs{j}){k}*10^6;
        
        TTName = strsplit(TList{k},'_');
        TTName = TTName{1};
        FileName = fullfile(sessInfo.mainDir,subDirs{j},strcat(TTName,'.ntt'));
        
        if ~length(ExtractionModeVector) == 0            
            [Timestamps, CellNumbers,Samples,Header] = ...
                Nlx2MatSpike(FileName, FieldSelectionFlags,...
                HeaderExtractionFlag, ExtractMode, ExtractionModeVector);
            
            if length(ExtractionModeVector) ~= length(Timestamps)
                error('Spike extrac number does not match')
            end
        else
            Header = [];
            Samples = [];
        end
        SpikeShape{k,j} = Samples;
        if size(Samples,1) > 0
            Spike_Whole{k} = cat(3,Spike_Whole{k},Samples);
        end
        
        TetID = str2num(TTName(3:end));
        % Get Frequency from Header
        if ~isempty(Header) && sum(contains(Header,'ADBitVolts'))>0
            % find the cell contains ADBitVolts and extract value
            a=contains(Header,'ADBitVolts');
            ADB = sscanf(Header{a},'%*s %f');
            ADBitVolts(TetID) = ADB(1);
        else
            warning('ADVolt not reported, spike input range default to 250');
            ADBitVolts(TetID) = 1/32767*250*10^-6;
        end

    end
end

% Get Frequency from Header
if ~isempty(Header) && sum(contains(Header,'SamplingFrequency'))>0
    % find the cell contains ADBitVolts and extract value
    a=contains(Header,'SamplingFrequency');
    Fs_Header = sscanf(Header{a},'%*s %d');
else
    warning('ADVolt not reported, sampling frequency default to 32000 Hz');
    Fs_Header = 32000;
end
    
SpikeShape2 = mat2cell(SpikeShape, size(SpikeShape, 1), ones(1, size(SpikeShape, 2)));
SpikeShape2 = mat2cell(SpikeShape2, size(SpikeShape2, 1), ones(1, size(SpikeShape2, 2)));

args = cell(length(SpikeShape2)*2, 1);
args(1:2:end) = subDirs;
args(2:2:end) = SpikeShape2;
Spike_Session = struct(args{:});
end