function FigSpikeShapeDist(inFile,AnalyzeSes)

spikeShape.Cell = [];
spikeShape.peak2Tr = [];
spikeShape.Avgrate = [];
spikeShape.Deriv_2 = [];

% Read in input information
sessInfo = SessInfoImport(inFile);

for i = AnalyzeSes(1:end)
    SpikeShapeFile = fullfile(sessInfo(i).mainDir,'Cell Property','SpikeShape.mat');
    load(SpikeShapeFile);
    
    for k = 1:length(SpikeProp.tFileNames)
        cellName = sprintf('%s%d%s%d%s','Rat-',SpikeProp.rat,'Day-',SpikeProp.day,SpikeProp.tFileNames{k});
        cellName2{1} = cellName;
        spikeShape.Cell = [spikeShape.Cell,cellName2{1}];
    end
    spikeShape.peak2Tr = [spikeShape.peak2Tr,SpikeProp.peak2tr];
    spikeShape.Avgrate = [spikeShape.Avgrate,SpikeProp.max_AvgRate];
    spikeShape.Deriv_2 = [spikeShape.Deriv_2,SpikeProp.Deriv_2];
end

plot3(spikeShape.peak2Tr,spikeShape.Avgrate,spikeShape.Deriv_2,'o')
xlabel('peak2tr')
ylabel('AvgRate')
zlabel('Trough Deriv2')
end