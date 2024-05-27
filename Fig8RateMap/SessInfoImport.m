function sessInfo = SessInfoImport(inFile)
% Define input parameters
% number of rows the program reads, increase if more then 700 rows are used up
max_number_of_i = 5000; 
worksheet = 'Sheet1';

% Read first column, to define all i numbers that occur
read_i = sprintf('A5:A%d',max_number_of_i);
[~, ~, SessIdx_Raw] = xlsread(inFile, worksheet, read_i);
SessIdx = reshape([SessIdx_Raw{:,1}],size(SessIdx_Raw));
SessIdx(isnan(SessIdx)) = [];

% Now import the rest of the data
read_total = sprintf('A5:T%d',max_number_of_i);
[~, ~, raw] = xlsread(inFile, worksheet ,read_total);
raw(cellfun(@(x) ~isempty(x) && isnumeric(x) && isnan(x),raw)) = {''};

% Transform the raw inpit into the right format
cellayerCH = cell(1,max(SessIdx));
diode = cell(1,max(SessIdx));

for row = 1:length(SessIdx_Raw)
    if ~isnan(SessIdx_Raw{row})  %only pick rows with valid i numbers
        i = SessIdx_Raw{row};
        
        if ischar(raw{row,2}); mainDir{i} = raw{row,2};
            oldFolder = cd;
            try cd(mainDir{i})  % check if paths are working
            catch, errordlg(sprintf('Path incorrect! Cannot read i %d, path: %s', i, mainDir{i})); end
            cd(oldFolder);
        end
        if isnumeric(raw{row,3});animal{i} = raw{row,3};end
        room(i) = raw(row,4);
        day(i) = raw(row,5);
        if isnumeric(raw{row,6});diode(i) = raw(row,6);end
        phaseSess = strrep(raw{row,7},' ','');
        phaseSess = strrep(phaseSess,'''','');
        bSess{i} = strsplit(phaseSess,',');   
        
        sleepSess = strrep(raw{row,8},' ','');
        sleepSess = strrep(sleepSess,'''','');
        sSess{i} = strsplit(sleepSess,','); 
        
        openSess = strrep(raw{row,9},' ','');
        openSess = strrep(openSess,'''','');
        oSess{i} = strsplit(openSess,',');
        
        TTList{i} = raw{row,10};
        PyrtList(i) = raw(row,11);
        IntList(i) = raw(row,12);
        region(i) = raw(row,13);
        cellayerCH(i) = raw(row,14);
        group(i) = raw(row,15);
        EEGch(i) = raw(row,16);
    end
end

fld = {
    'session'
    'mainDir'
    'animal'
    'room'
    'day'
    'diode'
    'sessDirs'
    'sleepDirs'
    'openDirs'
    'tList'
    'Pyr'
    'Int'
    'region'
    'cellLayerChann'
    'group'
    'EEGch'
    };

vars = {
    num2cell(1:max(SessIdx))
    mainDir
    animal
    room
    day
    diode
    bSess
    sSess
    oSess
    TTList
    PyrtList
    IntList
    region
    cellayerCH
    group
    EEGch
    };

args = cell(length(fld)*2, 1);
args(1:2:end) = fld;
args(2:2:end) = vars;

sessInfo = struct(args{:});
sessInfo = sessInfo(:);

end