function [Spikes_Session, tfiles] = loadSpikeTime_Treadmill(sessInfo,p)
SessDirs = sessInfo.sessDirs;
sleepDirs = sessInfo.sleepDirs;

subDirs = [SessDirs,sleepDirs];

% first to detect whether there is TTlist
if ~exist(fullfile(sessInfo.mainDir, sessInfo.tList), 'file')
    error('No TTlist exist for session: %s',sesInfo.mainDir);
end

% get cluster name from TTlist(name of .t file saved from mClust)
if ischar(sessInfo.tList)
    tfiles = ReadFileList2(fullfile(sessInfo.mainDir, sessInfo.tList));
elseif iscellstr(ttlist) % if ttlist itself is the tfiles
    tfiles = ttlist;
end

% read spikes from either subsession or whole session
if strcmpi(p.spikeMode,'subses')   
    % read in spikes with the unit of second
    % initiate cells to put spike time in
    tSp = cell(length(tfiles),length(subDirs));
    for ii = 1:length(subDirs)
        for jj = 1:length(tfiles)
            tFilePath = fullfile(sessInfo.mainDir,subDirs{ii},tfiles{jj}); % Full path to the .t file
            tSp{jj,ii} = readtfile2(tFilePath);
        end
    end
    
elseif strcmpi(extp.spike, 'wholeses')
    tSp = cell(1,length(tfiles));
    for jj = 1:length(tfiles)
            tFilePath = fullfile(sessInfo.mainDir,tfiles{jj}); % Full path to the .t file
            tSp{jj,1} = readtfile2(tFilePath);
    end
end

if ~ismatrix(tSp)
	error('tSp must be a 2-d array')
end

tSp = mat2cell(tSp, size(tSp, 1), ones(1, size(tSp, 2)));
tSp = mat2cell(tSp, size(tSp, 1), ones(1, size(tSp, 2)));

args = cell(length(tSp)*2, 1);
args(1:2:end) = subDirs;
args(2:2:end) = tSp;
Spikes_Session = struct(args{:});
end

function F = ReadFileList2(fn)

% F = ReadFileList(fn)
%
% INPUTS:
%   fn -- an ascii file of filenames, 1 filename per line
%
% OUTPUTS:
%   F -- a cell array of filenames suitable for use in programs
%        such as LoadSpikes
%
% Now can handle files with headers

% ADR 1998
% version L4.1
% status: PROMOTED

% v4.1 added ReadHeader

[fp,errmsg] = fopen(fn, 'rt');
if (fp == -1)
    error(['Could not open "', fn, '". ', errmsg]);
end

ifp = 1;
while (~feof(fp))
    F{ifp} = fgetl(fp);
    ifp = ifp+1;
end
fclose(fp);

F = F';
end

