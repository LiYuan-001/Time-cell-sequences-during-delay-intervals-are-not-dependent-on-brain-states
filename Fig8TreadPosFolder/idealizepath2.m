function [pathDataIdeal, locInfo, maze] = idealizepath2(pathData,mazeType,p)
if nargin<2
    mazeType = 'fig8rat';
end
[locInfo, maze] = mazeLocInfo(pathData,mazeType);
nPaths = length(pathData);
pathDataIdeal = pathData;
for path = 1:nPaths
    labelSeq = locInfo(path).labelSeq;
    locInds = locInfo(path).inds;
	locInfo(path).mazeType = mazeType;
    nSeq = length(labelSeq);
    x = pathData(path).x;
    y = pathData(path).y;
    t = pathData(path).t;
    for i = 1:nSeq
        loc = labelSeq{i};
        inds = locInds(i,1):locInds(i,2);
        if isempty(loc)
            x(inds) = nan;
            y(inds) = nan;
            continue;
        end
        if i == 1
            entryExitLoc = {labelSeq{i+1},labelSeq{i+1}};
        elseif i == nSeq
            entryExitLoc = {labelSeq{i-1},labelSeq{i-1}};
        else
            entryExitLoc = {labelSeq{i-1},labelSeq{i+1}};
        end
        xBnd = maze.locs.(loc)(:,1);
        yBnd = maze.locs.(loc)(:,2);
        switch loc
            case {'A16','A25','A34'}
                entryExit = {'t','b'};
            case {'A12','A23','A45','A56'}
                entryExit = {'l','r'};
            case 'N1'
                entryExit = {'b','r'};
            case 'N2'
                code = {'l','b','r'};
                adj = {'A12','A25','A23'};
                entryExit = code(ismember(adj,entryExitLoc));
            case 'N3'
                entryExit = {'l','b'};
            case 'N4'
                entryExit = {'t','l'};

            case 'N5'
                code = {'r','t','l'};
                adj = {'A45','A25','A56'};
                entryExit = code(ismember(adj,entryExitLoc));
            case 'N6'
                entryExit = {'r','t'};
        end
        [x(inds), y(inds)] = idealizesegment2(x(inds),y(inds),xBnd,yBnd,entryExit,loc);
    end
    
    [v vmean] = getVelocity2(x',y',t');
    speedIdx = (v > 100);
    while sum(speedIdx)>0
        [v vmean] = getVelocity2(x',y',t');
        speedIdx = (v > p.speed);
%         plot(x(speedIdx),y(speedIdx),'r.')
%         hold on
        x(speedIdx) = NaN;
        y(speedIdx) = NaN;
    end
%     [x,y] = interpolatePosition_sl(x,y,p.timeThreshold,p.sampRate);
    
    pathDataIdeal(path,1).x = x;
    pathDataIdeal(path,1).y = y;
    pathDataIdeal(path,1).v = v;
end