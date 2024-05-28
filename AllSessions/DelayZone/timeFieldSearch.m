% this is to detect single field
function [nFields,FieldBinX,fieldLabel] = timeFieldSearch(map,p)

% Counter for the number of fields
FieldBinX = [];

% Allocate memory to the arrays
M = length(map);
% Array that contain the bins of the map this algorithm has visited
visited = zeros(M,1);
fieldLabel = visited;

nanInd = isnan(map);
visited(nanInd) = 1;

% Array that will contain the bin positions to the current placefield

% Find the current maximum
[peak,maxId] = max(map);
visited(map<p.fieldTreshold*peak) = 1;
% sd = std(map,'omitnan');

% if peak >= p.peak && peak>= nanmean(map) + p.sdThreshold*sd
if peak >= p.peak
    % look for time field
    
    % Find the bins that construct the peak field
    [binsX,visited] = recursiveBins_Lin(map,visited,[],maxId,M);
    binsX = sort(binsX);
    if length(binsX)>= 1 % Minimum size of a placefield
        nFields = 1;
        % Total rate
        visited(binsX) = 1; % modified by Li
        % Bins in Field
        FieldBinX = binsX;
        fieldLabel(binsX) = 1;
    end
else
    nFields = 0;
end
end


function [binsX,visited] = recursiveBins_Lin(map,visited,binsX,ii,M)
% If outside boundaries of map -> return.
if ii<1 || ii>M
    return;
end
% If all bins are visited -> return.
if prod(prod(visited))
    return;
end
if visited(ii) % This bin has been visited before
    return;
else
    binsX = [binsX;ii];
    visited(ii) = 1;
    % Call this function again in each of the 24 neighbour bins
    [binsX,visited] = recursiveBins_Lin(map,visited,binsX,ii+1,M);
    [binsX,visited] = recursiveBins_Lin(map,visited,binsX,ii-1,M);
end
end
