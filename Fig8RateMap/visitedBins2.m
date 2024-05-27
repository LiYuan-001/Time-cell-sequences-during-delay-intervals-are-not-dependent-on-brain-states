% Calculates what area of the map that has been visited by the rat
function visited = visitedBins2(posx,posy,mapAxis1,mapAxis2)

binWidth = mapAxis1(2)-mapAxis1(1);

% Number of bins in each direction of the map
N = length(mapAxis1);
M = length(mapAxis2);
visited = zeros(M,N);

for ii = 1:N
    for jj = 1:M
        px = mapAxis1(ii);
        py = mapAxis2(jj);
        distance = sqrt( (px-posx).^2 + (py-posy).^2 );
        
        if min(distance) <= binWidth
            visited(jj,ii) = 1;
        end
    end
end
