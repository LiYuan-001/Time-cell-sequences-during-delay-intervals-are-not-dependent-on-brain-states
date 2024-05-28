function [spkx,spky,newTsp,spkPosInd] = spk2posInd(tsp,posx,posy,post)

% spkpost is to indicate spike time ts in postion time post idx
spkPosInd = [];

N = length(tsp);
spkx = nan(N,1);
spky = nan(N,1);
newTsp = nan(N,1);
count = 0;
for ii = 1:N
    tdiff = (post-tsp(ii)).^2;
    [m,ind] = min(tdiff);
    if m < 1/25
        count = count + 1;
        spkx(count) = posx(ind(1));
        spky(count) = posy(ind(1));
        newTsp(count) = tsp(ii);
        spkPosInd = [spkPosInd;ind(1)];
    end
end
spkx = spkx(1:count);
spky = spky(1:count);
newTsp = newTsp(1:count);
end