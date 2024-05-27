function plotpathspike(pdata, spkPosInd)

plot(pdata.x, pdata.y, 'color', ones(1,3)*.5);
hold on
scatter(pdata.x(spkPosInd), pdata.y(spkPosInd), 20, 'r', 'filled');