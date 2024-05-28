function plotTimeMap(map,pos,TITLEText,TICKText)
subplot(pos(1),pos(2),pos(3))
imagesc(map)
colormap(jet)
title(TITLEText,'Interpreter','None')
axis on
set(gca, 'xtick', [0 size(map,2)]);
set(gca, 'xticklabels', TICKText);
end