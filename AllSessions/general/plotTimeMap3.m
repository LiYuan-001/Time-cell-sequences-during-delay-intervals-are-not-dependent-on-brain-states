% sortMethod = 1, sort by map its own
% sort method = 2, map2 sort by map1
% no normalization
function plotTimeMap3(map1,map2,sortMethod,pos1,pos2,TITLEText1,TICKText1,TITLEText2,TICKText2)
if sortMethod == 1
    [~,peakIdx] = max(map1,[],2);
    [~,peak_Sort] = sort(peakIdx);
    cellMapTemp1_Sort = map1(peak_Sort,:);
    [~,peakIdx] = max(map2,[],2);
    [~,peak_Sort] = sort(peakIdx);
    cellMapTemp2_Sort = map2(peak_Sort,:);
else
    [~,peakIdx] = max(map1,[],2);
    [~,peak_Sort] = sort(peakIdx);
    cellMapTemp1_Sort = map1(peak_Sort,:);
    cellMapTemp2_Sort = map2(peak_Sort,:);
end

% cellMapTemp1_Sort_Norm = cellMapTemp1_Sort;
% cellMapTemp2_Sort_Norm = cellMapTemp2_Sort;

% cellMapTemp1_Sort_Norm = cellMapTemp1_Sort./max(cellMapTemp1_Sort,[],2);
% cellMapTemp2_Sort_Norm = cellMapTemp2_Sort./max(cellMapTemp2_Sort,[],2);

subplot(pos1(1),pos1(2),pos1(3))
imagesc(cellMapTemp1_Sort)
colormap(jet)
title(TITLEText1,'Interpreter','None')
axis on
set(gca, 'xtick', [0 size(cellMapTemp1_Sort,2)]);
set(gca, 'xticklabels', TICKText1);

subplot(pos2(1),pos2(2),pos2(3))
imagesc(cellMapTemp2_Sort)
colormap(jet)
title(TITLEText2,'Interpreter','None')
axis on
set(gca, 'xtick', [0 size(cellMapTemp2_Sort,2)]);
set(gca, 'xticklabels', TICKText2);

end