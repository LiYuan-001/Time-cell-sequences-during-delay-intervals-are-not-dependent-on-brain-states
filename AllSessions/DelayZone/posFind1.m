% start delay at around delay barrier
function  delayPos = posFind1(pathData,locInfo,startidx,endidx)

for i = 1:length(startidx)
    
    posInd1 = find(pathData.t == locInfo.tInt(startidx(i),1));
    posInd2 = find(pathData.t == locInfo.tInt(endidx(i),2));
    
    if isempty(posInd1) || isempty(posInd2)
        error('Cannot find such position')
    end
    posData = pathData.y(posInd1:posInd2);

    tempPosLabel = find(posData<=50);   
    start1Idx = tempPosLabel(1);
    
    posData = pathData.y(posInd1:posInd2);
    nanIdx = isnan(posData);
    
    delayIdx = posData > 10;
    delayIdx(nanIdx) = 1;
    a = find(diff(delayIdx)==-1);
    delayIdx = a(1);
    
    delayPos.start(i) = start1Idx+posInd1-1;
    delayPos.end(i) = delayIdx+posInd1-1;
    delayPos.startT(i) = pathData.t(start1Idx+posInd1-1);
    delayPos.endT(i) = pathData.t(delayIdx+posInd1-1);
end
end 
