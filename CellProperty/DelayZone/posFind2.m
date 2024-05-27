% start delay at around delay entrance
function  delayPos = posFind2(pathData,locInfo,startidx,endidx)

for i = 1:length(startidx)
    
    posInd1 = find(pathData.t == locInfo.tInt(startidx(i),1));
    posInd2 = find(pathData.t == locInfo.tInt(endidx(i),2));
    
    if isempty(posInd1) || isempty(posInd2)
        error('Cannot find such position')
    end
    
    posData = pathData.y(posInd1:posInd2);
    nanIdx = isnan(posData);
    delayIdx = posData > 0;
    delayIdx(nanIdx) = 1;
    a = find(diff(delayIdx)==-1);
    delayIdx = a(1);
    
    delayPos.start(i) = posInd1;
    delayPos.end(i) = posInd1+delayIdx;
    delayPos.startT(i) = pathData.t(posInd1);
    delayPos.endT(i) = pathData.t(posInd1+delayIdx);
end
end 
