function theta = thetaExtraction_Ratio(theta_Amp,theta_Phase,theta_Ts,tdr,tdr_timeStamps,p)

validIdx = [];
for i = 2:length(tdr)
    if tdr(i)>=p.tdrRatio
        [startTimeDifference,startIdx] = min(abs(theta_Ts-tdr_timeStamps(i-1)));
        [endTimeDifference,endIdx] = min(abs(theta_Ts-tdr_timeStamps(i)));
        if any([startTimeDifference,endTimeDifference] > tdr_timeStamps(2)-tdr_timeStamps(1))
                error('EEG could not find matching position')
        end
        validIdx = [validIdx,startIdx:endIdx];
    end
end
validIdx = unique(validIdx);
theta.LFPidx = validIdx;
theta.Amp = theta_Amp(validIdx);
theta.TimeStamp = theta_Ts(validIdx);
theta.phase = theta_Phase(validIdx);
end
