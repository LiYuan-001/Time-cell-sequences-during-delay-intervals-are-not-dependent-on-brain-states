% resample timestamp and EEG to make them match to each other
function EEGTs = timeStampResample(Samples,Timestamps,Fs)
    % ---------------------------------------------------------------------
    % resample timestamp to match EEG
    sampInsert = size(Samples,1)./size(Timestamps,1);
    
    % double check if Fs matchs sample time
    if (Timestamps(2)-Timestamps(1))/sampInsert ~= 10^6/Fs
        error('Sample Time is wrong')
    end
    
    insertTime = 0:10^6/Fs:(Timestamps(2)-Timestamps(1)-10^6/Fs);
    if size(insertTime,1)==1
        insertTime = insertTime';
    end
    
    insertTime2 = repmat(insertTime,1,size(Samples,2));
    TimeStamp = repmat(Timestamps,size(Samples,1),1);
    
    EEGTs = TimeStamp + insertTime2;
    EEGTs = reshape(EEGTs,length(Samples(:)),1);
    