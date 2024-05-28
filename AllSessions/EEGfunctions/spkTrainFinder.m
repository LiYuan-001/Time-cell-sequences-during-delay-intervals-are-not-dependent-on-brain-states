function [nTrains,TrainProp] = spkTrainFinder(ts,maxisi,minSpk)

% counter for the number of bursts
nTrains = 0;
nTrainsTemp = 0;
TrainProp.Trainstart =[];
TrainProp.Trainend = [];
TrainProp.TrainSpks= zeros(1,length(ts));
    
isi = diff(ts);
n = length(ts);

trainLabel = nan(length(ts),1);
if numel(isi)>0  
    if isi(1) <= maxisi
       trainLabel(1) = 1;
       nTrainsTemp = nTrainsTemp+1;
    else
       trainLabel(1) = 0;
    end
    
    for t = 2:n-1
    % start of the burst label as 1
       if (isi(t-1)>maxisi) && (isi(t)<=maxisi)
          trainLabel(t) = 1;
          nTrainsTemp = nTrainsTemp+1;
%      single spike label as 0
       elseif (isi(t-1)>maxisi) && (isi(t)>maxisi)
          trainLabel(t) = 0;
     % spikes inside bursts label as 2  
       elseif (isi(t-1)<=maxisi) && (isi(t)<=maxisi)
           trainLabel(t) = 2;
     % spikes at the end of the burst label as 3
       elseif (isi(t-1)<=maxisi) && (isi(t)> maxisi)
           trainLabel(t) = 3;
       end        
    end   
    
    if (isi(n-1)<=maxisi) 
        trainLabel(n) = 3;   
    else
        trainLabel(n) = 0; 
    end
else
    nTrainsTemp = 0;
    trainLabel(1)=0;
end

if sum((trainLabel==1))+sum((trainLabel==2))+sum((trainLabel==3))...
        +sum((trainLabel==0)) ~= length(ts)
    error('Burst Spike assignment is wrong')
end

Trainstartidx = find(trainLabel==1);
Trainendidx = find(trainLabel ==3);

% limit spike train by minium spike number criteria
if nTrainsTemp > 0     
    for i = 1:nTrainsTemp
        if Trainendidx(i)-Trainstartidx(i)+1 >= minSpk
            nTrains = nTrains+1;
            TrainProp.Trainstart(nTrains) = Trainstartidx(i);
            TrainProp.Trainend(nTrains) = Trainendidx(i);
            TrainProp.TrainSpks(Trainstartidx(i):Trainendidx(i)) = 1;
        end
    end
end

end