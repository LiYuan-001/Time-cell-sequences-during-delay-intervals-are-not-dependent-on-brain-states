function[nBursts,BurstProp,nSinglespikes,SinglespikeProp] = burstfinder(ts,maxisi)

% counter for the number of bursts
nBursts = 0;
nSinglespikes = 0;

BurstProp =[];
SinglespikeProp =[];

% For testing Buzsaki's 2001 burst paper. Silence before burst

burstsilence =[]; % time interval of previous one spike to burst
singlespikesilence =[];

isi = diff(ts);
n = length(ts);

if numel(isi)>0  
    if isi(1) <= maxisi
       burstlabel(1) = 1;
       nBursts = nBursts+1;
       burstsilence = [burstsilence;0];
    else
       burstlabel(1) = 0;
       nSinglespikes = nSinglespikes+1;
       singlespikesilence = [singlespikesilence;0];
    end
    
    for t = 2:n-1
    % start of the burst label as 1
       if (isi(t-1)>maxisi) && (isi(t)<=maxisi)
          burstlabel(t) = 1;
          nBursts = nBursts+1;
          burstsilence = [burstsilence; isi(t-1)];
%      single spike label as 0
       elseif (isi(t-1)>maxisi) && (isi(t)>maxisi)
          burstlabel(t) = 0;
          nSinglespikes = nSinglespikes+1;
          singlespikesilence = [singlespikesilence;isi(t-1)];
     % spikes inside bursts label as 2  
       elseif (isi(t-1)<=maxisi) && (isi(t)<=maxisi)
           burstlabel(t) = 2;
     % spikes at the end of the burst label as 3
       elseif (isi(t-1)<=maxisi) && (isi(t)> maxisi)
           burstlabel(t) = 3;
       end        
    end   
    
    if (isi(n-1)<=maxisi) 
        burstlabel(n) = 3;   
    else
        burstlabel(n) = 0; 
        singlespikesilence = [singlespikesilence;isi(n-1)];
        nSinglespikes = nSinglespikes+1;
    end
else
    nBursts = 0;
    nSinglespikes = 1;
    burstlabel(1)=0;
end

if sum((burstlabel==1))+sum((burstlabel==2))+sum((burstlabel==3))...
        +sum((burstlabel==0)) ~= length(ts)
    error('Burst Spike assignment is wrong')
end

Burststartidx = find(burstlabel==1);
Burstendidx = find(burstlabel ==3);
Singlespikeidx = find(burstlabel ==0);

BurststartTs = ts(Burststartidx);
BurstendTs = ts(Burstendidx);
SinglespikeTs = ts(Singlespikeidx);

if nBursts == 0
    BurstProp.Burststart = [];
    BurstProp.Burstend = [];
    BurstProp.BurststartTs = [];
    BurstProp.BurstendTs = [];
    BurstProp.BurstPreSilence = [];
    BurstProp.BurstSpks= [];
    
else
    BurstProp.Burststart = Burststartidx;
    BurstProp.Burstend = Burstendidx;
    BurstProp.BurststartTs = BurststartTs;
    BurstProp.BurstendTs = BurstendTs;
    BurstProp.BurstPreSilence = burstsilence;
    BurstProp.BurstSpks= find(burstlabel ~=0);
end

if nSinglespikes == 0
    SinglespikeProp.Singlespikeidx = [];
    SinglespikeProp.SinglespikeTs = [];
    SinglespikeProp.SinglepreSilence = []; 
else    
    SinglespikeProp.Singlespikeidx = Singlespikeidx;
    SinglespikeProp.SinglespikeTs = SinglespikeTs;
    SinglespikeProp.SinglepreSilence = singlespikesilence;
end