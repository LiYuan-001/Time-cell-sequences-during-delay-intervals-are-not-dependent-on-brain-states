function [rc,len]=rankseq(s1,s2)

  %compute rank order correlation between sequences
    
    
  % difference matrix
  d=ones(size(s1'))*s2 - s1'*ones(size(s2));
  
  % binary identity matrix
  d=(d==0);
  
  
  s=s1;
  s0=s2;
  if length(s1)<length(s2)
    s=s2;
    s0=s1;
    d=d';
  end
  % s0 is the shorter sequence
  
  
  
  %compute cell overlap (neurons contained in both)
  minseq=s(find(sum(d,2)>0));
  
  
  % delete neurons from the shorter sequence that are not in the minimal
  % sequence
  % 
  d0=ones(size(s0'))*minseq - s0'*ones(size(minseq));
  d0=(d0==0);
  s0=s0(sum(d0,2)>0);
  
  
  
  %find ordinal rank in the shorter sequence
  dd=minseq'*ones(size(s0))-ones(size(minseq'))*s0;
  [dum id]=min(abs(dd));
  len = length(id);
  %compute spearmans r
  if length(id)>1
  
    [r]=corrcoef(1:length(id),id);
    rc=r(1,2);

  else
    
    rc=nan;
  end
  
      
  