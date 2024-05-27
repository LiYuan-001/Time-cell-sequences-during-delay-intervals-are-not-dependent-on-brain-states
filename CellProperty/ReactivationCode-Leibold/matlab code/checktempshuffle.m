function [rvalR rvalL lenR lenL]=checktempshuffle(seqs, tmplR, tmplL, nrep);

  %rightward
  rvalR=[];
  lenR=[];
  
  %leftward
  rvalL=[];
  lenL=[];
  
  m=0;
  for nr=1:nrep
    
    
    for n=1:length(seqs)
      m=m+1;
      
      %shuffeld replay
      s1=seqs{n}(randperm(length(seqs{n})));
      
      
      %template Right
      s2=tmplR;
      
      [rc lenc] = rankseq(s1,s2);  
      
      rvalR(m)=rc;
      lenR(m)=lenc;
      
      %template left 
      s2=tmplL;
      
      [rc lenc] = rankseq(s1,s2);  
      
      rvalL(m)=rc;
      lenL(m)=lenc;
      
    end
  
  end