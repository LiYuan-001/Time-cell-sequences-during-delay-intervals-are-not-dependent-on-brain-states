function [rval,l]=checkmotseq_LiModified(seqs)

  %motive similarity
  %seqs contains all sequences in a session
  %
  %compute matrix of  spearmans rs for  all combinations of sequence pairs
  rval = zeros(length(seqs),length(seqs));
  l = zeros(length(seqs),length(seqs));
  for n=1:length(seqs)-1
    
    s1=seqs{n};

    for m=n+1:size(seqs,1)
      s2=seqs{m};
      [rc lenc]=rankseq(s1,s2);
      
 
      rval(n,m)=rc;
      l(n,m)=lenc;
      
      rval(m,n)=rc;
      l(m,n)=lenc;
        
    end
    
  end
  
  
  