function out=checkmotseq(seqs);

  %motive similarity
  %seqs contains all sequences in a session
  %
  %compute matrix of  spearmans rs for  all combinations of sequence pairs
  out=[];
  
  for n=1:length(seqs)
    
    s1=seqs(n).id;
    
    out(n,n).rval=0;
    out(n,n).l=0;
    
    for m=n+1:length(seqs)
      s2=seqs(m).id;
      [rc lenc]=rankseq(s1,s2);
      
      
      
     
      out(n,m).rval=rc;
      out(n,m).l=lenc;
      
      out(m,n).rval=rc;
      out(m,n).l=lenc;
      
      
    end
    
  end
  
  
  