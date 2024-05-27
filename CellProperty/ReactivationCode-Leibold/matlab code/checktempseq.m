function [rval,len]=checktempseq(seqs, tmpl);

  %example: seqs(n).id=[8 3 4 12 1 10 11 19], n=1:N
  
  %correlate sequence with templates
  rval=nan(size(seqs));
  len=nan(size(seqs));
   
  %template seqeunce
  s2=tmpl;
    
  for n=1:length(seqs)
    
    %sequence from pop burst
    s1=seqs{n};     

    [rc lenc] = rankseq(s1,s2);  
    
    rval(n)=rc;
    len(n)=lenc;
    
  end
  