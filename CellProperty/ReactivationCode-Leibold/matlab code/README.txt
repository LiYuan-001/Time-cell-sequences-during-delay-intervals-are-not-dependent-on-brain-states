Readme for sequence statistics

Christian Leibold, 20.03.2018

includes MATLAB(R) scripts

testsession: masterfunction that uses the spearman rank correlations
computed via rankseq, checktempseq, checktempshuffle and rtfun to
compute cdfs for the observed seqeunces and the shuffeld nullhypothesis.


rankseq.m, checktempseq.m, checktempshuffle.m, rtfun.m, testsession.m, 
checkmotseq.m, rtfun2.m

rankseq: computes rank order correlation between two sequences. Input:
Seqeunces s1 and s2 are row vectors of neuron numbers. output is is
the spearmans rank order correlations and the minimal overlapping
sequence length.

checktempseq: calls rankseq for a set of sequences qand a place
field template sequence

checktempshuffle: calls rankseq for a random shuffled replay seqeunec
and the left and the rightward tempates to produce rvalues for null
hypothesis.

rtfun: surrogate set of rvalues for random shuffles of different
lengths (input argument is a vector of lengths). ich radum sequnecs is
correates to two fixed permutations reflecting the testing against
template sequences.


testmot: same as testsession just for motifs

checkmotseq: same as checktempseq just for motifs (no tempaltes but r
values within the set of sequences)

rtfun2: same as rtfun for only one tempaltes seqeunces (not two
directional templates)


