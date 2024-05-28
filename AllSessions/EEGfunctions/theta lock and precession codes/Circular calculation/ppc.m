% Pairwise phase consistency (PPC)
% NeuroImage 2010. The pairwise phase consistency_A bias-free measure of
% rhythmic neuronal synchronization
% Input   alpha, spike phases (rad)
% Output  p, Pairwise phase consistency (PPC)
function p = ppc(alpha)

% number of spikes
N = length(alpha);

% calculate sum of dot product
s=0;
for i=1:N-1
    for j=i+1:N
        s = s + cos(alpha(i))*cos(alpha(j))+sin(alpha(i))*sin(alpha(j));
    end
end

% avarege
p = 2/(N*(N-1))*s;
end
