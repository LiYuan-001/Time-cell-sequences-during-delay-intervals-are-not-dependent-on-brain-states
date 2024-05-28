% Perform Sub sampling of spikes & calculate phase locking parameters
% 
% spkPhase, spike phase (rad)
% spkMin, Sub sample to spkMin
% r,     mean of resulatant length of Sub sampled phases
% mu,    mean of mu                of Sub sampled phases
% kappa, mean of kappa             of Sub sampled phases
function [r,mu,kappa] = subSmaple(spkPhase,spkMin)

% total spike number
N = length(spkPhase);

% 100 bootstrap (sub-sampling) samples without replacement
for i=1:100
    % Get random index for Sub sampling
    p = randperm(N);
    idx = p(1:spkMin);
    
    % Random Sub sampling
    D_spkPhase = spkPhase(idx);
    
    % Resultant vector length
    D_r(i) = circ_r(D_spkPhase);
    
    % Concentration parameter kappa
    [D_mu(i) D_kappa(i)] = circ_vmpar(D_spkPhase);
end

% Mean of Sub sampled parameters
r = mean(D_r);
mu = mean(D_mu);
kappa = mean(D_kappa);
end