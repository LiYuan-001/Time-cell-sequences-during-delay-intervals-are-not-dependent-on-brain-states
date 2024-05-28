function [p U UC] = circ_raotest(alpha)

% [p U UC] = circ_raotest(alpha)
%   Calculates Rao's spacing test by comparing distances between points on
%   a circle to those expected from a uniform distribution.
%
%   H0: Data is distributed uniformly around the circle. 
%   H1: Data is not uniformly distributed around the circle.
%
%   Alternative to the Rayleigh test and the Omnibus test. Less powerful
%   than the Rayleigh test when the distribution is unimodal on a global
%   scale but uniform locally.
%
%   Due to the complexity of the distributioin of the test statistic, we
%   resort to the tables published by 
%       Russell, Gerald S. and Levitin, Daniel J.(1995)
%       'An expanded table of probability values for rao's spacing test'
%       Communications in Statistics - Simulation and Computation
%   Therefore the reported p-value is the smallest alpha level at which the
%   test would still be significant. If the test is not significant at the
%   alpha=0.1 level, we return the critical value for alpha = 0.05 and p =
%   0.5.
%
%   Input:
%     alpha     sample of angles
%
%   Output:
%     p         smallest p-value at which test would be significant
%     U         computed value of the test-statistic u
%     UC        critical value of the test statistic at sig-level
%
%
%   References:
%     Batschelet, 1981, Sec 4.6
%
% Circular Statistics Toolbox for Matlab

% By Philipp Berens, 2009
% berens@tuebingen.mpg.de

alpha = alpha(:);

% for the purpose of the test, convert to angles
alpha = circ_rad2ang(alpha);
n = length(alpha);
alpha = sort(alpha);

% compute test statistic
U = 0;
lambda = 360/n;
for j = 1:n-1
    ti = alpha(j+1) - alpha(j);
    U = U + abs(ti - lambda);
end

tn = (360 - alpha(n) + alpha(1));
U = U + abs(tn-lambda);

U = (1/2)*U;

% get critical value from table
[p UC] = getVal(n,U);
end
