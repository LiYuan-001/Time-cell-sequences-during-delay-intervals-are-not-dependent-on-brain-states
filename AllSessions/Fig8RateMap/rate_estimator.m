% Calculate the rate for one position value
function [r,edge_corrector] = rate_estimator(spkx,spky,x,y,invh,posx,posy,post)
% edge-corrected kernel density estimator

invhy = invh; % for linear maps

conv_sum = sum(gaussian_kernel(((spkx-x)*invh),((spky-y)*invh)));
edge_corrector =  trapz(post,gaussian_kernel(((posx-x)*invh),((posy-y)*invh)));
% edge_corrector(edge_corrector<0.01) = NaN;

r = (conv_sum / (edge_corrector + 0.1)) + 0.1; % regularised firing rate for "wellbehavedness"
% i.e. no division by zero or log of zero
% Gaussian kernel for the rate calculation
function r = gaussian_kernel(x,y)
% k(u) = ((2*pi)^(-length(u)/2)) * exp(u'*u)
r = 0.15915494309190 * exp(-0.5*(x.*x + y.*y));