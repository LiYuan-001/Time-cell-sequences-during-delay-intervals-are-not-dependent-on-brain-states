function [x, y] = polycenter(x, y)
%POLYCENTER Find center point of enclosing square for a polygon
%
% [x, y] = POLYCENTER(x, y)

% Siavash Ahmadi
% 4/6/2017

x = mean([min(x), max(x)]);

if nargout > 1
	y = mean([min(y), max(y)]);
end