% +------------------------------------------------------+
% |           Angular Distance and Similarity            |
% |              with MATLAB Implementation              | 
% |                                                      |
% | Author: Ph.D. Eng. Hristo Zhivomirov        07/31/21 | 
% +------------------------------------------------------+
% 
% function: [ad, as] = angmetric(x, y)
%
% Input:
% x - signal in the time domain; x could be vector or 
%     matrix with time across columns and indexes across rows
% y - signal in the time domain; y could be vector or 
%     matrix with time across columns and indexes across rows
% 
% Output:
% ad - angular distance between two vectors or matrices columns. Range: [0 1] 
% as - angular similarity between two vectors or matrices columns. Range: [0 1] 

function [ad, as] = angmetric(x, y)

% check if x is vector and if it 
% is - convet it to a column-vector
if isvector(x), x = x(:); end 

% check if y is vector and if it 
% is - convet it to a column-vector
if isvector(y), y = y(:); end 

% calculate the cosine similarity
xcell = num2cell(x, 1);
ycell = num2cell(y, 1);
xnorm = cellfun(@(x) norm(x), xcell);
ynorm = cellfun(@(x) norm(x), ycell);
cs = dot(x, y)./xnorm./ynorm;

% calculate the angular distance
% Note: the range of ad is [0 1]
ad = acos(cs)/pi;

% calculate the angular similarity
% Note: the range of as is [0 1]
as = 1 - ad;

end