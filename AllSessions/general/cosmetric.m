% +------------------------------------------------------+
% |            Cosine Distance and Similarity            |
% |              with MATLAB Implementation              | 
% |                                                      |
% | Author: Ph.D. Eng. Hristo Zhivomirov        07/31/21 | 
% +------------------------------------------------------+
% 
% function: [cd, cs] = cosmetric(x, y)
%
% Input:
% x - signal in the time domain; x could be vector or 
%     matrix with time across columns and indexes across rows
% y - signal in the time domain; y could be vector or 
%     matrix with time across columns and indexes across rows
% 
% Output:
% cd - cosine distance between two vectors or matrices columns. Range: [0 2] 
% cs - cosine similarity between two vectors or matrices columns. Range: [-1 1]

function [cd, cs] = cosmetric(x, y)

% check if x is vector and if it 
% is - convet it to a column-vector
if isvector(x), x = x(:); end 

% check if y is vector and if it 
% is - convet it to a column-vector
if isvector(y), y = y(:); end 

% calculate the cosine similarity
% Note: the range of cs is [-1 1]
xcell = num2cell(x, 1);
ycell = num2cell(y, 1);
xnorm = cellfun(@(x) norm(x), xcell);
ynorm = cellfun(@(x) norm(x), ycell);
cs = dot(x, y)./xnorm./ynorm;

% calculate the cosine distance
% Note: the range of cd is [0 2]
cd = 1 - cs;

end