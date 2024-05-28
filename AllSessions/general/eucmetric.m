% +------------------------------------------------------+
% |          Euclidean Distance and Similarity           |
% |              with MATLAB Implementation              | 
% |                                                      |
% | Author: Ph.D. Eng. Hristo Zhivomirov        07/31/21 | 
% +------------------------------------------------------+
% 
% function: [ed, es] = eucmetric(x, y)
%
% Input:
% x - signal in the time domain; x could be vector or 
%     matrix with time across columns and indexes across rows
% y - signal in the time domain; y could be vector or 
%     matrix with time across columns and indexes across rows
% 
% Output:
% ed - Euclidean distance between two vectors or matrices columns. Range: [0 Inf] 
% es - Euclidean similarity between two vectors or matrices columns. Range: [0 1]

function [ed, es] = eucmetric(x, y)

% check if x is vector and if it 
% is - convet it to a column-vector
if isvector(x), x = x(:); end 

% check if y is vector and if it 
% is - convet it to a column-vector
if isvector(y), y = y(:); end 

% calculate the Euclidean distance
% Note: the range of ed is [0 Inf]
xcell = num2cell(x, 1);
ycell = num2cell(y, 1);
ed = cellfun(@(x, y) norm(x-y), xcell, ycell);

% calculate the Euclidean similarity
% Note: the range of es is [0 1]. This is similarity based not on  
% the square of Euclidean distance (as is shown in the [1]) but on
% the true Euclidean distance, i.e. it not emphasizes the distances
% so boldly.
es = 1./(1 + ed);

end