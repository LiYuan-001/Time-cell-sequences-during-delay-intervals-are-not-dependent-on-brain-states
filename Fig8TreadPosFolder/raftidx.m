function L = raftidx(logicalArray)

% L = repmat(find(lau.rton(logicalArray)), 2,1);
% if isempty(L)
% 	L = [];
% 	return
% end
% L(2,:) = L(2,:) + lau.raftsize(logicalArray) - 1;

% uncommented above line and wrote the following @date 12/06/2015 1:50 AM
% accomplishes the same objective
if isempty(logicalArray)
	L = logicalArray;
	return;
end

L = [find(rton(logicalArray(:)'==1)); find(rtoff(logicalArray(:)'==1))];
L = L(:);

if isempty(L)
	L = [];
end

L = reshape(L, 2, length(L)/2);