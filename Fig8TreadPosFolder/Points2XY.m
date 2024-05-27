function [x, y] = Points2XY(Points)

ij = cellfun(@(p) myfind(p), column2cell(Points));
[x, y] = cellfun(@(p,i) convert(p(1:i)), column2cell(Points), num2cell(ij));

end

function [x, y] = vt_bit2xy(bitfield_value)
	x = bitshift(bitshift(bitfield_value, 20, 'uint32'), -20, 'uint32');
	y = bitshift(bitshift(bitfield_value, 4, 'uint32'), -20, 'uint32');
end

function [x, y] = convert(x)
[x, y] = vt_bit2xy(x);
x = mean(x);
y = mean(y);
end

function C = column2cell(S)
C = mat2cell(S, size(S, 1), ones(1, size(S, 2)));
end

function i = myfind(p)
	i = find(p, 1, 'last');
	if isempty(i)
		i = 1;
	end
end