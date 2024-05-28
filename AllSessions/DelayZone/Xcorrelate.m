% Xcorrelate calculates the pearson correlation between the rate map from two
% cells. Pixels that are visited less than 150 ms in the maps are also omitted
% from both cells.
function [corrValue,pix] = Xcorrelate(map1,map2)
% Both the rate maps and the time maps will have equal size
[M,N] = size(map1);
[M2,N2] = size(map2);

if M~=M2 || N~=N2
    error('Correlation sizes do not match!');
else
% % Exclude bins that was visited less than 150 ms in either room
% for ii = 1:N
%     for jj = 1:M
%         if timeMap1(ii,jj) < 0.150 | timeMap2(ii,jj) < 0.150
%             map1(ii,jj) = NaN;
%             map2(ii,jj) = NaN;
%         end
%     end
% end

% Transform the 2-D maps to 1-D arrays by assembling the columns from the
% maps
    A = reshape(map1,M*N,1);
    B = reshape(map2,M*N,1);
    % Find the pixels containing NaN in A
    An = isnan(A);
    index = find(An==0);
    % Remove the bins with NaN in A from both maps
    A = A(index);
    B = B(index);
    % Find the pixels that still contain NaN in B
    Bn = isnan(B);
    index = find(Bn==0);
    % Remove the bins with NaN in B from both maps
    A = A(index);
    B = B(index);
    % Number of bins used in the correlation
    pix = length(A);

    % Calculate the correlation for the two rate maps.
    if length(A) ~= length(B)
        error('Correlation sizes do not match!');
    elseif isempty(A) || isempty(B)
            corrValue = NaN; % Cannot do calculation on empty arrays
        else
            corrC = corrcoef(A,B,'Rows','pairwise');
            corrValue = corrC(1,2);
    end
end
end
