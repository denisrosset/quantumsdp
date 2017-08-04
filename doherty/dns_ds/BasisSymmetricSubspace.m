function [B B01] = BasisSymmetricSubspace(d, k)
% ProjectorSymmetricSubspace Returns a basis of the symmetric subspace
%
% INPUT
% d   Dimension of single copy
% k   Number of copies
%
% OUTPUT
% B   Orthonormal basis of the symmetric subspace
% B01 Orthogonal but not normalized basis with 0,1 coefficients
    B01 = sparse(0, d^k);
    for ind = 1:d^k
        subs = cell(1, k);
        [subs{:}] = ind2sub(d*ones(1, k), ind);
        subs = cell2mat(subs);
        if all(subs(2:end) - subs(1:end-1) >= 0)
            allsubs = unique(perms(subs), 'rows');
            element = sparse(1, d^k);
            for r = 1:size(allsubs, 1)
                subs = allsubs(r, :);
                subs = num2cell(subs);
                newind = sub2ind(d*ones(1, k), subs{:});
                element(newind) = 1;
            end
            B01(end + 1, :) = element;
        end
    end
    B = diag(1./sqrt(sum(B01, 2))) * B01;
end
