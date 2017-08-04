function [F D indPT] = GeneralizedGellMann(d)
% GeneralizedGellMann Orthogonal generalized Gell-Mann matrices with integer coefficients
%
% INPUTS
% d      Dimension
%
% OUTPUTS
% F      Cell array of size (d^2, 1) containing a basis of the Hermitian matrices in
%        C^{d x d} such that trace(F_i F_j) = D(j)*(i == j); thus the matrices are orthogonal 
%        with respect to the Hilbert-Schmidt inner product, but not normalized.
%        F{1} = eye(d).
%
%
% D      Normalization factor
%
% indPT  Indices of the basis elements that take a minus sign under partial transposition 
    F = cell(d, d);
    F{1,1} = sparse(eye(d));
    indPT = zeros(d, d);
    B = OrthoIntegerBasis(d);
    for k = 1:d
        for j = 1:d
            if k < j
                F{k,j} = sparse([j k], [k j], [1 1], d, d);
            elseif k > j
                F{k,j} = sparse([j k], [k j], [-1i 1i], d, d);
                indPT(k,j) = 1;
            else % k == j
                F{k,j} = diag(B(:,k));
            end
        end
    end
    F = F(:);
    indPT = logical(indPT(:));
    for k = 1:d*d
        for j = 1:d*d
            D(k,j) = trace(F{k}*F{j});
        end
    end
    assert(isequal(D, diag(diag(D))));
    D = diag(D);
end
