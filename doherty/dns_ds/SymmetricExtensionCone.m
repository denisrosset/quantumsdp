function [Cons MainCons PPTCons] = SymmetricExtensionCone(coeffs, k, useSym, usePPT)
% SymmetricExtensionCone Compute SDP constraints corresponding to symmetric extension cones
%
% INPUTS
% coeffs     (dA^2)x(dB^2) matrix of type double/sdpvar containing the coefficients
%            of the Hermitian operator to constrain in the symmetric extension cone
%            H = sum_ij coeffs(i,j) * kron(FA{i}, FB{j})
%            where FA, FB are given by GeneralizedGellMann
% k          Number of copies
% useSym     Whether to use the symmetric subspace
% usePPT     Whether to use the PPT constraints
    dA = sqrt(size(coeffs, 1));
    dB = sqrt(size(coeffs, 2));
    [FA DA indPTA] = GeneralizedGellMann(dA);
    [FB DB indPTB] = GeneralizedGellMann(dB);
    factor = DB(1)^(k-1); % not used (because we define a cone), but the symmetric extension has
                          % larger trace by this factor (because the basis is not normalized)
    [~, B01] = BasisSymmetricSubspace(dB, k);
    DupB = B01'*B01;
    
    
    % indices of the rows/columns to preserve, because the
    % matrix we consider has support and range in the symmetric
    % subspace
    %
    % we consider the same for the PPT matrices
    symIndices = [];
    symIndicesPPT = cell(1, k);
    fieldDim = 2; % 1 if complex, 2 if real
    if useSym
        indicesB = SymmetricCanonicalIndices(dB, k);
        symIndices = indicesB(:)';
        dBext = length(symIndices);
        dBPPT = zeros(1, k);
        if usePPT
            % order is B1 ... Bk1 Bk1+1 ... Bk1+k2
            % but warning: kron reverses the order compared to ind2sub
            for k1 = 1:k
                k2 = k - k1;
                % thus indices1 has factor
                indices1 = (SymmetricCanonicalIndices(dB, k1) - 1) * dB^k2;
                indices2 = SymmetricCanonicalIndices(dB, k2);
                indicesB = bsxfun(@plus, indices1', indices2);
                symIndicesPPT{k1} = indicesB(:)';
                dBPPT(k1) = length(symIndicesPPT{k1});
            end
        end
    else
        dBext = dB^k;
        dBPPT = ones(1, k) * dB^k;
    end
    % order is A B1 ... Bk
    S = sparse(dA*dBext*fieldDim, dA*dBext*fieldDim);
    for j = 1:k
        PPT{j} = sparse(dA*dBPPT(j)*fieldDim, dA*dBPPT(j)*fieldDim);
    end
    binds = SymmetricCanonicalIndices(dB^2, k);
    for bind = binds
        bs = cell(1, k);
        [bs{:}] = ind2sub(dB^2*ones(1, k), bind);
        bs = cell2mat(bs);
        switch sum(bs > 1)
          case 0 % Alice marginal
            b = 1;
          case 1
            b = bs(find(bs > 1));
          otherwise
            b = 0;
        end
        allbs = unique(perms(bs), 'rows');
        B = sparse(dB^k, dB^k);
        for r = 1:size(allbs, 1)
            El = 1;
            for c = 1:k
                El = kron(El, FB{allbs(r, c)});
            end
            B = B + El;
        end
        for a = 1:dA^2
            if b == 0
                coeff = sdpvar;
            else
                coeff = coeffs(a, b);
            end
            A = FA{a};
            if useSym
                S = S + coeff * ComplexToReal(kron(A, B(symIndices, symIndices)));
            else
                S = S + coeff * ComplexToReal(kron(A, B));
            end
            if usePPT
                PT = full(B); % partial transpositions are added one by one
                for j = 1:k
                    d1 = dB^(j-1);
                    d2 = dB;
                    d3 = dB^(k-j);
                    PTsplit = reshape(PT, [d1 d2 d3 d1 d2 d3]);
                    PTsplit = permute(PTsplit, [1 5 3 4 2 6]);
                    PT = reshape(PTsplit, [dB^k dB^k]);
                    if useSym
                        PPT{j} = PPT{j} + coeff * ...
                                 ComplexToReal(kron(A, PT(symIndicesPPT{j}, symIndicesPPT{j})));
                    else
                        PPT{j} = PPT{j} + coeff * ComplexToReal(kron(A, PT));
                    end
                end
            end
        end
    end
    MainCons = [S >= 0];
    % PPT constraints
    PPTCons = [];
    if usePPT
        for k1 = 1:k
            PPTCons = [PPTCons; PPT{k1} >= 0];
        end
    end
    Cons = [MainCons
            PPTCons];
end
