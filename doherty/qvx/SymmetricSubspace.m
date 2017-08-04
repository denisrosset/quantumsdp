classdef SymmetricSubspace
    properties(SetAccess = immutable)
        d;
        n;
        dim;
        cumProdFull;
        dimTable;
    end
    methods
        function obj = SymmetricSubspace(d, n)
            obj.d = d;
            obj.n = n;
            obj.dimTable = zeros(d, n);
            for i = 1:d
                for j = 1:n
                    obj.dimTable(i, j) = SymmetricSubspace.computeDimension(i, j);
                end
            end
            if n == 0
                obj.dim = 1;
            else
                obj.dim = obj.dimTable(d, n);
            end
            cp = cumprod(d*ones(1, n));
            obj.cumProdFull = [1 cp(1:end-1)];
        end
        function sub = indToSubSym(obj, ind)
        % indToSubSym Convert from an index of the symmetric basis to subindices
        %
        % Let dim be the dimension of this SymmetricSubspace(d, n)
        %
        % INPUT  ind is a vector of m indices of the symmetric basis, values in 1...dim
        % OUTPUT sub is a (m x n) matrix, the r-th row is defined with respect to ind(r)
        %        sub(r, :) = [c1 c2 ... cn] with 1 <= c1 <= <= c2 <= ... <= cn <= d
            ind = ind(:);
            nRows = length(ind);
            switch obj.n
              case 0
                sub = zeros(nRows, 0);
              case 1
                sub = ind;
              otherwise
                [sortedInd, index] = sort(ind);
                nRows = length(ind);
                sortedSub = SymmetricSubspace.indToSubSymHelper(sortedInd, obj.d, obj.n, obj.dimTable);
                sub = zeros(nRows, obj.n);
                sub(index, :) = sortedSub;
            end
        end
        function ind = subToIndSym(obj, sub)
        % subToIndSym Converts from subindices to the index in the symmetric basis
        %
        % Let dim be the dimension of this SymmetricSubspace(d, n)
        %
        % INPUT:  sub is a (m x n) matrix, whose r-th row is composed of subindices
        %         1 <= c1 ... cn <= d
        %
        % OUTPUT  ind is a vector of length m, whose r-th element correspond to the 
        %         index of sub(r, :) in the symmetric basis
            nRows = size(sub, 1);
            switch obj.n
              case 0
                ind = ones(nRows, 1);
              case 1
                ind = sub;
              otherwise
                sub = sort(sub')'; % put each row of subindices in increasing order
                                   % sort treats columns individually, so transpose
                [sortedSub, index] = sortrows(sub); % sorts the subindices for speed
                                                    % we have sortedSub = sub(index, :)
                rStart = 1;
                sortedInd = SymmetricSubspace.subToIndSymHelper(sortedSub, obj.d, obj.n, obj.dimTable);
                ind = zeros(nRows, 1);
                ind(index) = sortedInd;
            end
        end
        function sub = indToSubFull(obj, ind)
            ind = ind(:) - 1;
            sub = zeros(length(ind), obj.n);
            for i = 1:obj.n
                sub(:, i) = mod(ind, obj.d);
                ind = ind - sub(:, i);
                sub(:, i) = sub(:, i) + 1;
                ind = ind ./ obj.d;
            end
        end
        function ind = subToIndFull(obj, sub);
            assert(size(sub, 2) == obj.n);
            ind = (obj.cumProdFull * (sub' - 1) + 1)';
        end
    end
    methods(Static)        
        function sub = indToSubSymHelper(ind, d, n, dimTable)
        % helper function for indToSubSym
        % assumes that ind is already sorted
            if n == 1 % handle the trivial case
                sub = ind;
                return
            end
            nRows = length(ind);
            rStart = 1;
            nElementsBefore = 0;
            sub = zeros(nRows, n);
            % for speed, we select the block of indices whose first subindex value
            % is the same; then we call indToSubSymHelper recursively for the remaining
            % columns of the block
            for firstCol = 1:d
                % when the first column index is i, the elements of the remaining
                % columns can be chosen from i to d, thus there are (d - i + 1)
                % choices for these (n - 1) columns
                sizeOfBlock = SymmetricSubspace.computeDimension(d - firstCol + 1, n - 1);
                startIndex = nElementsBefore + 1;
                endIndex = startIndex + sizeOfBlock - 1;
                if rStart <= nRows && ind(rStart) <= endIndex
                    % we have a block
                    rNextStart = rStart + 1;
                    while rNextStart <= nRows && ind(rNextStart) <= endIndex
                        rNextStart = rNextStart + 1;
                    end
                    rEnd = rNextStart - 1;
                    remainingColsInd = ind(rStart:rEnd) - nElementsBefore;
                    remainingColsSub = SymmetricSubspace.indToSubSymHelper(remainingColsInd, d - firstCol + 1, n - 1, dimTable);
                    sub(rStart:rEnd, 1) = firstCol;
                    sub(rStart:rEnd, 2:end) = remainingColsSub + firstCol - 1;
                    rStart = rNextStart;
                end
                nElementsBefore = nElementsBefore + sizeOfBlock;
            end            
        end
        function ind = subToIndSymHelper(sub, d, n, dimTable)
        % helper function for subToIndSym
        % assumes that sub is already sorted, i.e.
        % each subindex is in the canonical form (increasing)
        % and the rows are sorted lexicographically
            if n == 1
                ind = sub;
                return
            end
            rStart = 1;
            nRows = size(sub, 1);
            ind = zeros(nRows, 1);
            % for speed, we treat all rows with the same "first column value" as a group
            % then, the subspace spanned by the remaining column is also a symmetric subspace
            % whose dimension is reduced by the index of the first column value (canonical indices
            % are increasing), and number of copies is n - 1.
            while rStart <= nRows
                rNextStart = rStart + 1;
                firstCol = sub(rStart, 1);
                while rNextStart <= nRows && sub(rNextStart, 1) == firstCol
                    rNextStart = rNextStart + 1;
                end
                rEnd = rNextStart - 1;
                nElementsBefore = 0;
                for firstColBefore = 1:(firstCol-1)
                    % when the first column index is i, the elements of the remaining
                    % columns can be chosen from i to d, thus there are (d - i + 1)
                    % choices for these (n - 1) columns
                    sizeOfBlock = dimTable(d - firstColBefore + 1, n - 1);
                    nElementsBefore = nElementsBefore + sizeOfBlock;
                end
                remainingCols = sub(rStart:rEnd, 2:end) - firstCol + 1;
                remainingColsInd = SymmetricSubspace.subToIndSymHelper(remainingCols, d - firstCol + 1, n - 1, dimTable);
                ind(rStart:rEnd) = nElementsBefore + remainingColsInd;
                rStart = rNextStart;
            end
        end
    end
    methods(Static)
        function rho = herm(d)
            cvx_begin set sdp
            variable rho(d, d) hermitian
            rho >= 0
            cvx_end
        end
        function dim = computeDimension(d, n)
            dim = nchoosek(n + d - 1, n);
        end
        function p = uniquePermutations(a)
        % Returns all the unique permutations of the vector a
        %
        % taken from https://www.mathworks.com/matlabcentral/newsreader/view_thread/164470
            [u, ~, J] = unique(a);
            p = u(up(J, length(a)));
            function p = up(J, n)
                ktab = histc(J,1:max(J));
                l = n;
                p = zeros(1, n);
                s = 1;
                for i=1:length(ktab)
                    k = ktab(i);
                    c = nchoosek(1:l, k);
                    m = size(c,1);
                    [t, ~] = find(~p.');
                    t = reshape(t, [], s);
                    c = t(c,:)';
                    s = s*m;
                    r = repmat((1:s)',[1 k]);
                    q = accumarray([r(:) c(:)], i, [s n]);
                    p = repmat(p, [m 1]) + q;
                    l = l - k;
                end
            end            
        end       
    end
end
