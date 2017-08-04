function indices = SymmetricCanonicalIndices(d, k, indexd)
% SymmetricCanonicalIndices Returns canonical index representative for a k-d symmetric tensor
%
% Consider a k-dimensional tensor T(i1,i2...ik), where i1..ik = 1...d
% This function generates all k-tuples (i1,i2...ik) such that i1 <= i2 <= ... <= ik
% and returns a list of the corresponding indices by sub2ind(d*ones(1,k), i1,i2,...,ik)
    if nargin == 2
        indexd = d;
    end
    switch k
      case 0
        indices = 1;
      case 1
        indices = 1:d;
      otherwise
        indices = [];
        for lastIndex = 1:d
            subindices = SymmetricCanonicalIndices(lastIndex, k - 1, indexd);
            lastIndexWeight = (lastIndex-1)*indexd^(k-1);
            newIndices = subindices + lastIndexWeight;
            indices = [indices newIndices];
        end
    end
end
