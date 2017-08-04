function B = OrthoIntegerBasis(d)
% OrthoIntegerBasis(d)
%
% For an integer d, returns an integer matrix M such that M is full rank,
% and M'*M is diagonal.
    if d == 1
        B = [1];
        return
    end
    f = factor(d);
    B = PrimitiveBasis(f(1));
    for i = 2:length(f)
        B = kron(B, PrimitiveBasis(f(i)));
    end
    function Bp = PrimitiveBasis(p)
        Bp = zeros(p, p);
        Bp(:,1) = 1;
        for i = 1:p-1
            Bp(i,i+1) = p-i;
            Bp(i+1:end,i+1) = -1;
        end
    end
end
