function coeffs = CoeffsFromOperator3(H, dA, dB, dC)
    assert(size(H, 1) == dA*dB*dC);
    [FA DA] = GeneralizedGellMann(dA);
    [FB DB] = GeneralizedGellMann(dB);
    [FC DC] = GeneralizedGellMann(dC);
    switch class(H)
      case 'sdpvar'
        coeffs = sdpvar(dA^2, dB^2, dC^2, 'full', 'real');
      case 'double'
        coeffs = zeros(dA^2, dB^2, dC^2);
      case 'sym'
        coeffs = sym(zeros(dA^2, dB^2, dC^2));
      otherwise
        error(['Unsupported input class for H: ' class(H)]);
    end
    for a = 1:dA^2
        for b = 1:dB^2
            for c = 1:dC^2
                coeffs(a, b, c) = trace(kron(kron(FA{a}, FB{b}), FC{c}) * H)/DA(a)/DB(b)/DC(c);
            end
        end
    end
end
