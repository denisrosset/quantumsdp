function coeffs = CoeffsFromOperator2(H, dA, dB)
    assert(size(H, 1) == dA*dB);
    [FA DA] = GeneralizedGellMann(dA);
    [FB DB] = GeneralizedGellMann(dB);
    switch class(H)
      case 'sdpvar'
        coeffs = sdpvar(dA^2, dB^2, 'full', 'real');
      case 'double'
        coeffs = zeros(dA^2, dB^2);
      case 'sym'
        coeffs = sym(zeros(dA^2, dB^2));
      otherwise
        error(['Unsupported input class for H: ' class(H)]);
    end
    for a = 1:dA^2
        for b = 1:dB^2
            coeffs(a, b) = trace(kron(FA{a}, FB{b}) * H)/DA(a)/DB(b);
        end
    end
end
