function H = OperatorFromCoeffs2(coeffs)
    dA = sqrt(size(coeffs, 1));
    dB = sqrt(size(coeffs, 2));
    n = (dA^2)*(dB^2);
    H = zeros(dA*dB, dA*dB);
    FA = GeneralizedGellMann(dA);
    FB = GeneralizedGellMann(dB);
    for a = 1:dA^2
        for b = 1:dB^2
            H = H + coeffs(a, b) * kron(FA{a}, FB{b});
        end
    end
end
