function R = ComplexToReal(C)
% For a complex Hermitian matrix C of dimension dxd, returns a real
% matrix R of dimension d'xd' with d' = 2d, where R has the same
% eigenvalues as C.
    R = kron(real(C), [1 0; 0 1]) + kron(imag(C), [0 -1; 1 0]);
end
