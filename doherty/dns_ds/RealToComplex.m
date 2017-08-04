function C = RealToComplex(R)
    d2 = size(R, 1);
    assert(size(R, 2) == d2);
    assert(mod(d2, 2) == 0);
    d = d2/2;
    R = reshape(R, [2 d 2 d]);
    R = permute(R, [2 4 1 3]);
    R = reshape(R, [d*d 2*2]);
    CR = R*[1;0;0;1]/2;
    CI = R*[0;-1;1;0]/2;
    C = reshape(CR, [d d]) + 1i*reshape(CI, [d d]);
end
