Doherty qutrit-qutrit PPT entangled state
=========================================

Tests of Doherty hierarchy implementations on the 3x3 state given in Eq. 68 of PRA 69 022308.

- The implementation "DNS" is the one proposed in Doherty's paper, except that we use an orthogonal basis with integer coefficients (so that Eq. 15 is not satisfied exactly).

- The implementation "DS" applies the Bose symmetry trick to remove columns from the matrix of the DNS approach.

- The implementation "QVX" is inspired by the implementation in QETLAB, but uses integer coefficients, and uses the symmetric basis also on the PPT constraints.
