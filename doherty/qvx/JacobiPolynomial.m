classdef JacobiPolynomial
% JacobiPolynomial Represents a Jacobi polynomial
%
% Modern object-oriented formulation with parts taken from:
%
% - the LGPL library by John Burkardt
%   https://people.sc.fsu.edu/~jburkardt/m_src/jacobi_polynomial/jacobi_polynomial.html
% - QETLAB (2-Clause BSD)
%
% The LGPL of Burkardt prevails.
   
    properties(SetAccess = immutable)
        n;     % degree of the Jacobi polynomial
        alpha; % first real parameter
        beta;  % second real parameter
    end
    
    methods
        function J = JacobiPolynomial(n, alpha, beta)
            assert(n >= 0);
            J.n = n;
            J.alpha = alpha;
            J.beta = beta;
        end
        function c = coeffs(J)
        % QETLAB
            switch J.n
              case 0
                c = 1;
              case 1
                c = [J.alpha+J.beta+2, J.alpha-J.beta]/2;
              otherwise
                a = J.alpha;
                b = J.beta;
                n = J.n;
                prev_c1 = JacobiPolynomial(n - 1, a, b).coeffs;
                prev_c2 = JacobiPolynomial(n - 2, a, b).coeffs;
                c = ((2*n+a+b-1)*(a^2-b^2)*[0,prev_c1] + (2*n+a+b-1)*(2*n+a+b)*(2*n+a+b-2)*[prev_c1,0] ...
                     - 2*(n+a-1)*(n+b-1)*(2*n+a+b)*[0,0,prev_c2]) / (2*n*(n+a+b)*(2*n+a+b-2));
            end
        end
        function x = zeros(J)
            if J.alpha > -1 && J.beta > -1
                x = zerosByEigenvalues(J);
            else 
                x = roots(J.coeffs);
            end
            x = sort(x);
        end
        function x = zerosByRoots(J)
            x = roots(J.coeffs);
        end
        function x = zerosByEigenvalues(J)
            assert(J.alpha > -1 && J.beta > -1);
            % Reference of original code
            % Sylvan Elhay, Jaroslav Kautsky,
            % Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
            % Interpolatory Quadrature,
            % ACM Transactions on Mathematical Software,
            % Volume 13, Number 4, December 1987, pages 399-415.
            
            n = J.n;
            alpha = J.alpha;
            beta = J.beta;
            
            ab = alpha + beta;
            abi = 2.0 + ab;
            %
            %  Define the zero-th moment.
            %
            zemu = 2.0^( ab + 1.0 ) * gamma ( alpha + 1.0 ) ...
                   * gamma ( beta + 1.0 ) / gamma ( abi );
            %
            %  Define the Jacobi matrix.
            %
            x = zeros ( n, 1 );
            bj = zeros ( n, 1 );

            x(1) = ( beta - alpha ) / abi;
            bj(1) = 4.0 * ( 1.0 + alpha ) * ( 1.0 + beta ) ...
                    / ( ( abi + 1.0 ) * abi * abi );
            a2b2 = beta * beta - alpha * alpha;

            for i = 2 : n
                abi = 2.0 * i + ab;
                x(i) = a2b2 / ( ( abi - 2.0 ) * abi );
                abi = abi^2;
                bj(i) = 4.0 * i * ( i + alpha ) * ( i + beta ) * ( i + ab ) ...
                        / ( ( abi - 1.0 ) * abi );
            end
            bj(1:n) =  sqrt ( bj(1:n) );

            w = zeros ( n, 1 );
            w(1) = sqrt ( zemu );
            %
            %  Diagonalize the Jacobi matrix.
            %
            [ x, w ] = JacobiPolynomial.imtqlx ( n, x, bj, w );
        end
    end
   
    methods(Static)
        
        function [ d, z ] = imtqlx ( n, d, e, z )
        %*****************************************************************************80
        %
        %% IMTQLX diagonalizes a symmetric tridiagonal matrix.
        %
        %  Discussion:
        %
        %    This routine is a slightly modified version of the EISPACK routine to
        %    perform the implicit QL algorithm on a symmetric tridiagonal matrix.
        %
        %    The authors thank the authors of EISPACK for permission to use this
        %    routine.
        %
        %    It has been modified to produce the product Q' * Z, where Z is an input
        %    vector and Q is the orthogonal matrix diagonalizing the input matrix.
        %    The changes consist (essentialy) of applying the orthogonal transformations
        %    directly to Z as they are generated.
        %
        %  Licensing:
        %
        %    This code is distributed under the GNU LGPL license.
        %
        %  Modified:
        %
        %    04 January 2010
        %
        %  Author:
        %
        %    Original FORTRAN77 version by Sylvan Elhay, Jaroslav Kautsky.
        %    MATLAB version by John Burkardt.
        %
        %  Reference:
        %
        %    Sylvan Elhay, Jaroslav Kautsky,
        %    Algorithm 655: IQPACK, FORTRAN Subroutines for the Weights of
        %    Interpolatory Quadrature,
        %    ACM Transactions on Mathematical Software,
        %    Volume 13, Number 4, December 1987, pages 399-415.
        %
        %    Roger Martin, James Wilkinson,
        %    The Implicit QL Algorithm,
        %    Numerische Mathematik,
        %    Volume 12, Number 5, December 1968, pages 377-383.
        %
        %  Parameters:
        %
        %    Input, integer N, the order of the matrix.
        %
        %    Input, real D(N), the diagonal entries of the matrix.
        %
        %    Input, real E(N), the subdiagonal entries of the
        %    matrix, in entries E(1) through E(N-1). 
        %
        %    Input, real Z(N,1), a vector to be operated on.
        %
        %    Output, real D(N,1), the diagonal entries of the diagonalized matrix.
        %
        %    Output, real Z(N,1), the value of Q' * Z, where Q is the matrix that 
        %    diagonalizes the input symmetric tridiagonal matrix.
        %
            itn = 30;

            prec = eps;

            if ( n == 1 )
                return
            end

            e(n) = 0.0;

            for l = 1 : n

                j = 0;

                while ( 1 )

                    for m = l : n

                        if ( m == n )
                            break
                        end

                        if ( abs ( e(m) ) <= prec * ( abs ( d(m) ) + abs ( d(m+1) ) ) )
                            break
                        end

                    end

                    p = d(l);

                    if ( m == l )
                        break
                    end

                    if ( j == itn )
                        fprintf ( 1, '\n' );
                        fprintf ( 1, 'IMTQLX - Fatal error!\n' );
                        fprintf ( 1, '  Iteration limit exceeded.\n' );
                        error ( 'IMTQLX - Fatal error!' );
                    end

                    j = j + 1;
                    g = ( d(l+1) - p ) / ( 2.0 * e(l) );
                    r =  sqrt ( g * g + 1.0 );
                    g = d(m) - p + e(l) / ( g + JacobiPolynomial.r8_sign ( g ) * abs ( r ) );
                    s = 1.0;
                    c = 1.0;
                    p = 0.0;
                    mml = m - l;

                    for ii = 1 : mml

                        i = m - ii;
                        f = s * e(i);
                        b = c * e(i);

                        if ( abs ( f ) >= abs ( g ) )
                            c = g / f;
                            r =  sqrt ( c * c + 1.0 );
                            e(i+1) = f * r;
                            s = 1.0 / r;
                            c = c * s;
                        else
                            s = f / g;
                            r =  sqrt ( s * s + 1.0 );
                            e(i+1) = g * r;
                            c = 1.0 / r;
                            s = s * c;
                        end

                        g = d(i+1) - p;
                        r = ( d(i) - g ) * s + 2.0 * c * b;
                        p = s * r;
                        d(i+1) = g + p;
                        g = c * r - b;
                        f = z(i+1);
                        z(i+1) = s * z(i) + c * f;
                        z(i) = c * z(i) - s * f;

                    end

                    d(l) = d(l) - p;
                    e(l) = g;
                    e(m) = 0.0;

                end

            end

            for ii = 2 : n

                i = ii - 1;
                k = i;
                p = d(i);

                for j = ii : n
                    if ( d(j) < p )
                        k = j;
                        p = d(j);
                    end
                end

                if ( k ~= i )
                    d(k) = d(i);
                    d(i) = p;
                    p = z(i);
                    z(i) = z(k);
                    z(k) = p;
                end

            end

            return
        end

        function value = r8_sign ( x )
        %*****************************************************************************80
        %
        %% R8_SIGN returns the sign of an R8.
        %
        %  Discussion:
        %
        %    The value is +1 if the number is positive or zero, and it is -1 otherwise.
        %
        %  Licensing:
        %
        %    This code is distributed under the GNU LGPL license.
        %
        %  Modified:
        %
        %    21 March 2004
        %
        %  Author:
        %
        %    John Burkardt
        %
        %  Parameters:
        %
        %    Input, real X, the number whose sign is desired.
        %
        %    Output, real VALUE, the sign of X.
        %
            if ( 0 <= x )
                value = +1.0;
            else
                value = -1.0;
            end

            return
        end
    end
    
end
