classdef SeparableConeDef
    
    properties(SetAccess = immutable)
        
        dims;      % = [dA dB], dimensions of the subsystems composing the symmetric cone        

        dA;        % Dimension of the first subsystem A
        
        dB;        % Dimension of the second subsystem B
        
        approx;    % Either 'inner', 'outer' or 'approx'
                   %
                   % 'outer' is the Doherty-type symmetric extension
                   % 'exact' is valid only for dA*dB <= 6, and sets
                   %         k = 1 and 'ppt' = [1]
                   % 'inner' is the Navascues inner approximation, valid
                   %         when 'ppt' = []
                   %         or 'ppt' = 'navascues', that is [ceil(k/2)]

        k;         % Number of copies of subsystem B in the extension
        
        ppt;       % PPT constraints to add. Can be of the following form:
                   %  = []            no PPT constraints
                   %  = [t_1 ... t_n] with 1 <= t_j <= k (duplicates ignored)
                   %                  for each t_j, adds a PPT constraint such that a number t_j of copies
                   %                  of B is transposed
                   %  = 'doherty'     equivalent to [1 2 ... k], 
                   %                  PPT conditions in the original 2004 Doherty paper
                   %  = 'navascues'   equivalent to [ceil(k/2)], 
                   %                  PPT condition in the 2009 Navascues paper,
                   %                  also the way it is implemented in QETLAB, as of end 2016
                   % Default: [] (do not add PPT constraints)
        
        pptCuts;   % integer vector enumerating the cuts corresponding to the PPT constraints
                   % When ppt = 'doherty' or 'navascues', see interpretation above,
                   % otherwise pptCuts = ppt.
    end
    
    methods
        
        function def = SeparableConeDef(dims, approx, k, varargin)
        % SeparableConeDef Defines a exact/approx. formulation of the separable cone
        %
        % INPUTS
        %
        % dims       = [dA dB] Dimension of subsystems A and B
        % approx     = 'exact', 'outer', 'inner', see property definitions above
        % k          = number of copies of subsystem B in the extension
        % 
        % Additional properties are given each by a key/value pair. Possible keys
        % are 'ppt', as explained in the 'properties' section of
        % 'SeparableConeDef.m'.
        %
        % (Right now, we don't have additional properties; we could
        % add the possibility to work in a nonsymmetrized basis for
        % easier verification).
        %
        % Examples
        %
        % The symmetric extension considered in Doherty 2004, Section VII.A, consists
        % of the symmetric extension of a qutrit-qutrit state, using two copies of B, 
        % with two additional PPT constraints
        %
        % def = SeparableConeDef([3 3], 'outer', 2, 'ppt', 'doherty')
        %
        % which is equivalent to
        %
        % def = SeparableConeDef([3 3], 'outer', 2, 'ppt', [1 2])
        %
        % To use a formulation without symmetry handling, write:
        %   
        % def = SeparableConeDef([3 3], 'outer', 2, 'ppt', 'doherty')
        %
        % To obtain the PPT exact formulation for dA * dB <= 6, write:
        % def = SeparableConeDef([dA dB], 'exact')
        %
        % References
        % Navascues 2009, DOI: 10.1103/PhysRevLett.103.160404
        % Doherty 2004, DOI: 10.1103/PhysRevA.69.022308
            assert(length(dims) == 2);
            def.dims = dims(:)';
            def.dA = def.dims(1);
            def.dB = def.dims(2);
            if nargin < 2
                if def.dA * def.dB > 6 && def.dA > 1 && def.dB > 1
                    error('No exact formulations exist for dA*dB > 6. Specify approximation type.');
                else
                    approx = 'exact';
                end
            end
            switch approx
              case 'exact'
                if def.dA * def.dB > 6 && def.dA > 1 && def.dB > 1
                    error('Exact formulations of the symmetric cone are valid for 2x2 and 2x3 states');
                end
                if nargin < 3
                    k = 1;
                end
                def.k = k;
                def.approx = approx;
                def.ppt = [1];
                def.pptCuts = [1];
              otherwise
                def.k = k;
                def.approx = approx;
                def.ppt = [];
                def.pptCuts = [];
            end
            assert(def.k >= 1);
            
            % Interpret varargin

            assert(mod(length(varargin), 2) == 0, 'Each key should have an associated value');

            toString = @(var) evalc(['disp(var)']);

            for keyIndex = 1:2:length(varargin)
                key = varargin{keyIndex};
                value = varargin{keyIndex + 1};
                switch key
                  case 'ppt'
                    def.ppt = value;
                    if isequal(value, 'doherty')
                        def.pptCuts = 1:def.k;
                    elseif isequal(value, 'navascues')
                        def.pptCuts = ceil(def.k/2);
                    else
                        def.pptCuts = unique(value(:)');
                        assert(all(def.pptCuts >= 1) && all(def.pptCuts <= def.k));
                    end
                  otherwise
                    warning(['Unsupported option: ' toString(key)]);
                end
            end

            % verifies the sanity of the definition
            switch def.approx
              case 'exact'
                if length(def.pptCuts) == 0
                    error('Exact formulations of the symmetric cone require PPT constraints.');
                end
              case 'inner'
                if length(def.pptCuts) == 0
                    if def.k == 1
                        error('The inner approximation requires either PPT constraints or k > 1.')
                    end
                elseif isequal(def.pptCuts, ceil(def.k/2))
                    % good
                else
                    error('The inner approximation is only valid for ''ppt'' = [] or ''navascues''');
                end
            end
        end
        
        function epsN = innerFactor(def)
        % Correction of Eqs. (3) and (4) of DOI:10.1103/PhysRevLett.103.160404
            if length(def.pptCuts) == 0
                epsN = 1/(def.k + def.dB);
            else
                Jn = floor(def.k/2) + 1;
                Ja = def.dB - 2;
                Jb = mod(def.k, 2);
                epsN = (1 - max(JacobiPolynomial(Jn, Ja, Jb).zeros))*def.dB/(2*(def.dB-1));
            end            
        end
        
        function [AtauSym Arho ArhoAIdB] = ConstraintRepresentsSym(def)
        % equality constraints that express that
        % tauFull(dB^k*dA, dB^k*dA) == rho(dB*dA, dB*dA)
        % with tau represented in the symmetric subspace as tauSym
        %
        % the constraint is Arho * rho(:) == AtauSym * tauSym(:)
        %
        % ArhoAIdB is such that ArhoAIdB * rho(:) is the same as
        % Arho * (rho_A (x) Id_B)
            dA = def.dA;
            dB = def.dB;
            d = dA*dB;
            k = def.k;
            tauS = SymmetricSubspace(dB, k);
            dBsym = tauS.dim;

            % (full) index subspaces AB and R = B^(k-1)
            B_A = MultiIndex([dB dA]);
            AB_AB = MultiIndex([d d]);
            B_A_B_A = MultiIndex([dB dA dB dA]);
            Bk = MultiIndex(dB*ones(1, k));
            B_R = MultiIndex([dB dB^(k-1)]);
            Bsym_A = MultiIndex([dBsym dA]);
            Bsym_A_Bsym_A = MultiIndex([dBsym dA dBsym dA]);
            
            nRest = dB^(k-1); % dimension over which we perform the partial trace
            traceOver = (1:nRest)';
            nEqs = (d+1)*d/2; % upper triangle dimension
            Arho = sparse(nEqs, d^2);
            ArhoAIdB = sparse(nEqs, d^2);
            AtauSymT = sparse(dBsym*dA*dBsym*dA, nEqs); % stores the transpose, MATLAB uses
                                                        % compressed column storage
            i = 1;
            for r = 1:d
                rBA = B_A.indToSub(r);
                rB = rBA(:,1);
                rA = rBA(:,2);
                rBfull = Bk.indToSub(B_R.subToInd([rB*ones(nRest, 1) traceOver])); % format B_R to B_B.._B
                rBsym = tauS.subToIndSym(rBfull);
                for c = r:d % restrict constraints to upper triangle
                    cBA = B_A.indToSub(c);
                    cB = cBA(:,1);
                    cA = cBA(:,2);
                    cBfull = Bk.indToSub(B_R.subToInd([cB*ones(nRest, 1) traceOver])); % format B_R to B_B.._B
                    cBsym = tauS.subToIndSym(cBfull);
                    
                    Arho(i, AB_AB.subToInd([r c])) = 1;
                    if rB == cB
                        indices = B_A_B_A.subToInd([(1:dB)' rA*ones(dB,1) (1:dB)' cA*ones(dB,1)]); 
                        ArhoAIdB(i, indices) = 1;
                    end
                    tauInd = Bsym_A_Bsym_A.subToInd([rBsym rA*ones(nRest, 1) cBsym cA*ones(nRest, 1)]);
                    tauOnes = ones(size(tauInd));
                    col = sparse(tauInd, tauOnes, tauOnes, dBsym*dA*dBsym*dA, 1);
                    AtauSymT(:, i) = col;
                    i = i + 1;
                end
            end
            AtauSym = AtauSymT.';
        end
        
        function [Appt Arho] = ConstraintPPT(def)           
        % equality constraints that express that
        % rho(dB*dA, dB*dA)^TB == ppt(dB*dA, dB*dA)
            dA = def.dA;
            dB = def.dB;
            d = dA*dB;
            matRCSym = SymmetricSubspace(d, 2); % parameterization of upper triangle
            nEqs = matRCSym.dim;
            rc = matRCSym.indToSubSym((1:nEqs)');
            BA_BA = MultiIndex([d d]);
            B_A = MultiIndex([dB dA]);
            B_A_B_A = MultiIndex([dB dA dB dA]);
            rBA = B_A.indToSub(rc(:,1));
            cBA = B_A.indToSub(rc(:,2));
            rhoInd = BA_BA.subToInd(rc);
            pptInd = B_A_B_A.subToInd([cBA(:,1) rBA(:,2) rBA(:,1) cBA(:,2)]);
            Arho = sparse(1:nEqs, rhoInd', ones(1, nEqs), nEqs, d*d);
            Appt = sparse(1:nEqs, pptInd', ones(1, nEqs), nEqs, d*d);
        end
        
        function [ApptSym AtauSym] = ConstraintPPTCutSym(def, p)
        % equality constraints that express that
        % tauFull(dB^k*dA, dB^k*dA)^T1...^Tp == sigmaFull(dB^k*dA, dB^k*dA)
        % in their respective symmetric subspaces
            dA = def.dA;
            dB = def.dB;
            d = dA*dB;
            k = def.k;
            k1 = p;
            k2 = k - p;
            sym = SymmetricSubspace(dB, k);
            dBsym = sym.dim;
            sym1 = SymmetricSubspace(dB, k1);
            sym2 = SymmetricSubspace(dB, k2);
            dB1 = sym1.dim;
            dB2 = sym2.dim;
            dPPT = dA * dB1 * dB2;
            B2_B1_A = MultiIndex([dB2 dB1 dA]);
            B2B1A_B2B1A = MultiIndex([dB2*dB1*dA dB2*dB1*dA]);
            Bsym_A_Bsym_A = MultiIndex([dBsym dA dBsym dA]);                       
            nEqs = (dPPT+1)*dPPT/2;
            ApptSymT = sparse(dPPT^2, nEqs);
            AtauSymT = sparse(dBsym*dA*dBsym*dA, nEqs);
            i = 1;
            for r = 1:dPPT
                cols = (r:dPPT)'; % upper triangle
                nRows = length(cols);
                rows = r*ones(nRows,1);
                rB2B1A = B2_B1_A.indToSub(rows);
                cB2B1A = B2_B1_A.indToSub(cols);
                rA = rB2B1A(:,3);
                cA = cB2B1A(:,3);
                rBFullInd = [sym2.indToSubSym(cB2B1A(:,1)) sym1.indToSubSym(rB2B1A(:,2))];
                cBFullInd = [sym2.indToSubSym(rB2B1A(:,1)) sym1.indToSubSym(cB2B1A(:,2))];
                rBsym = sym.subToIndSym(rBFullInd);
                cBsym = sym.subToIndSym(cBFullInd);
                blockRows = i:i+nRows-1;
                ApptSymT(B2B1A_B2B1A.subToInd([rows cols]), blockRows) = eye(nRows);
                AtauSymT(Bsym_A_Bsym_A.subToInd([rBsym rA cBsym cA]), blockRows) = eye(nRows);
                i = i + nRows;
            end
            ApptSym = ApptSymT.';
            AtauSym = AtauSymT.';
        end
        
    end % methods
    
    methods(Static)
        
    end
    
end % classdef
