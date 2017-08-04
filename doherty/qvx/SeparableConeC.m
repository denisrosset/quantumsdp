function rho = SeparableConeC(def)
% SeparableConeC Exact/approx. formulation of the cone of separable operators
%
% INPUTS
% def        Separable cone exact/approx. definition, see SeparableConeDef
%
% OUTPUTS
% rho        Complex (dA*dB)x(dA*dB) matrix representing the AB system
%            The A,B basis ordering is such that a product state
%            rho = rhoA (x) rhoB = kron(rhoA, rhoB)
    
    dims = def.dims;
    dA = dims(1);
    dB = dims(2);
    d = dA*dB;
    k = def.k;
    
    if k == 1 % handle separately the PPT criterion alone, without extension
        d = dA*dB;
        if length(def.pptCuts) == 0
            warning('The symmetric 1-extension without PPT constraint is trivial');
            cvx_begin set sdp
            variable rho(d, d) hermitian % main variable
            rho >= 0 % semidefinite positive
            cvx_end
        else
            cvx_begin set sdp
            variable rho(d, d) hermitian % main variable
            variable ppt1(d, d) hermitian % partial transpose
            rho >= 0
            ppt1 >= 0
            [Arho Appt] = def.ConstraintPPT;
            Arho * rho(:) == Appt * ppt1(:)
            cvx_end
        end
        return
    end
    
    tauS = SymmetricSubspace(dB, k);
    dBsym = tauS.dim;
    
    pptvar = cell(1, k);
    cvx_begin set sdp
    variable rho(dB*dA, dB*dA) hermitian % main variable
    variable tauSym(dBsym*dA, dBsym*dA) hermitian  % symmetric extension in symmetric basis
    for p = def.pptCuts
        k1 = p;
        k2 = k - p;
        dBsym1 = SymmetricSubspace(dB, k1).dim;
        dBsym2 = SymmetricSubspace(dB, k2).dim;
        varName = ['ppt' num2str(p)];
        varDecl = [varName '(dBsym1*dBsym2*dA, dBsym1*dBsym2*dA)'];
        variable(varDecl, 'hermitian');
        pptvar{p} = eval(varName);
    end
    for p = def.pptCuts
        pptvar{p} >= 0
        ppt = pptvar{p};
        [ApptSym AtauSym] = def.ConstraintPPTCutSym(p);
        ApptSym * ppt(:) == AtauSym * tauSym(:)
    end
    %    rho >= 0
    tauSym >= 0
    [AtauSym Arho ArhoAIdB] = def.ConstraintRepresentsSym;
    if isequal(def.approx, 'inner')
        epsN = def.innerFactor;
        (dB*(1-epsN) * Arho + epsN * ArhoAIdB) * rho(:) == dB * AtauSym * tauSym(:)
    else
        Arho * rho(:) == AtauSym * tauSym(:)
    end
    cvx_end
end
