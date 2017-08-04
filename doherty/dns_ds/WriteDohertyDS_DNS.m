% requires YALMIP and SeDuMi

proj = @(x) x(:)*x(:)';
psiplus = [1 0 0 0 1 0 0 0 1]';
psiplus = proj(psiplus);
C1 = [0 1 0
      0 0 0
      0 0 0];
C2 = [0 0 0
      0 0 1
      0 0 0];
C3 = [0 0 0
      0 0 0
      1 0 0];
sigmaplus = proj(C1(:)) + proj(C2(:)) + proj(C3(:));
C1 = C1';
C2 = C2';
C3 = C3';
VsigmaplusV = proj(C1(:)) + proj(C2(:)) + proj(C3(:));
options = sdpsettings('verbose', 1, 'solver', 'sedumi');
for ext = 2:3
for useSym = 0:1
    yalmip clear;
    v = sdpvar;
    rhov = 2*psiplus + v * sigmaplus + (5-v)*VsigmaplusV;
    coeffsv = CoeffsFromOperator2(3*rhov, 3, 3);
    Cons = SymmetricExtensionCone(coeffsv, ext, useSym, true);
    [model,recoverymodel,diagnostic,internalmodel] = export(Cons, v, options);
    At = model.A;
    b = model.b;
    c = model.C;
    K = model.K;
    if useSym
        save(['../ds_doherty_' num2str(ext) '.mat'], 'At', 'b', 'c', 'K');
    else
        save(['../dns_doherty_' num2str(ext) '.mat'], 'At', 'b', 'c', 'K');
    end
end
end

