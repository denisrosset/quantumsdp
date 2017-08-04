% requires CVX, SeDuMi and the shim "diagout"

% we do not normalize the state to get integer coefficients
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

for ext = 2:7
    def = SeparableConeDef([3 3], 'outer', ext, 'ppt', 'doherty');
    cvx_solver diagout
    cvx_solver_settings('dumpfile', ['../qvx_doherty_' num2str(ext) '.mat'])
cvx_begin
variable v
    rhov = 2*psiplus + v * sigmaplus + (5-v)*VsigmaplusV;
    minimize v
    subject to
    rhov == SeparableConeC(def);
    cvx_problem
cvx_end
    
end
