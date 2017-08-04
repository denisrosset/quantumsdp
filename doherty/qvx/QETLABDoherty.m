proj = @(x) x(:)*x(:)';
psiplus = [1 0 0 0 1 0 0 0 1]';
psiplus = proj(psiplus)/3;
C1 = [0 1 0
      0 0 0
      0 0 0];
C2 = [0 0 0
      0 0 1
      0 0 0];
C3 = [0 0 0
      0 0 0
      1 0 0];
sigmaplus = proj(C1(:))/3 + proj(C2(:))/3 + proj(C3(:))/3;
C1 = C1';
C2 = C2';
C3 = C3';
VsigmaplusV = proj(C1(:))/3 + proj(C2(:))/3 + proj(C3(:))/3;

for ext = 2:7
    def = SeparableConeDef([3 3], 'outer', ext, 'ppt', 'doherty');
    cvx_solver sedumi %diagout
    %    cvx_solver_settings('dumpfile', ['cvx_good_doherty_' num2str(ext) ...
    %                        '.mat'])
cvx_begin
variable v
    rhov = 2*psiplus/7 + v * sigmaplus/7 + (5-v)*VsigmaplusV/7;
    minimize v
    subject to
    rhov == SeparableConeC(def);
    cvx_problem
cvx_end
    
end
