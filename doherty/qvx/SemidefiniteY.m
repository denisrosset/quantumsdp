function C = SemidefiniteY(lme)
% SemidefiniteY Constructs the constraint lme >=_SDP 0 in a robust manner
%
% INPUT
%
% lme     (Hermitian) matrix, either double or sdpvar
%
% OUTPUT
%
% C       Constraint lme >= 0
%
% Emits a warning if YALMIP could misinterpret the constraint as element-wise instead of SDP
    if ~ishermitian(lme)
        warning('Input is not Hermitian/symmetric, verify that YALMIP formulates matrix inequalities.');
    end
    C = [lme >= 0];
end
