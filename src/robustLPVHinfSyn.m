function out = robustLPVHinfSyn(ussList,gam)
%robustLPVHinfSyn
%   Detailed explanation goes here

% dimensions
dim.LMI39 = ????;
dim.LMI40 = ????;
dim.LMI42 = ????;
dim.LMI43 = ????;

% begin semi-definite program
cvx_begin sdp
    
    % create CVX variables
    variable X (n,n) semidefinite
    for ii = 1:N
        eval(['variable G' num2str(ii) '(????,????)'])
        eval(['G(:,:,ii) = G' num2str(ii) ';']);
        for jj = 1:M
            eval(['variable Y' num2str(ii) num2str(jj) '(????,????) semidefinite'])
            eval(['Y(:,:,ii,jj) = Y' num2str(ii) num2str(jj) ';']);
        end
    end
    
    % objective variables
    variable delSq nonnegative

    % optimization
    minimize delSq
    subject to
        
        LMIS...

cvx_end

% results
out.???? = ;

end

function M = LMI39(sysij,X,alphakl,betakl,betalk)

Aij = sysij.A;
Buj = sysij.Bu;

Uij = Aij*X + Buj*Gi;

M = alphakl*X + betakl*Uij + betalk*uij';

end

function M = LMI40(sysij,X,Gi,gam,dim)

nw = dim.nw;
nzInf = dim.nzInf;

Aij = sysij.A;
Bwij = sysij.Bw;
Buj = sysij.Bu;
CzInfij = sysij.CzInf;
DzInfuij = sysij.DzInfu;
DzInfwij = sysij.DzInfw;

Uij = Aij*X + Buj*Gi;
Vij = CzInfij*X + DzInfuij*Gi;

M = [Uij+Uij',Bwij,Vij';
     Bwij',-eye(nw),DzInfwij';
     Vij,DzInfwij,-gam^2*eye(nzInf)];

end

function M = LMI42(Yij)

M = trace(Yij);

end

function M = LMI43(sysij,X,Yij,Gi)

Cz2ij = sysij.Cz2;
Dz2uj = sysij.Dz2u;

Wij = Cz2ij*X + Dz2uj*Gi;

M = [Yij,Wij;
     Wij',X];

end