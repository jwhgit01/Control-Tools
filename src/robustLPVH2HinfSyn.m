function out = robustLPVH2HinfSyn(sys,opts)
%robustLPVHinfSyn
%   Detailed explanation goes here
%
% Inputs:
%
%   syslist     An array of structs (one for each vertex of the LPV
%               polytope) each with the fields "Vertex" and "Sys"
%               describing the uncertainty polytope. "Sys" is a cell array
%               where each element is a state space object corresponding to
%               the each column of "Vertex".
%
% Outputs:
%
%   out         A struct with the following fields:
%                   field1  blah blah
%                   field2  blah blah blah

%%%%
% TODO: add support for LMI regions
%%%%

% sufficiently small number
epsil = 100*eps;

% Get robustLPVH2HinfOptions properties
dim = opts.Dim;
hinfPerf = opts.HinfPerf;
h2Perf = opts.H2Perf;

% dimensions
nx = dim.x;
nu = dim.u;
nw = dim.w;
ninf = dim.zinf;
n2 = dim.z2;
dim.LMI39 = nx;
dim.LMI40 = nx + nw + ninf;
dim.LMI43 = n2 + nx;
dim.LMI69 = 2*nx;
dim.LMI70 = 2*nx + nw + ninf;
dim.LMI73 = n2 + nx;
dim.LMI74 = 2*nx + nw;

if strcmp(opts.LyapunovFunction,'Stationary')

    % begin semi-definite program
    cvx_begin sdp
        
        % Objective variables
        if isempty(hinfPerf) && ~isempty(h2Perf)
            variable gammaSq nonnegative
            deltaSq = h2Perf^2;
            minimize gammaSq
        elseif ~isempty(hinfPerf) && isempty(h2Perf)
            gammaSq = hinfPerf^2;
            variable deltaSq nonnegative
            minimize deltaSq
        elseif ~isempty(hinfPerf) && ~isempty(h2Perf)
            gammaSq = hinfPerf^2;
            deltaSq = h2Perf^2;
        elseif isempty(hinfPerf) && isempty(h2Perf) && n2 == 0
            variable gammaSq nonnegative
            minimize gammaSq
        elseif isempty(hinfPerf) && isempty(h2Perf) && ninf == 0
            variable deltaSq nonnegative
            minimize deltaSq
        else
            error('Currently, cannot minimize both hinfPerf and h2Perf simultaneously.')
        end

        % Optimization variables
        variable X(nx,nx) semidefinite
        for ii = 1:dim.N
            eval(['variable G' num2str(ii) '(nu,nx)'])
            eval(['G(:,:,ii) = G' num2str(ii) ';']);
            if n2 > 0
                for jj = 1:dim.M(ii)
                    eval(['variable Y' num2str(ii) num2str(jj) '(n2,n2) symmetric'])
                    eval(['Y(:,:,ii,jj) = Y' num2str(ii) num2str(jj) ';']);
                end
            end
        end
    
        % Optimization constraints
        for ii = 1:dim.N
            for jj = 1:dim.M(ii)
                sysij = sys.Systems{ii}.Systems{jj};
                Gi = G(:,:,ii);
                if ninf > 0
                    opts.LMI40(sysij,X,Gi,gammaSq) <= -epsil*eye(dim.LMI40)
                end
                if n2 > 0
                    Yij = Y(:,:,ii,jj);
                    opts.LMI42(Yij) - deltaSq <= -epsil
                    opts.LMI43(sysij,X,Yij,Gi) >= epsil*eye(dim.LMI43)
                end   
            end
        end
    
    cvx_end
    
    % LPV vertex gains
    K = zeros(nu,nx,dim.N);
    for ii = 1:dim.N
        K(:,:,ii) = G(:,:,ii)/X;
    end
    
    % results
    if ninf > 0
        out.HinfPerf = sqrt(gammaSq);
    end
    if n2 > 0
        out.H2Perf = sqrt(deltaSq);
        out.Y = Y;
    end
    out.K = K;
    out.X = X;
    out.G = G;

elseif strcmp(opts.LyapunovFunction,'ParameterDependent')

    % begin semi-definite program
    cvx_begin sdp
        
        % Objective variables
        if isempty(hinfPerf) && ~isempty(h2Perf)
            variable gammaSq nonnegative
            deltaSq = h2Perf^2;
            minimize gammaSq
        elseif ~isempty(hinfPerf) && isempty(h2Perf)
            gammaSq = hinfPerf^2;
            variable deltaSq nonnegative
            minimize deltaSq
        elseif ~isempty(hinfPerf) && ~isempty(h2Perf)
            gammaSq = hinfPerf^2;
            deltaSq = h2Perf^2;
        elseif isempty(hinfPerf) && isempty(h2Perf) && n2 == 0
            variable gammaSq nonnegative
            minimize gammaSq
        elseif isempty(hinfPerf) && isempty(h2Perf) && ninf == 0
            variable deltaSq nonnegative
            minimize deltaSq
        else
            error('Currently, cannot minimize both hinfPerf and h2Perf simultaneously.')
        end

        % Optimization variables
        variable H(nx,nx)
        for ii = 1:dim.N
            eval(['variable S' num2str(ii) '(nu,nx)'])
            eval(['S(:,:,ii) = S' num2str(ii) ';']);
            for jj = 1:dim.M(ii)
                %eval(['variable XD' num2str(ii) num2str(jj) '(nx,nx) semidefinite'])
                %eval(['XD(:,:,ii,jj) = XD' num2str(ii) num2str(jj) ';']);
                if ninf > 0
                    eval(['variable Xinf' num2str(ii) num2str(jj) '(nx,nx) semidefinite'])
                    eval(['Xinf(:,:,ii,jj) = Xinf' num2str(ii) num2str(jj) ';']);
                end
                if n2 > 0
                    eval(['variable X2' num2str(ii) num2str(jj) '(nx,nx) semidefinite'])
                    eval(['X2(:,:,ii,jj) = X2' num2str(ii) num2str(jj) ';']);
                    eval(['variable Y' num2str(ii) num2str(jj) '(n2,n2) symmetric'])
                    eval(['Y(:,:,ii,jj) = Y' num2str(ii) num2str(jj) ';']);
                end
            end
        end
    
        % Optimization constraints
        for ii = 1:dim.N
            for jj = 1:dim.M(ii)
                sysij = sys.Systems{ii}.Systems{jj};
                Si = S(:,:,ii);   
                %XDij = XD(:,:,ii,jj);
                if ninf > 0
                    Xinfij = Xinf(:,:,ii,jj);
                    opts.LMI70(sysij,Xinfij,H,Si,gammaSq) <= -epsil*eye(dim.LMI70)
                end
                if n2 > 0
                    X2ij = X2(:,:,ii,jj);
                    Yij = Y(:,:,ii,jj);
                    opts.LMI42(Yij) - deltaSq <= -epsil
                    opts.LMI73(sysij,X2ij,H,Si,Yij) >= epsil*eye(dim.LMI73)
                    opts.LMI74(sysij,X2ij,H,Si) <= -epsil*eye(dim.LMI74)
                end   
            end
        end
    
    cvx_end
    
    % LPV vertex gains
    K = zeros(nu,nx,dim.N);
    for ii = 1:dim.N
        K(:,:,ii) = S(:,:,ii)/H;
    end
    
    % results
    if ninf > 0
        out.HinfPerf = sqrt(gammaSq);
        out.Xinf = Xinf;
    end
    if n2 > 0
        out.H2Perf = sqrt(deltaSq);
        out.X2 = X2;
        out.Y = Y;
    end 
    out.K = K;
    %out.XD = XD;
    out.H = H;
    out.S = S;

else
    error('''LyapunovFunction'' option must be either ''Stationary'' or ''ParameterDependent''');
end

end
