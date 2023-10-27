classdef robustLPVH2HinfOptions
    %robustLPVH2HinfOptions Summary of this class goes here
    %
    % Copyright (c) 2023 Jeremy W. Hopwood. All rights reserved.
    %
    % This class ...
    %


    properties
        HinfPerf
        H2Perf
        LyapunovFunction
        LMIRegion
        Dim
    end

    methods
        function obj = robustLPVH2HinfOptions()
            %robustLPVH2HinfOptions Construct an instance of this class
            
            % Default properties
            obj.HinfPerf = [];
            obj.H2Perf = [];
            obj.LyapunovFunction = 'Stationary'; % or 'ParameterDependent'
            obj.LMIRegion.Type = 'None'; % or 'DL' or 'DR'
            obj.LMIRegion.alpha = []; % m x m symmetric matrix
            obj.LMIRegion.beta = []; % m x m matrix
            obj.LMIRegion.Chi = []; % m x m symmetric matrix
            obj.Dim.x = [];
            obj.Dim.u = [];
            obj.Dim.w = [];
            obj.Dim.zinf = [];
            obj.Dim.z2 = [];
            obj.Dim.N = []; % number of LPV vertices
            obj.Dim.M = []; % 1xN array of numbers of uncertainty vertices

        end

        function M = LMI39(obj,sysij,X,k,l)
            alphakl = obj.LMIRegion.alpha(k,l);
            betakl = obj.LMIRegion.beta(k,l);
            betalk = obj.LMIRegion.beta(l,k);
            Aij = sysij.A;
            Buj = sysij.Bu;
            Uij = Aij*X + Buj*Gi;
            M = alphakl*X + betakl*Uij + betalk*Uij';
        end
        
        function M = LMI40(obj,sysij,X,Gi,gamSq) 
            nu = obj.Dim.u;
            nw = obj.Dim.w;
            ninf = obj.Dim.zinf;
            Aij = sysij.A;
            Buj = sysij.B(:,1:nu);
            Bwij = sysij.B(:,nu+1:nu+nw);
            Czinfij = sysij.C(1:ninf,:);
            Dzinfuij = sysij.D(1:ninf,1:nu);
            Dzinfwij = sysij.D(1:ninf,nu+1:nu+nw);
            Uij = Aij*X + Buj*Gi;
            Vij = Czinfij*X + Dzinfuij*Gi;
            M = [Uij+Uij',Bwij,Vij';
                 Bwij',-eye(nw),Dzinfwij';
                 Vij,Dzinfwij,-gamSq*eye(ninf)];
        end
        
        function M = LMI42(obj,Yij)
            % Also LMI (72)
            M = trace(Yij);
        end
        
        function M = LMI43(obj,sysij,X,Yij,Gi)
            nu = obj.Dim.u;
            ninf = obj.Dim.zinf;
            n2 = obj.Dim.z2;   
            Cz2ij = sysij.C(ninf+1:ninf+n2,:);
            Dz2uj = sysij.D(ninf+1:ninf+n2,1:nu);
            Wij = Cz2ij*X + Dz2uj*Gi;
            M = [Yij,Wij;
                 Wij',X];
        end

        function M = LMI69(obj,sysij,XDij,H,Si,k,l)
            alphakl = obj.LMIRegion.alpha(k,l);
            betakl = obj.LMIRegion.beta(k,l);
            betalk = obj.LMIRegion.beta(l,k);
            Chikl = obj.LMIRegion.Chi(k,l);
            Aij = sysij.A;
            Buj = sysij.Bu;
            Uij = Aij*H + Buj*Si;
            M = [alphakl*XDij+betakl*Uij+betalk*Uij',betalk*(XDij-H')+Chikl*Uij;
                 betakl*(XDij-H)+Chikl*Uij',Chikl*(XDij-H-H')];
        end

        function M = LMI70(obj,sysij,Xinfij,H,Si,gamSq)
            nx = obj.Dim.x;
            nu = obj.Dim.u;
            nw = obj.Dim.w;
            nzinf = obj.Dim.zinf;
            Aij = sysij.A;
            Buj = sysij.B(:,1:nu);
            Bwij = sysij.B(:,nu+1:nu+nw);
            Czinfij = sysij.C(1:nzinf,:);
            Dzinfuij = sysij.D(1:nzinf,1:nu);
            Dzinfwij = sysij.D(1:nzinf,nu+1:nu+nw);
            Uij = Aij*H + Buj*Si;
            Vij = Czinfij*H + Dzinfuij*Si;
            M = [Uij+Uij',-Xinfij+H'-Uij,Bwij,Vij';
                 -Xinfij+H-Uij',-(H+H'),zeros(nx,nw),-Vij';
                 Bwij',zeros(nw,nx),-eye(nw),Dzinfwij';
                 Vij,-Vij,Dzinfwij,-gamSq*eye(nzinf)];
        end

        function M = LMI73(obj,sysij,X2ij,H,Si,Yij)
            nu = obj.Dim.u;
            ninf = obj.Dim.zinf;
            n2 = obj.Dim.z2;   
            Cz2ij = sysij.C(ninf+1:ninf+n2,:);
            Dz2uj = sysij.D(ninf+1:ninf+n2,1:nu);
            Wij = Cz2ij*H + Dz2uj*Si;
            M = [Yij,Wij;
                 Wij',H+H'-X2ij];
        end

        function M = LMI74(obj,sysij,X2ij,H,Si)
            nx = obj.Dim.x;
            nu = obj.Dim.u;
            nw = obj.Dim.w;
            Aij = sysij.A;
            Buj = sysij.B(:,1:nu);
            Bwij = sysij.B(:,nu+1:nu+nw);
            Uij = Aij*H + Buj*Si;
            M = [Uij+Uij',-X2ij+H'-Uij,Bwij,;
                 -X2ij+H-Uij',-(H+H'),zeros(nx,nw);
                 Bwij',zeros(nw,nx),-eye(nw)];
        end

    end
end