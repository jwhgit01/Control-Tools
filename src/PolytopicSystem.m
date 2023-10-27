classdef PolytopicSystem
    %PolytopicSystem Summary of this class goes here
    %
    % Copyright (c) 2023 Jeremy W. Hopwood. All rights reserved.
    %
    % This class ...
    %
    % Properties:
    %
    %   Systems     A 1xN cell array of state space objects, uncertain
    %               state space objects, or polytopic system objects.
    %
    %   Vetices
    %
    %   Dimensions
    %
    %   ParameterNames
    %
    %   Legacy
    %  
    % Methods:
    %
    %   method1    The first output...
    %

    properties
        Systems
        Vertices
        LinearizationState
        LinearizationInput
        ParameterNames
        Legacy
    end

    methods
        function obj = PolytopicSystem(Systems,Vertices)
            %PolytopicSystem Construct an instance of this class
            %   Detailed explanation goes here

            % Create an empty object if no inputs are given
            if nargin < 1
                return
            end

            % Assign just the vertices if Systems is empty.
            if isempty(Systems)
                obj.Systems = {};
                obj.Vertices = Vertices;
                return
            end

            % Ensure inputs are given correctly.
            if ~iscell(Systems)
                error('Systems must be given as a 1xN cell array.');
            end
            if all(size(Systems)~=1)
                error('Systems cell array must be a one-dimensional array.');
            end
            if length(Systems) ~= size(Vertices,2)
                error('Number of systems must equal number of vertices.');
            end
            
            % Set object properties.
            obj.Systems = Systems;
            obj.Vertices = Vertices;

        end % PolytopicSystem

        function obj = AddVertex(obj,System,Vertex)
            %AddVertex
            obj.Systems{end+1} = System;
            obj.Vertices(:,end+1) = Vertex;
        end % AddVertex

        function [ny,nu,nx,np,N] = size(obj)
            % size Determine dimensions of the system based on its type.
            mc = metaclass(obj.Systems{1});
            if strcmp(mc.Name,'PolytopicSystem')
                sys0 = obj.Systems{1}.Systems{1};
            elseif strcmp(mc.Name,'ss')
                sys0 = obj.Systems{1};
            else
                error('Invalid vertex system type.')
            end
            ny = size(sys0.C,1);
            nu = size(sys0.B,2);
            nx = size(sys0.A,1);
            np = size(obj.Vertices,1);
            N = size(obj.Vertices,2);
        end % size

        function sys = Evaluate(obj,method,arg)
            %Evaluate
            %
            % If method = 'Nominal', loop through each LPV vertex and
            % evaluate at the center of the uncertainty polytope.
            %
            % If method = 'Gaussian', loop through each LPV vertex and
            % sample each uncertain parameter from a Gaussian distribution
            % with standard deviation equal to 1/arg the 'plus or minus'
            % bound on the parameter.
            %
            % If method = 'Uniform', loop through each LPV vertex and
            % sample each uncertain parameter from a uniform distribution
            % with bounds equal to arg the parameter bounds.
            %
            % If method = 'Parameter', evaluate the system for a parameter
            % value of arg. If the vertex systems are polytopic systems
            % themselves, convert each vertex to a uncertain state space
            % object.


            % Evaluate the system based on the given method.
            if strcmp(method,'Nominal')

                % Initialize the output
                sys = obj;

                % Get the dimensions of the system.
                [ny,nu,nx,~,N] = size(obj);

                % Loop through the LPV polytope vertices.
                for ii = 1:N

                    % Initialize the LTI vertex system
                    A = zeros(nx,nx);
                    B = zeros(nx,nu);
                    C = zeros(ny,nx);
                    D = zeros(ny,nu);

                    % Polytopic uncertain vertex system and 
                    sysi = obj.Systems{ii};

                    % If this vertex is not uncertain, skip it.
                    mc = metaclass(sysi);
                    if strcmp(mc.Name,'ss')
                        continue
                    end

                    % Number of uncertainty polytope vertices.
                    Mi = size(sysi.Vertices,2);
        
                    % Coordinates of the center of the polytope.
                    c = (1/Mi)*ones(1,Mi);

                    % Evaluate the system using polytope coordinates.
                    for jj = 1:Mi
                        sysj = sysi.Systems{jj};
                        A = A + c(jj)*sysj.A;
                        B = B + c(jj)*sysj.B;
                        C = C + c(jj)*sysj.C;
                        D = D + c(jj)*sysj.D;
                    end
                    sys.Systems{ii} = ss(A,B,C,D);

                end

            elseif strcmp(method,'Gaussian') || strcmp(method,'Uniform')

                % Initialize the output
                sys = obj;

                % Get the dimensions of the system.
                [ny,nu,nx,~,N] = size(obj);

                % Loop through the LPV polytope vertices.
                for ii = 1:N

                    % Initialize vertex system
                    A = zeros(nx,nx);
                    B = zeros(nx,nu);
                    C = zeros(ny,nx);
                    D = zeros(ny,nu);

                    % Polytopic uncertain vertex system and 
                    sysi = obj.Systems{ii};

                    % If this vertex is not uncertain, skip it.
                    mc = metaclass(sysi);
                    if strcmp(mc.Name,'ss')
                        continue
                    end

                    % Number of uncertainty polytope vertices and
                    % uncertain parameters.
                    Mi = size(sysi.Vertices,2);
                    nth = size(sysi.Vertices,1);

                    % Obtain min and max values of each parameter.
                    pmin = min(sysi.Vertices,[],2);
                    pmax = max(sysi.Vertices,[],2);
        
                    % Loop through parameters to get polytopic coordinates
                    for jj = 1:nth

                        % parameter bounds of jth parameter
                        pjbar = (pmax(jj)+pmin(jj))/2;
                        pjrange = pmax(jj)-pmin(jj);

                        % Uncertain parameter value
                        if strcmp(method,'Gaussian')
                            pj = (1/arg)*(pjrange/2)*randn + pjbar;
                        else
                            pj = arg*pjrange*rand - pjbar;
                        end
                               
                        % Linear interpolation to obtain polytopic
                        % coordinates.
                        z = (pj-pmin(jj))/(pmax(jj)-pmin(jj));
                        if jj == 1
                            c = [1-z, z];
                        else
                            c = [c*(1-z), c*z];
                        end
        
                    end

                    % Evaluate vertex system using polytope coordinates.
                    for jj = 1:Mi
                        sysj = sysi.Systems{jj};
                        A = A + c(jj)*sysj.A;
                        B = B + c(jj)*sysj.B;
                        C = C + c(jj)*sysj.C;
                        D = D + c(jj)*sysj.D;
                    end
                    sys.Systems{ii} = ss(A,B,C,D);

                end

            elseif strcmp(method,'Parameter')

                % Get the dimensions of the system.
                [ny,nu,nx,np,N] = size(obj);

                % Loop through the LPV polytope vertices.
                for ii = 1:N

                    % Vertex system and number of vertices.
                    sysi = obj.Systems{ii};

                    % If the vertex system is a state space or uncertain
                    % state space object, get matrices. If it is polytopic
                    % uncertain, convert to uss object.
                    mc = metaclass(sysi);
                    if strcmp(mc.Name,'ss') || strcmp(mc.Name,'uss')

                        % Get matrices.
                        A = sysi.A;
                        B = sysi.B;
                        C = sysi.C;
                        D = sysi.D;

                    elseif strcmp(mc.Name,'PolytopicSystem')
                        warning('For a large number of vertices, this may take a while...');
                        warning('Your should evaluate the vertex system polytopes before evaluating the LPV polytope.');

                        % Number of uncertainty polytope vertices and
                        % parameters.
                        Mi = size(sysi.Vertices,2);
                        nth = size(sysi.Vertices,1);

                        % Initialize vertex system
                        A = umat(zeros(nx,nx));
                        B = umat(zeros(nx,nu));
                        C = umat(zeros(ny,nx));
                        D = umat(zeros(ny,nu));

                        % Compute parameter bounds and nominal value.
                        pmin = min(sysi.Vertices,[],2);
                        pmax = max(sysi.Vertices,[],2);
                        pbar = (pmax+pmin)/2;
                        plusminus = (pmax-pmin)/2;

                        % Loop through parameters and linearly interpolate
                        % using the ureal parameter to obtain the polytopic
                        % coordinates.
                        for jj = 1:nth

                            % Create ureal parameter.
                            pname = ['theta_' num2str(ii) '_' num2str(jj)];
                            pj = ureal(pname,pbar(jj),'PlusMinus',plusminus(jj));

                            % Linearly interpolate
                            z = (pj-pmin(jj))/(pmax(jj)-pmin(jj));
                            if jj == 1
                                c = [1-z, z];
                            else
                                c = [c*(1-z), c*z];
                            end

                        end

                        % Evaluate vertex system using polytope coordinates.
                        for jj = 1:Mi
                            sysj = sysi.Systems{jj};
                            A = A + c(jj)*sysj.A;
                            B = B + c(jj)*sysj.B;
                            C = C + c(jj)*sysj.C;
                            D = D + c(jj)*sysj.D;
                        end

                    end

                    % Construct the vertex LTI state space system.
                    obj.Systems{ii} = ss(A,B,C,D);

                end

                % With the vertex systems of the LPV polytope specificed as
                % either state space or uncertain state space objects,
                % compute the LPV poyltope coordinates for the given
                % parameter value, p=arg.
                pmin = min(obj.Vertices,[],2);
                pmax = max(obj.Vertices,[],2);
                for ii = 1:np
                    z = (arg(ii)-pmin(ii))/(pmax(ii)-pmin(ii));
                    if ii == 1
                        c = [1-z, z];
                    else
                        c = [c*(1-z), c*z];
                    end
                end

                % Evaluate vertex system using polytope coordinates.
                x0 = zeros(nx,1);
                u0 = zeros(nu,1);
                A = zeros(nx,nx);
                B = zeros(nx,nu);
                C = zeros(ny,nx);
                D = zeros(ny,nu);
                for ii = 1:N
                    sysi = obj.Systems{ii};
                    x0 = x0 + c(ii)*obj.LinearizationState(:,ii);
                    u0 = u0 + c(ii)*obj.LinearizationInput(:,ii);
                    A = A + c(ii)*sysi.A;
                    B = B + c(ii)*sysi.B;
                    C = C + c(ii)*sysi.C;
                    D = D + c(ii)*sysi.D;
                     
                end
                sys = ss(A,B,C,D);
                sys.UserData.x0 = x0;
                sys.UserData.u0 = u0;

            else
                error('Invalid evaluation method.');
            end

        end % Evaluate

    end % public methods

    methods (Static)
        
        function obj = uss2pol(usys)
            %uss2pol
            %
            % Factory for creating a polytopic system from a uss object.
            %
            % Rules:
            %   - system must be affine in the parameters
            %   - more than one parameter may appear in a matrix element
            
            % Dimensions
            nx = size(usys.A,1);
            nu = size(usys.B,2);
            ny = size(usys.C,1);
            np = size(fieldnames(usys.Uncertainty),1);
            
            % Get the names and ranges of these parameters
            unc = struct2cell(usys.Uncertainty);
            pNames = cell(np,1);
            range = cell(np,1);
            for ii = 1:np
                pNames{ii} = unc{ii}.Name;
                range{ii} = unc{ii}.Range;
            end
            
            % Convert ureal matrices to affine form
            Ap = umat2aff(usys.A,pNames);
            Bp = umat2aff(usys.B,pNames);
            Cp = umat2aff(usys.C,pNames);
            Dp = umat2aff(usys.D,pNames);
            
            % Loop though parameters and create ltisys "objects"
            affSysCell = cell(np+1,1);
            affSysCell(:,:) = {zeros(nx+ny+1,nx+nu+1)};
            affSys = cell2struct(affSysCell,cat(1,'o',pNames));
            affSys.o = ltisys(Ap.o,Bp.o,Cp.o,Dp.o);
            for rr = 1:np
                affSys.(pNames{rr}) = ltisys(Ap.(pNames{rr}),...
                                            Bp.(pNames{rr}),...
                                            Cp.(pNames{rr}),...
                                            Dp.(pNames{rr}));
            end
            
            % Specify range of uncertain parameters for legacy functions
            pv = pvec('box',cell2mat(range));
            
            % Create affine system using legacy functions
            affSysCell = struct2cell(affSys);
            affSysList = [affSysCell{:}];
            affSys = psys(pv,affSysList);
            
            % Transform to polytopic system using legacy functions
            Legacy.PV = pv;
            Legacy.Vertex = polydec(pv);
            Legacy.AffSys = affSys;
            Legacy.PolSys = aff2pol(affSys);
            
            % Put the legacy result into a readable format
            Vertices = Legacy.Vertex;
            N = size(Vertices,2);
            Systems = cell(1,N);
            for rr = 1:N
                sysMat = real(Legacy.PolSys(1:nx+ny,(2*rr)+(rr-1)*(nx+nu)+1:(2*rr)+(rr)*(nx+nu)));
                Ar = sysMat(1:nx,1:nx);
                Br = sysMat(1:nx,nx+1:nx+nu);
                Cr = sysMat(nx+1:nx+ny,1:nx);
                Dr = sysMat(nx+1:nx+ny,nx+1:nx+nu);
                Systems{1,rr} = ss(Ar,Br,Cr,Dr);
            end

            % Call constructor
            obj = PolytopicSystem(Systems,Vertices);

            % Also store legacy desciption
            obj.Legacy = Legacy;

        end % uss2pol

    end % static methods

end % classdef