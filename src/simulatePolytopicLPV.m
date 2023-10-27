classdef simulatePolytopicLPV
    %simulatePolytopicLPV Summary of this class goes here
    %
    % Copyright (c) 2023 Jeremy W. Hopwood. All rights reserved.
    %
    % This class ...
    %
    % It is constructed given an struct which stores a polytopic
    % representation of an LPV system. It has the folowing structure:
    %
    %


    properties
        LPVSystem
        InitialCondition
        Input
        ParameterTrajectory
        ParameterDynamics
        LinearizationPointJacobian
        TimeSpan
        HinfPerfFunction
        H2PerfFunction
        Dimensions
    end

    methods
        function obj = simulatePolytopicLPV(LPVSystem)
            %simulatePolytopicLPV Construct an instance of this class
            
            % Set properties.
            obj.LPVSystem = LPVSystem;

            % Get dimensions.
            obj.Dimensions.x = size(LPVSystem.Systems{1}.A,1);
            obj.Dimensions.u = size(LPVSystem.Systems{1}.B,2);
            obj.Dimensions.y = size(LPVSystem.Systems{1}.C,1);
            obj.Dimensions.N = size(LPVSystem.Vertices,2);
            obj.Dimensions.rho = size(LPVSystem.Vertices,1);

            % Default LinearizationPointJacobian
            obj.LinearizationPointJacobian = zeros(obj.Dimensions.x,obj.Dimensions.rho);

        end % simulatePolytopicLPV

        function [t,x,y] = Simulate(obj)
            %Simulate
            
            % Initial condition.
            x0 = obj.InitialCondition;
            if isempty(x0)
                x0 = zeros(obj.Dimensions.x,1);
            end

            % Time span.
            tspan = obj.TimeSpan;
            if isempty(tspan)
                tspan = [0 10];
            end

            % Simulate dynamics using ode45
            [t,x] = ode45(@obj.Dynamics,tspan,x0);

            % Compute output
            y = zeros(length(t),obj.Dimensions.y);
            for k = 1:length(t)
                y(k,:) = Output(obj,t(k,1),x(k,:).');
            end

        end  % Simulate

        function dxdt = Dynamics(obj,t,x)
            %Dynamics

            % Parameter trajectory.
            if isempty(obj.ParameterTrajectory)
                rho = zeros(obj.Dimensions.rho,1);
            elseif isnumeric(obj.ParameterTrajectory)
                rho = obj.ParameterTrajectory;
            else
                rho = obj.ParameterTrajectory(t,x);
            end

            % Parameter dynamics.
            if isempty(obj.ParameterDynamics)
                drhodt = zeros(obj.Dimensions.rho,1);
            elseif isnumeric(obj.ParameterDynamics)
                drhodt = obj.ParameterDynamics;
            else
                drhodt = obj.ParameterDynamics(t,x);
            end
            
            % System inputs.
            if isempty(obj.Input)
                u = zeros(obj.Dimensions.u,1);
            elseif isnumeric(obj.Input)
                u = obj.Input;
            else
                u = obj.Input(t,x);
            end

            % Evaluate polytopic system
            sys = obj.LPVSystem.Evaluate('Parameter',rho);

            % State and input perturbation
            x0 = sys.UserData.x0;
            u0 = sys.UserData.u0;
            dx = x - x0;
            du = u - u0;

            % Linearization point dynamics
            dx0drho = obj.LinearizationPointJacobian;
            dx0dt = dx0drho*drhodt;

            % LPV system dynamics
            dxdt = sys.A*dx + sys.B*du + dx0dt;

        end % Dynamics

        function y = Output(obj,t,x)
            %Output

            % Parameter trajectory.
            if isempty(obj.ParameterTrajectory)
                rho = zeros(obj.Dimensions.rho,1);
            elseif isnumeric(obj.ParameterTrajectory)
                rho = obj.ParameterTrajectory;
            else
                rho = obj.ParameterTrajectory(t,x);
            end
            
            % System inputs.
            if isempty(obj.Input)
                u = zeros(obj.Dimensions.u,1);
            elseif isnumeric(obj.Input)
                u = obj.Input;
            else
                u = obj.Input(t,x);
            end

            % Evaluate polytopic system
            sys = obj.LPVSystem.Evaluate('Parameter',rho);

            % State and input perturbation
            x0 = sys.UserData.x0;
            u0 = sys.UserData.u0;
            dx = x - x0;
            du = u - u0;

            % LPV system output
            y = sys.C*dx + sys.D*du;

        end % Output
    end
end