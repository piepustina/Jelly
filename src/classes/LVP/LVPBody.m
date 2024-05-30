%#codegen
classdef LVPBody < Body
    %LVPBODY Class that models a locally volume preserving body
    
    properties (Access = public)
        Nodes                       (3, :) double
        Elements                    (10, :) double
        Primitives                  (:, 1) cell
        Centroids                   (3, :) double
        VertexElementsPositions     (3, 4, :) double
        ElementsVolume              (1, :) double
        Backbone                    (1, 1)
        ConnectorPoint              (3, 1) double% Point where the successor body is connected
        n = 0 % Number of DOFs
        MassDensity                 (1, 1) double
        YoungModulus                (1, 1) double
        PoissonRatio                (1, 1) double
        DampingFactor               (1, 1) double
    end

    properties
        Parameters = 0
    end

    properties (Access = private) % Make this private
        NNodes                          (1, 1) = 0
        NElements                       (1, 1) = 0
        NPrimitives                     (1, 1) = 0
        QIdx                            (:, 2) % Store the index of the start and end configuration for each primitive
        PrimitiveUsesBackbone           (:, 1) logical
        rCentroids                      (3, :)
        drCentroids                     (3, :)
        ddrCentroids                    (3, :)
        JrCentroids                     (3, :, :) % Jacobian of the centroids w.r.t. q
        CentroidMass                    (1, :)
        CentroidMassTensor              (1, 1, :) % Same as CentroidMass but in tensorized form for some computations
        DeformationGradientCentroid     (3, 3, :) % Deformation gradient at the centroids
        JDeformationGradientCentroid    (9, :, :) % Jacobian of the deformation gradient at the centroids w.r.t. q
        JCentroid                       (3, :, :) % Jacobian of the function f w.r.t. q
        % dDeformationGradientCentroid    (3, 3, :) % Time derivative of the deformation gradient at the centroids
        nGaussPoints                    (1, 1)
    end

    properties (Access = public, Constant)
        %Maximum number of allowed primitives for code generation support.
        MaxPrimitivesNumber                 = 10;
    end

    methods (Static)
        % Overload the loadobj method from the Body class
        function obj = loadobj(S)
            % Load the primitives
            Primitives  = cell(LVPBody.MaxPrimitivesNumber, 1);

            NP          = length(S.BodyPrimitives);
            for i=1:LVPBody.MaxPrimitivesNumber
                if i <= NP
                    loadobjPrimitive     = str2func(string(deblank(S.BodyPrimitivesClass(i, :))) + ".loadobj");
                    Primitives{i}        = loadobjPrimitive(S.BodyPrimitives(i));
                else
                    Primitives{i} = 0;
                end
            end

            % Build the LBPBody object
            obj = LVPBody(S.Nodes, S.Elements, Primitives, S.NGaussPoints, S.MassDensity, S.YoungModulus, S.PoissonRatio, S.DampingFactor);
        end

    end

    methods
        % Overload the saveobj method from the Body class
        function S = saveobj(obj)
            % NOTE: Since MATLAB does not support nested struct with cell arrays, for the primitives we have to use struct arrays all with fileds that are primitive datatypes.
            % Create a cell array containig a struct representation of the primitives            
            PrimitivesStuctArray = cellfun(@saveobj, obj.Primitives, "UniformOutput", false);
            PrimitivesClass      = cellfun(@class, obj.Primitives, "UniformOutput", false);

            % Build the struct representing the LVPBody
            S = struct('BodyType', class(obj), ...
                       'BodyDoF', obj.n, ...
                       'Nodes', obj.Nodes, ...
                       'Elements', obj.Elements, ...
                       'NGaussPoints', obj.nGaussPoints, ...
                       'MassDensity', obj.MassDensity, ...
                       'YoungModulus', obj.YoungModulus, ...
                       'PoissonRatio', obj.PoissonRatio, ...
                       'DampingFactor', obj.DampingFactor, ...
                       'BodyPrimitivesClass', char(PrimitivesClass), ...
                       'BodyPrimitives', {cell2mat(PrimitivesStuctArray)});
        end
    end

    methods
        function obj = LVPBody(Nodes, Elements, Primitives, NGaussPoints, MassDensity, YoungModulus, PoissonRatio, DampingFactor)
            %LVPBODY Construct a locally volume preserving body. Nodes and elements follow the same notation described in
            % - https://nl.mathworks.com/help/pde/ug/pde.femesh-properties.html#bup5xv9-Elements
            % - https://nl.mathworks.com/help/pde/ug/mesh-data.html
            arguments (Input)
                Nodes           (3, :)   double
                Elements        (10, :)  double % The elements are assumed to be tetrahedra
                Primitives      (:, 1)   cell
                NGaussPoints    (1, 1) = 10 %Number of Gaussian points for the computations along the backbone
                MassDensity     (1, 1) = 0
                YoungModulus    (1, 1) = 0
                PoissonRatio    (1, 1) = 0
                DampingFactor   (1, 1) = 0
            end
            % Assign the inputs
            obj.Nodes           = Nodes;
            obj.Elements        = Elements;
            obj.MassDensity     = MassDensity;
            obj.YoungModulus    = YoungModulus;
            obj.PoissonRatio    = PoissonRatio;
            obj.DampingFactor   = DampingFactor;
            obj.nGaussPoints    = NGaussPoints;
            % Compute the number of nodes, elements, primitives and centroids
            obj.NNodes      = size(obj.Nodes, 2);
            obj.NElements   = size(obj.Elements, 2);
            
            % Store the primitives and their number   
            obj.Primitives  = Primitives;
            obj.PrimitiveUsesBackbone = false(length(Primitives), 1);
            obj.NPrimitives = 0;
            obj.QIdx        = zeros(LVPBody.MaxPrimitivesNumber, 2);
            obj.n           = 0;
            for i = 1:length(Primitives)
                if isa(Primitives{i}, "LVPPrimitive")
                    obj.NPrimitives   = obj.NPrimitives + 1;
                    % Compute the start and end index for the current primitive
                    obj.QIdx(i, 1:2)  = [obj.n + 1, obj.n + Primitives{i}.n];
                    % Update the total number of DOFs
                    obj.n             = obj.n + Primitives{i}.n;
                    % Check if the current primitive uses the backbone
                    if isa(Primitives{i}, "BackbonePrimitive")
                        obj.PrimitiveUsesBackbone(i) = true;
                    end
                end
            end
            
            % Compute the centroids and the position of the elements
            obj.Centroids                   = zeros(3, obj.NElements);
            obj.VertexElementsPositions     = zeros(3, 4, obj.NElements);
            obj.ComputeCentroids();
            % Preallocate the data structures that will store the status of the centroids w.r.t. the center of mass
            obj.rCentroids                      = obj.Centroids;
            obj.drCentroids                     = zeros(3, obj.NElements);
            obj.ddrCentroids                    = zeros(3, obj.NElements);
            obj.JrCentroids                     = zeros(3, obj.n, obj.NElements);
            obj.DeformationGradientCentroid     = zeros(3, 3, obj.NElements);
            obj.JCentroid                       = zeros(3, obj.n, obj.NElements);
            % obj.dDeformationGradientCentroid    = zeros(3, 3, obj.NElements);
            obj.JDeformationGradientCentroid    = zeros(9, obj.n, obj.NElements);

            % Compute the volume of each element.
            % Since volume is locally preserved, this function is called only once for all.
            obj.ElementsVolume = zeros(1, obj.NElements);
            obj.ComputeElementsVolume();
            obj.CentroidMass        = obj.MassDensity*obj.ElementsVolume;
            obj.CentroidMassTensor  = reshape(obj.CentroidMass, 1, 1, []);

            % Compute the position of the connector point for the following body.
            % It is assumed that this point lies along the body backbone
            obj.ConnectorPoint = [0; ...
                                  0; ...
                                  max(obj.Nodes(3, 1:obj.NNodes))];

            % Construct the backbone
            obj.Backbone = LVPBackbone(max(obj.Nodes(3, 1:obj.NNodes)), NGaussPoints, Primitives, obj.Centroids(3, :));
            
            % Allocate the variables inherited from the Body class whose dimension depends on the number of DOFs
            obj.v_par_              = zeros(3, obj.n);
            obj.omega_par_          = zeros(3, obj.n);
            obj.int_dr_X_pv_r_      = zeros(obj.n, 3);
            obj.int_pv_r_O_dd_r_    = zeros(obj.n, 1);
            obj.grad_int_dr_        = zeros(obj.n, 3);
            obj.grad_int_r_X_dr_    = zeros(obj.n, 3);
            obj.grad_J_             = zeros(3, 3, obj.n);
            obj.grad_v_com_         = zeros(obj.n, 3);
            obj.K_                  = zeros(obj.n, 1);
            obj.D_                  = zeros(obj.n, 1);
        end
        
        % Getters methods
        function C = get.Centroids(obj)
            C = obj.Centroids;
        end

        function N = get.Nodes(obj)
            N = obj.Nodes;
        end

        function E = get.Elements(obj)
            E = obj.Elements;
        end

        function B = get.Backbone(obj)
            B = obj.Backbone;
        end

        % Override the Body method Update to update the status of the LVPBody
        function Update(obj, q, dq, ddq, options)
            arguments
                obj (1, 1) LVPBody
                q   (:, 1) {mustBeVector}
                dq  (:, 1) {mustBeVector}
                ddq (:, 1) {mustBeVector}
                options.EvaluateKinematicTerms (1, 1) logical = true
                options.EvaluateInertialTerms (1, 1) logical  = true
                options.EvaluateExternalForces (1, 1) logical = true
            end

            % Update the kinematics
            if options.EvaluateKinematicTerms == true
                obj.UpdateKinematics(q, dq, ddq);
            end
            
            % Update the interital terms
            % NOTE: the configuration variables are not required since all the computations are done using "state" variables of the class
            if options.EvaluateInertialTerms == true
                obj.p_com_              = obj.p_com();
                obj.v_com_rel_          = obj.v_com_rel();
                obj.a_com_rel_          = obj.a_com_rel();
                obj.I_                  = obj.I();
                obj.m_                  = obj.m();
                obj.J_                  = obj.J();
                obj.int_dr_             = obj.int_dr(q, dq);%TODO: To be removed since it is not required
                obj.int_ddr_            = obj.int_ddr(q, dq, ddq);%TODO: To be removed since it is always zero
                obj.int_r_X_dr_         = obj.int_r_X_dr();
                obj.int_r_X_ddr_        = obj.int_r_X_ddr();
                obj.int_dr_X_pv_r_      = obj.int_dr_X_pv_r();
                obj.int_pv_r_O_dd_r_    = obj.int_pv_r_O_dd_r();
                obj.int_dr_O_dr_        = obj.int_dr_O_dr(q, dq);%TODO: To be removed since it is always zero
                obj.grad_int_dr_        = obj.grad_int_dr(q);%TODO: To be removed since it is always zero
                obj.grad_int_r_X_dr_    = obj.grad_int_r_X_dr();
                obj.grad_J_             = obj.grad_J();
                obj.grad_v_com_         = obj.grad_v_com();
            end

            % Evaluate the external generalized forces
            if options.EvaluateExternalForces == true
                obj.K_                  = obj.K(q, "UpdateKinematics", ~options.EvaluateKinematicTerms);
                obj.D_                  = obj.D(q, dq, "UpdateKinematics", ~options.EvaluateKinematicTerms);
            end
        end
        
        %% Kinematics terms
        % Transformation from base to body tip
        function T_ = T(obj, q)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                T_ (4, 4)
            end
            T_ = obj.Backbone.T(q);
        end
        
        % Relative velocity of the body tip
        function v_rel_ = v_rel(obj, q, dq)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
                dq  (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                v_rel_ (3, 1)
            end
            v_rel_ = obj.Backbone.v_rel(q, dq);
        end
        % Relative angular velocity
        function omega_rel_ = omega_rel(obj, q, dq)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
                dq  (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                omega_rel_ (3, 1)
            end
            omega_rel_ = obj.Backbone.omega_rel(q, dq);
        end
        % Relative linear acceleration of the tip
        function a_rel_ = a_rel(obj, q, dq, ddq)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
                dq  (:, 1) = zeros(obj.n, 1)
                ddq (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                a_rel_ (3, 1)
            end
            a_rel_ = obj.Backbone.a_rel(q, dq, ddq);
        end
        % Relative angular acceleration
        function domega_rel_ = domega_rel(obj, q, dq, ddq)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
                dq  (:, 1) = zeros(obj.n, 1)
                ddq (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                domega_rel_ (3, 1)
            end
            domega_rel_ = obj.Backbone.domega_rel(q, dq, ddq);
        end
        % Jacobian of the linear velocity of the tip with respect to dq in the tip frame
        function v_par_ = v_par(obj, q)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
            end

            arguments (Output)
                v_par_ (3, :)
            end
            v_par_ = obj.Backbone.v_par(q);
        end
        % Jacobian of the angular velocity with respect to dq in the tip frame
        function omega_par_ = omega_par(obj, q)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
            end

            arguments (Output)
                omega_par_ (3, :)
            end
            omega_par_ = obj.Backbone.omega_par(q);
        end
        
        % Center of mass position
        function p_com_  = p_com(obj, q)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1)  = zeros(obj.n, 1)
            end
            arguments (Output)
                p_com_ (3, 1)
            end
            
            % Update the kinematics
            if nargin == 2
                obj.UpdateKinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            end
            % Return the CoM position comuted in Kinematics
            p_com_  = obj.p_com_;
        end
        % Relative velocity of the center of mass
        function v_com_rel_ = v_com_rel(obj, q, dq)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1)  = zeros(obj.n, 1)
                dq   (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                v_com_rel_ (3, 1)
            end
            
            % Update the kinematics
            if nargin == 2
                obj.UpdateKinematics(q, dq, zeros([obj.n, 1], 'like', q));
            end
            % Return the CoM velocity
            v_com_rel_   = obj.v_com_rel_;
        end
        % Relative acceleration of the center of mass
        function a_com_rel_  = a_com_rel(obj, q, dq, ddq)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
                dq  (:, 1) = zeros(obj.n, 1)
                ddq (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                a_com_rel_ (3, 1)
            end
            % Update the kinematics
            if nargin == 2
                obj.UpdateKinematics(q, dq, ddq);
            end
            % Compute the acceleration
            a_com_rel_ = obj.a_com_rel_;
        end

        %% Dynamic terms
        
        % Inertia matrix
        function I_ = I(obj, q)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
            end

            arguments (Output)
                I_ (3, 3)
            end
            
            % Update the DK if q was provided
            if nargin == 2
                obj.UpdateKinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            end

            % Compute the inertia
            Sr_     = skew(obj.rCentroids);
            Sr_T_   = pagetranspose(Sr_);
            I_      = sum(obj.CentroidMassTensor.*pagemtimes(Sr_, Sr_T_), 3);
        end

        % Mass
        function m_ = m(obj)
            arguments (Input)
                obj (1, 1) LVPBody
            end
            arguments (Output)
                m_ (1, 1)
            end
            m_ = sum(obj.CentroidMass);
        end

        % Time derivative of the inertia matrix
        function J_ = J(obj, q, dq)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
               dq   (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                J_ (3, 3)
            end
            
            % Update the kinematics
            if nargin >= 2
                obj.UpdateKinematics(q, dq, zeros([obj.n, 1], 'like', q));
            end

            %Compute the time derivative of the inertia
            Sr_     = skew(obj.rCentroids);
            Sr_T_   = pagetranspose(Sr_);
            dSr_    = skew(obj.drCentroids);
            dSr_T_  = pagetranspose(dSr_);
            J_      = sum(obj.CentroidMassTensor.*(pagemtimes(dSr_T_, Sr_) + pagemtimes(Sr_T_, dSr_)), 3);
        end
        

        % Integral of \dot{r}
        % TODO: This is always zero, remove
        function int_dr_ = int_dr(obj, q, dq)
            int_dr_ = zeros(3, 1, 'like', q);
        end

        % Integral of \ddot{r}
        % TODO: This is always zero, remove
        function int_ddr_ = int_ddr(obj, q, dq, ddq)
            int_ddr_ = zeros(3, 1, 'like', q);
        end

        % Integral of \cross(r, \dot{r})
        function int_r_X_dr_ = int_r_X_dr(obj, q, dq)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
               dq   (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                int_r_X_dr_ (3, 1)
            end

            % Update the kinematics
            if nargin >= 2
                obj.UpdateKinematics(q, dq, zeros([obj.n, 1], 'like', q));
            end
            
            % Compute the term
            int_r_X_dr_ = sum(obj.CentroidMassTensor.*pagemtimes(skew(obj.rCentroids), reshape(obj.drCentroids, 3, 1, [])), 3);
       end
        
        % Integral of \cross(r, \ddot{r})
        function int_r_X_ddr_ = int_r_X_ddr(obj, q, dq, ddq)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
               dq   (:, 1) = zeros(obj.n, 1)
              ddq   (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                int_r_X_ddr_ (3, 1)
            end
            % Update the kinematics
            if nargin >= 2            
                obj.UpdateKinematics(q, dq, ddq);
            end

            % Compute the term
            int_r_X_ddr_    = sum(obj.CentroidMassTensor.*pagemtimes(skew(obj.rCentroids), reshape(obj.ddrCentroids, 3, 1, [])), 3);
        end
        
        
        % Integral of \cross(\dor{r}, \jacobian{r}{q})
        function int_dr_X_pv_r_ = int_dr_X_pv_r(obj, q, dq)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
               dq   (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                int_dr_X_pv_r_ (:, 3)
            end
            if nargin >= 2 % Update the kinematics
                obj.UpdateKinematics(q, dq, zeros([obj.n, 1], 'like', q));
            end

            % Compute the term
            int_dr_X_pv_r_ = sum(pagetranspose(obj.CentroidMassTensor.*pagemtimes(skew(obj.drCentroids), obj.JrCentroids)), 3);
        end


        %Integral of \dot(\jacobian{r}{q}, \ddot{r})
        function int_pv_r_O_dd_r_ = int_pv_r_O_dd_r(obj, q, dq, ddq)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
               dq   (:, 1) = zeros(obj.n, 1)
              ddq   (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                int_pv_r_O_dd_r_ (:, 1)
            end

            % Update the kinematics
            if nargin >= 2
                obj.UpdateKinematics(q, dq, ddq);
            end
            % Compute the term
            int_pv_r_O_dd_r_ = sum(obj.CentroidMassTensor.*pagemtimes(pagetranspose(obj.JrCentroids), reshape(obj.ddrCentroids, 3, 1, [])), 3);
        end

        % Integral of \dot(\dot{r}, \dot{r})
        % TODO: This is always zero, remove
        function int_dr_O_dr_ = int_dr_O_dr(obj, q, dq)
            int_dr_O_dr_ = zeros(1, 1, 'like', q);
        end
        
        % Jacobian of the integral of \dot{r}
        % TODO: This is always zero, remove
        function grad_int_dr_ = grad_int_dr(obj, q)
            grad_int_dr_ = zeros(obj.n, 3, 'like', q);
        end

        % Jacobian of the integral of \cross{r, \dot{r}}
        function grad_int_r_X_dr_ = grad_int_r_X_dr(obj, q)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                grad_int_r_X_dr_ (:, 3)
            end
            
            % Update the DK
            if nargin == 2
                obj.UpdateKinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            end

            % Compute the term
            grad_int_r_X_dr_ = sum(obj.CentroidMassTensor.*pagetranspose(pagemtimes(skew(obj.rCentroids), obj.JrCentroids)), 3);
        end
        
        %Jacobian of the time derivative of the inertia
        function grad_J_ = grad_J(obj, q)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                grad_J_ (3, 3, :)
            end
            
            % Update the DK
            if nargin == 2
                obj.UpdateKinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            end
            % Compute the term
            grad_J_ = zeros(3, 3, obj.n, 'like', obj.rCentroids(1:3, 1));
            % Store useful variables
            Sr      = skew(obj.rCentroids);
            Sr_T    = pagetranspose(Sr);
            % Iterate over all DoFs
            for j = 1:obj.n
                % Get the required terms
                SJr     = skew(squeeze(obj.JrCentroids(1:3, j, :)));
                SJr_T   = pagetranspose(SJr);
                % Compute the term for the j-th configuration variable of the body
                grad_J_(1:3, 1:3, j) = sum(obj.CentroidMassTensor.*(pagemtimes(SJr_T, Sr) + pagemtimes(Sr_T, SJr)), 3);
            end
        end
        
        % Jacobian of the center of mass velocity
        function grad_v_com_ = grad_v_com(obj, q)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                grad_v_com_ (:, 3)
            end
            
            %Update the DK
            if nargin == 2 
                obj.UpdateKinematics(q, zeros([obj.n, 1], 'like', q), zeros([obj.n, 1], 'like', q));
            end
            % Return the output compute during the update of the kinematics
            grad_v_com_ = obj.grad_v_com_;
        end

        % Compute the deformation gradient
        function F = DeformationGradient(obj, q)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                F (3, 3, :)
            end
            
            % Update the kinematics
            if nargin == 2
                obj.UpdateKinematics(q, zeros(obj.n, 1, "like", q), zeros(obj.n, 1, "like", q));
            end

            F = obj.DeformationGradientCentroid;
        end

        function Jq = BodyJacobian(obj, q)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
            end
            arguments (Output)
                Jq (3, :, :)
            end
            
            % Update the kinematics
            if nargin == 2
                obj.UpdateKinematics(q, zeros(obj.n, 1, "like", q), zeros(obj.n, 1, "like", q));
            end

            Jq = obj.JCentroid;
        end
        
        % Compute the stiffness force
        function Kq = K(obj, q, options)
            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1) = zeros(obj.n, 1)
                options.UpdateKinematics (1, 1) = false
            end

            arguments (Output)
                Kq  (:, 1)
            end

            % Update the Jacobians through the kinematics
            if options.UpdateKinematics == true
                obj.UpdateKinematics(q, zeros(obj.n, 1, "like", q), zeros(obj.n, 1, "like", q));
            end
            
            % Get the deformation gradient
            F       = repmat(reshape(obj.DeformationGradientCentroid, 3, 3, 1,obj.NElements), 1, 1, obj.n, 1);
            % Get the Jacobian of the deformation gradient w.r.t. q
            JF      = reshape(obj.JDeformationGradientCentroid, 3, 3, obj.n, obj.NElements);

            % Compute the projection of the right Caucy strain tensor in the configuration space
            gradC   = pagemtimes(pagetranspose(JF), F) + pagemtimes(pagetranspose(F), JF);

            
            % Extract the diagonal elements, i.e., the principal stretches
            if obj.n ~= 1
                GradCStretches = squeeze(gradC(1, 1, 1:obj.n, 1:obj.NElements)) ...
                               + squeeze(gradC(2, 2, 1:obj.n, 1:obj.NElements)) ...
                               + squeeze(gradC(3, 3, 1:obj.n, 1:obj.NElements));
            else
                GradCStretches = reshape(squeeze(gradC(1, 1, 1:obj.n, 1:obj.NElements)) ...
                               + squeeze(gradC(2, 2, 1:obj.n, 1:obj.NElements)) ...
                               + squeeze(gradC(3, 3, 1:obj.n, 1:obj.NElements)), 1, []);
            end

            % Compute the Lamè constants
            %lambda = obj.YoungModulus*obj.PoissonRatio/((1+obj.PoissonRatio)*(1-2*obj.PoissonRatio));
            mi     = obj.YoungModulus/(2*(1+obj.PoissonRatio));

            % Compute the elastic energy
            Kq          = zeros(obj.n, 1, "like", q);
            % if norm(q) <= 1e-4
            %     return;
            % end
            Kq          = mi*sum(GradCStretches.*obj.ElementsVolume, 2);

        end

        % Compute the damping force
        function Dq = D(obj, q, dq, options)

            arguments (Input)
                obj (1, 1) LVPBody
                q   (:, 1)                  = zeros(obj.n, 1)
                dq  (:, 1)                  = zeros(obj.n, 1)
                options.UpdateKinematics    = false
            end

            arguments (Output)
                Dq  (:, 1)
            end

            % Update the Jacobians through the kinematics
            if options.UpdateKinematics == true
                obj.UpdateKinematics(q, dq, zeros(obj.n, 1, "like", q));
            end

            %Dq = obj.DampingFactor*dq;return;

            % Get the time derivative of the deformation gradient
            dF      = repmat(reshape(pagemtimes(obj.JDeformationGradientCentroid, dq), 3, 3, 1, obj.NElements), 1, 1, obj.n, 1);

            % Get the Jacobian of the deformation gradient w.r.t. q
            JF      = reshape(obj.JDeformationGradientCentroid, 3, 3, obj.n, obj.NElements);
                        
            % Compute the projection of the right Caucy strain tensor in the configuration space
            dgradC  = pagemtimes(pagetranspose(JF), dF) + pagemtimes(pagetranspose(dF), JF);

            % Extract the diagonal elements, i.e., the time derivative of the principal stretches
            if obj.n ~= 1
                dGradCStretches = squeeze(dgradC(1, 1, 1:obj.n, 1:obj.NElements)) ...
                                + squeeze(dgradC(2, 2, 1:obj.n, 1:obj.NElements)) ...
                                + squeeze(dgradC(3, 3, 1:obj.n, 1:obj.NElements));
            else
                dGradCStretches = reshape(squeeze(dgradC(1, 1, 1:obj.n, 1:obj.NElements)) ...
                                + squeeze(dgradC(2, 2, 1:obj.n, 1:obj.NElements)) ...
                                + squeeze(dgradC(3, 3, 1:obj.n, 1:obj.NElements)), 1 , []);
            end
            % Compute the Lamè constants
            %lambda = obj.YoungModulus*obj.PoissonRatio/((1+obj.PoissonRatio)*(1-2*obj.PoissonRatio));
            mi     = obj.YoungModulus/(2*(1+obj.PoissonRatio));

            % Compute the elastic energy
            Dq          = zeros(obj.n, 1, "like", q);
            Dq          = obj.DampingFactor*mi*sum(dGradCStretches.*obj.ElementsVolume, 2);  
        end
        
        %% General utility functions
        % Volume computation the elements
        function ComputeElementsVolume(obj)
            for i = 1:obj.NElements
                d         = obj.VertexElementsPositions(1:3, 4*(i-1)+4);
                a_minus_d = obj.VertexElementsPositions(1:3, 4*(i-1)+1) - d;
                b_minus_d = obj.VertexElementsPositions(1:3, 4*(i-1)+2) - d;
                c_minus_d = obj.VertexElementsPositions(1:3, 4*(i-1)+3) - d;
                obj.ElementsVolume(i) = abs(a_minus_d'*cross(b_minus_d, c_minus_d))/6;
            end
        end
        
        %% Plot functions
        % Override the plot function
        function p = plot(obj, q, options)
            % Plot the body in the current configuration
            arguments
                obj                           (1, 1) LVPBody
                q                             (:, 1) double = zeros(obj.n, 1)
                options.BaseTransformation    (4, 4) double = eye(4)
                options.Color                 (1, 3) double = [0 160 219]./256;
                options.FaceAlpha             (1, 1) double = 1
                options.LineStyle             (1, 1) = "none";% Possible values: "-" | "--" | ":" | "-." | "none"
            end
            % Set the nodes as the target computation of the LVPBackbone
            obj.Backbone.setPoints(obj.Nodes(3, :));
            % Compute the current position of the body points as a function of q
            [Nodesq, ~, ~, ~, ~, ~] = obj.UpdateKinematicsAtX(q, zeros(size(q)), zeros(size(q)), obj.Nodes);
            % Get the pose of the backbone
            BackbonePose = pagemtimes(options.BaseTransformation, obj.Backbone.UniquePose);
            % Set back the centroids as target points for the computations of the LVPBackbone
            obj.Backbone.setPoints(obj.Centroids(3, :));
            % Apply the base transformation to all the bodies points
            TNodes = options.BaseTransformation*[Nodesq; ones(1, obj.NNodes)];
            TNodes = TNodes(1:3, :)';
            % Build the faces. The number of faces is 4*obj.NElements
            Faces  = zeros(4*obj.NElements, 4);
            for i = 1:obj.NElements
                Faces(4*(i-1)+1, :) = [obj.Elements(1, i), obj.Elements(2, i), obj.Elements(3, i), obj.Elements(1, i)];
                Faces(4*(i-1)+2, :) = [obj.Elements(1, i), obj.Elements(2, i), obj.Elements(4, i), obj.Elements(1, i)];
                Faces(4*(i-1)+3, :) = [obj.Elements(1, i), obj.Elements(3, i), obj.Elements(4, i), obj.Elements(1, i)];
                Faces(4*(i-1)+4, :) = [obj.Elements(2, i), obj.Elements(3, i), obj.Elements(4, i), obj.Elements(2, i)];
            end
            % Plot the points using the patch object
            p = patch('Faces', Faces, 'Vertices', TNodes, ...
                        'FaceColor', options.Color, ...
                        'FaceAlpha', options.FaceAlpha, ...
                        'LineStyle', options.LineStyle);

            % Plot also the backbone
            %[~, idxZ] = sort(squeeze(BackbonePose(3, 4, :)));
            %plot3(squeeze(BackbonePose(1, 4, idxZ)), squeeze(BackbonePose(2, 4, idxZ)), squeeze(BackbonePose(3, 4, idxZ)), "LineWidth", 2, "Color", "r");
            plot3(squeeze(BackbonePose(1, 4, :)), squeeze(BackbonePose(2, 4, :)), squeeze(BackbonePose(3, 4, :)), "LineWidth", 2, "Color", "r");
        end
        
        % Evaluate the kinematics in the current configuration
        function [xq, dxq, ddxq, Jq, Jx_ref, JJq] = UpdateKinematicsAtX(obj, q, dq, ddq, x)
            arguments (Input)
                obj         (1, 1) LVPBody
                q           (:, 1) = zeros(obj.n, 1)
                dq          (:, 1) = zeros(obj.n, 1)
                ddq         (:, 1) = zeros(obj.n, 1)
                x           (3, :) = zeros(3, 1)
            end
            
            arguments (Output)
                xq      (3, :)
                dxq     (3, :)
                ddxq    (3, :)
                Jq      (3, :, :)
                Jx_ref  (3, 3, :)
                % Jdx_ref (3, 3, :)
                JJq     (9, :, :)
            end

            % Preallocate the output
            Nx          = size(x, 2);
            xq          = x;
            dxq         = zeros(3, Nx, "like", x);
            ddxq        = zeros(3, Nx, "like", x);
            Jq          = zeros(3, obj.n, Nx, "like", x);
            Jx_ref      = repmat(eye(3, "like", x), 1, 1, Nx);
            % Jdx_ref     = repmat(eye(3, "like", x), 1, 1, Nx);
            JJq         = zeros(9, obj.n, Nx, "like", x);% Jacobian of the primitive Jacobian w.r.t. x_ref w.r.t. q
            I3          = repmat(eye(3, "like", x), 1, 1, Nx);

            % Create a configuration vector for the backbone that is progressively update as the primitives are applied.
            % This way the backbone and the overall body keep a consistent status over the computations.
            qProgressive    = zeros(obj.n, 1, "like", q);
            dqProgressive   = zeros(obj.n, 1, "like", q);
            ddqProgressive  = zeros(obj.n, 1, "like", q);
            
            % Update the kinematics by applying progressively the primitives to the body
            for i = 1:obj.MaxPrimitivesNumber
                if i <= obj.NPrimitives
                    % Check that we are dealing with a primitive (required for code generation support)
                    if ~isa(obj.Primitives{i}, "LVPPrimitive")
                        continue;
                    end
                    % Extract the configuration associated with the primitive
                    qIdx            = obj.QIdx(i, 1):obj.QIdx(i, 2);
                    qPrimitive      = q(qIdx, 1);
                    dqPrimitive     = dq(qIdx, 1);
                    ddqPrimitive    = ddq(qIdx, 1);
                    
                    % % If the configuration variables and time derivatives are zero for the current primitive, we do not have to evaluate the primitive.
                    % if all(qPrimitive == 0) && all(dqPrimitive == 0) && all(ddqPrimitive == 0)
                    %     continue;
                    % end

                    % qProgressive    = zeros(obj.n, 1, "like", q);
                    % dqProgressive   = zeros(obj.n, 1, "like", q);
                    % ddqProgressive  = zeros(obj.n, 1, "like", q);
    
                    % Update the configuration to account for the current deformation primitive
                    qProgressive(qIdx, 1)      = qPrimitive;
                    dqProgressive(qIdx, 1)     = dqPrimitive;
                    ddqProgressive(qIdx, 1)    = ddqPrimitive;
    
                    % Update the status of the backbone if the primitive uses the backbone
                    if obj.PrimitiveUsesBackbone(i) == true
                        obj.Backbone.Update(qProgressive, dqProgressive, ddqProgressive);
                    end
                    
                    % Apply the primitive and its time derivatives
                    % Note that the order of update of accelerations, velocities and positions is relevant because we use single variables
                    % [xq, dxq, ddxq, JfxqPrimitive_vec, Jfxx_vec, Jfxx_ref_vec, Jdfxx_vec, Jdfxx_ref_vec, Jdfxdx_vec, JJfx_q_vec, JJfx_ref_x_vec, JJfx_ref_q_vec] = obj.Primitives{i}.Update(obj.Backbone, qPrimitive, dqPrimitive, ddqPrimitive, xq, dxq, ddxq);

                    [xq, dxq, ddxq, JfxqPrimitive_vec, Jfxx_vec, Jfxx_ref_vec, JJfx_q_vec, JJfx_ref_x_vec, JJfx_ref_q_vec] = obj.Primitives{i}.Update(obj.Backbone, qPrimitive, dqPrimitive, ddqPrimitive, xq, dxq, ddxq);

                    
                    % Convert the Jacobians to 3D tensors
                    Jfxq            = reshape(JfxqPrimitive_vec, 3, obj.n, Nx); 
                    Jfxx            = reshape(Jfxx_vec, 3, 3, Nx);
                    Jfxx_ref        = reshape(Jfxx_ref_vec, 3, 3, Nx);
                    % Jdfxx           = reshape(Jdfxx_vec, 3, 3, Nx);
                    % Jdfxx_ref       = reshape(Jdfxx_ref_vec, 3, 3, Nx);
                    % Jdfxdx          = reshape(Jdfxdx_vec, 3, 3, Nx);
                    JJfx_q          = reshape(JJfx_q_vec, 9, obj.n, Nx);
                    JJfx_ref_x      = reshape(JJfx_ref_x_vec, 9, 3, Nx);
                    JJfx_ref_q      = reshape(JJfx_ref_q_vec, 9, obj.n, Nx) + pagemtimes(JJfx_ref_x, Jq);
                    % Update the Jacobian w.r.t. q
                    if i == 1
                        JJq     = JJfx_ref_q + JJfx_q;
                        Jq      = Jfxq;
                        Jx_ref  = Jfxx_ref + Jfxx;
                        % Jdx_ref = Jdfxx_ref + Jdfxx;
                    else
                        JJq     = JJfx_ref_q + ...
                                + pagemtimes(pagemkron(I3, Jfxx), JJq) ...
                                + pagemtimes(pagemkron(pagetranspose(Jx_ref), I3), JJfx_q);
                        % Jdx_ref = Jdfxx_ref + pagemtimes(Jdfxx, Jx_ref) + pagemtimes(Jdfxdx, Jdx_ref);
                        Jq      = Jfxq + pagemtimes(Jfxx, Jq);
                        Jx_ref  = Jfxx_ref + pagemtimes(Jfxx, Jx_ref);
                    end
                end
            end
        end
    
        function [xq, dxq, ddxq, Jq, Jx_ref, JJq] = UpdateKinematics(obj, q, dq, ddq)
            arguments (Input)
                obj         (1, 1) LVPBody
                q           (:, 1) = zeros(obj.n, 1)
                dq          (:, 1) = zeros(obj.n, 1)
                ddq         (:, 1) = zeros(obj.n, 1)
            end
            % Update the status of the centroids
            % [xq, dxq, ddxq, Jq, Jx_ref, Jdx_ref, JJq] = obj.UpdateKinematicsAtX(q, dq, ddq, obj.Centroids);
            [xq, dxq, ddxq, Jq, Jx_ref, JJq] = obj.UpdateKinematicsAtX(q, dq, ddq, obj.Centroids);

            % Store the deformation gradient
            obj.DeformationGradientCentroid     = Jx_ref;
            % obj.dDeformationGradientCentroid    = pagemtimes(JJq, dq);
            obj.JDeformationGradientCentroid    = JJq;
            % Store the gradient of f(x, q) w.r.t. q
            obj.JCentroid                       = Jq;
            
            % We assume that the successor body is attached to center of the body, i.e., at the end of the backbone
            obj.T_                  = obj.Backbone.T_;
            obj.v_rel_              = obj.Backbone.v_rel_;
            obj.omega_rel_          = obj.Backbone.omega_rel_;
            obj.a_rel_              = obj.Backbone.a_rel_;
            obj.domega_rel_         = obj.Backbone.domega_rel_;
            obj.v_par_              = obj.Backbone.v_par_;
            obj.omega_par_          = obj.Backbone.omega_par_;

            % Compute the position of the body elements in the body (tip) frame
            RT                      = obj.T_(1:3, 1:3)';
            d                       = obj.T_(1:3, 4);
            dd                      = obj.Backbone.EtaGauss(4:6, end);
            ddd                     = obj.Backbone.dEtaGauss(4:6, end);
            Jd                      = obj.Backbone.JEtaGauss(4:6, 1:obj.n, end);
            omega                   = obj.Backbone.EtaGauss(1:3, end);
            Jomega                  = obj.Backbone.JEtaGauss(1:3, 1:obj.n, end);
            domega                  = obj.Backbone.dEtaGauss(1:3, end);
            dRT                     = -RT*skew(omega);
            ddRT                    = -RT*skew(domega) -dRT*skew(omega);
            
            % Compute the relative position of the centroids in the body frame
            pElements               = RT*(xq - d);
            dpElements              = dRT*(xq - d) + RT*(dxq  - dd);
            JdpElements             = pagemtimes(RT, Jq - Jd + pagemtimes(skew(xq-d), Jomega));
            ddpElements             = ddRT*(xq -d) + 2*dRT*(dxq  - dd) + RT*(ddxq  - ddd);

            % Compute the position, velocity and acceleration of the center of mass
            m                       = obj.m();% Mass of the entire body
            dm                      = obj.CentroidMass; % Mass of each infinitesimal centroid
            dm_tensor               = obj.CentroidMassTensor; % Tensorided version of dm for the computation of the Jacobians
            % Com position
            obj.p_com_              = zeros(3, 1, 'like', q);
            obj.p_com_              = (sum(dm.*pElements, 2))/m;
            % CoM velocity
            obj.v_com_rel_          = zeros(3, 1, 'like', q);
            obj.v_com_rel_          = (sum(dm.*dpElements, 2))/m;
            % CoM acceleration
            obj.a_com_rel_          = zeros(3, 1, 'like', q);
            obj.a_com_rel_          = (sum(dm.*ddpElements, 2))/m;
            % Position Jacobian of the CoM
            J_p_comi                = zeros(3, obj.n, 'like', q);
            J_p_comi(1:3, 1:obj.n)  = (sum(pagemtimes(dm_tensor, JdpElements), 3))/m;
            
            % Compute the relative position, first, second order time derivative and Jacobian of the body points in the body frame and w.r.t the Center of Mass
            obj.rCentroids          = pElements     - obj.p_com_;
            obj.JrCentroids         = JdpElements   - J_p_comi;
            obj.drCentroids         = dpElements    - obj.v_com_rel_;
            obj.ddrCentroids        = ddpElements   - obj.a_com_rel_;

            % Update the gradient of time derivative of the CoM position in the local frame (not the overall CoM velocity!)
            obj.grad_v_com_         = (RT*(skew(RT'*obj.p_com_)*Jomega - Jd + sum(pagemtimes(dm_tensor, Jq), 3)/m))';
        end
    
    end

    methods (Access = private)
        % Compute the positione of the centroid of the elements. The centroid position is used to approximate all the volume integrals.
        function ComputeCentroids(obj)
            % Get the position of the four vertices of the elements
            obj.VertexElementsPositions = reshape(obj.Nodes(1:3, obj.Elements(1:4, 1:obj.NElements)), 3, 4, obj.NElements);

            % Compute the centroids
            obj.Centroids = squeeze(mean(obj.VertexElementsPositions(1:3, 1:4, 1:obj.NElements), 2));
            
            % Sort the centroids by the z position and accordingly the tethrahedra
            [~, idxSort]  = sort(obj.Centroids(3, 1:obj.NElements));
            obj.Centroids(1:3, 1:obj.NElements)                     = obj.Centroids(1:3, idxSort);
            % Permute accordingly the elements and the position of their vertices
            obj.Elements(1:10, 1:obj.NElements)                     = obj.Elements(1:10, idxSort);
            obj.VertexElementsPositions(1:3, 1:4, 1:obj.NElements)  = obj.VertexElementsPositions(1:3, 1:4, idxSort);
        end
    end
end

