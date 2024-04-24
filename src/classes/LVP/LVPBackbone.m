classdef LVPBackbone < Body
    %Class that models the backbone of a Locally Volume Preserving (LVP) body.

    properties (Access = public)
        % Pose of the backbone at the backbone points backbone
        Pose        (4, 4, :) double
        % Velocity twist of the backbone at the backbone points
        dPose       (6, :) double
        % Acceleration twist of the backbone at the backbone points
        ddPose      (6, :) double
        % Jacobian of the pose
        JPose       (6, :, :) double
        % Strain at the points along the backbone
        Strain      (6, :) double
        % First order time derivative of the strain along the backbone points
        dStrain     (6, :) double
        % Second order time derivative of the strain along the backbone points
        ddStrain    (6, :) double
        % Jacobian of the strain at the points along the backbone
        JStrain     (6, :, :) double
        % Derivative of the strain w.r.t. s (curvilinear abscissa)
        Strain_ds   (6, :) double
        % Derivative of the first order time derivative of the strain w.r.t. s (curvilinear abscissa)
        dStrain_ds  (6, :) double
        % Derivative w.r.t. s of the strain Jacobian
        dJStrain    (6, :, :) double
        % Elongation at the points along the backbone
        Elongation  (:, 1) double
        % First order time derivative of the elongation at the points along the backbone
        dElongation  (:, 1) double
        % Second order time derivative of the elongation at the points along the backbone
        ddElongation  (:, 1) double
        % Jacobian of the elongation
        JElongation  (:, :) double
        % Cell array of body primitives
        Primitives (:, 1) cell
    end

    properties (Access = protected)
        NPrimitives                 (1, 1) = 0
        PrimitiveModifiesBackbone   (:, 1) logical % Vector of length NPrimitives such that PrimitiveModifiesBackbone(i) = true if the primitive modifies the backbone
        QIdx                        (:, 2)% Store the starting and end indexes of the configuration vector for each primitive
        % Gauss points and weights to compute the integrals from [0, Points(i)] 
        GaussPointsBackbonePoints   (1, :)
        GaussWeightsBackbonePoints  (1, :)
        xiGauss                     (6, :)
        dxiGauss                    (6, :)
        ddxiGauss                   (6, :)
        JxiGauss                    (6, :, :) % Strain basis
        xiGauss_ds                  (6, :) % Derivative of the strain w.r.t. s
        dxiGauss_ds                 (6, :) % Derivative of the first order time derivative of the strain w.r.t. s
        dJxiGauss                   (6, :, :) % Derivative w.r.t. s of the Jacobian strain
        ElongationGauss             (1, :)
        JElongationGauss            (:, :)
        dElongationGauss            (1, :)
        ddElongationGauss           (1, :)
        JxiGaussZanna               (6, :, :, 2) % Strain basis evaluated at the Zanna quadrature points for the Gauss points
        dJxiGaussZanna              (6, :, :, 2) % Strain basis evaluated at the Zanna quadrature points for the Gauss points
        
        % Points along the backbone (curvilinear abscissa) where the quantities must be evaluated
        Points                              (1, :) double
        % Number of points
        nPoints                             (1, 1) = 0
    end

    properties (Access = public)
        % Rest length of the backbone
        RestLength                  (1, 1) double   = 0
        % Reference value of the backbone strain, i.e., strain in the stress-free configuration
        ReferenceStrain             (6, 1) double   = [0;0;0;0;0;1];
        % Number of Gaussian points provided by the user for the numerical computation
        NGaussPoints                (1, 1) = 0;
        % Gaussian points
        GaussPoints                 (1, :) double
        % Gaussian weights
        GaussWeights                (1, :) double
        % Number of internal Gaussian points. This is the number of Gaussian points provided by the user with, in addition, the rest length (RestLength)
        NGaussPointsInt             (1, 1) = 0;
        % Homogeneous transformation matrix from Gaussian points to base
        gGauss                      (4, 4, :) double
        % Velocity twist of the Gaussian points expressed in the base frame
        EtaGauss                    (6, :)    double
        % Acceleration twist of the Gaussian points expressed in the base frame
        dEtaGauss                   (6, :)    double
        % Geometric Jacobian of the Gaussian points expressed in the base frame
        JEtaGauss                   (6, :, :) double
        % Acceleration twist of the Gaussian points expressed in the local frame
        dEtaGaussLocal              (6, :)    double
        % Geometric Jacobian of the Gaussian points epxressed in the local frame
        JEtaGaussLocal              (6, :, :) double
    end

    properties (Access = private)
        % Vector storing the unique points among the points along the backbone provided by the user
        UniquePoints                        (1, :) double
        % Mapping from the unique points to the points along the backbone provided by the user
        UniquePointsToPoints                (1, :) double
        % Number of unique points along the the backbone
        nUniquePoints                       (1, 1) = 0
        % Mapping from the unique points along the backbone to the Gauss points
        UniquePointsToGaussPoints           (:, 2) double
        % Expanded version of UniquePointsToGaussPoints to compute faster some terms
        UniquePointsToGaussPointsExpanded   (:, 1) double
        % Pose of the backbone at the backbone points backbone
        UniquePose        (4, 4, :) double
        % Velocity twist of the backbone at the backbone points
        UniquedPose       (6, :) double
        % Acceleration twist of the backbone at the backbone points
        UniqueddPose      (6, :) double
        % Jacobian of the pose
        UniqueJPose       (6, :, :) double
        % Strain at the points along the backbone
        UniqueStrain      (6, :) double
        % First order time derivative of the strain along the backbone points
        UniquedStrain     (6, :) double
        % Second order time derivative of the strain along the backbone points
        UniqueddStrain    (6, :) double
        % Jacobian of the strain at the points along the backbone
        UniqueJStrain     (6, :, :) double
        % Derivative of the strain w.r.t. s (curvilinear abscissa)
        UniqueStrain_ds   (6, :) double
        % Derivative of the first order time derivative of the strain w.r.t. s (curvilinear abscissa)
        UniquedStrain_ds  (6, :) double
        % Derivative w.r.t. s of the strain Jacobian
        UniquedJStrain    (6, :, :) double
        % Elongation at the points along the backbone
        UniqueElongation  (:, 1) double
        % First order time derivative of the elongation at the points along the backbone
        UniquedElongation  (:, 1) double
        % Second order time derivative of the elongation at the points along the backbone
        UniqueddElongation  (:, 1) double
        % Jacobian of the elongation
        UniqueJElongation  (:, :) double
        % Jacobian of the strain at the Zanna points of the specified unique points
        UniqueJStrainZanna1 (6, :, :)
        UniqueJStrainZanna2 (6, :, :)
        % Derivative of the Jacobian of the strain w.r.t. s at the Zanna points of the specified unique points
        UniquedJStrainZanna (6, :, :, 2)
    end
    
    % Implement the abstract properties inherited from the Body class
    properties
        % By default a LVPBackbone has no DOFs
        n                          = 0
        % Parameters
        Parameters
    end
    
    methods
        function obj = LVPBackbone(RestLength, NGaussPoints, Primitives, Points)
            %Build an instance of a LVB backbone.
            %Args:
            %   n ([double, int32])     : Number of DOFs of the backbone
            %   RestLength ([double])   : Rest length of the backbone
            %   NGaussPoints            : Number of Gaussian points used for the numerical computation of the integrals
            
            arguments 
                RestLength          (1, 1)
                NGaussPoints        (1, 1)
                Primitives          (:, 1) cell
                Points              (1, :) = []
            end

            
            % Call the superclass constructor.
            obj = obj@Body();
            
            % Store the geometric parameters of the backbone
            obj.RestLength      = RestLength;
            obj.NGaussPoints    = NGaussPoints;
            obj.Parameters      = [obj.RestLength; obj.NGaussPoints];
            
            % Store the primitives and their number. These primitives are shared with those of the body
            obj.NPrimitives                 = 0;
            obj.n                           = 0;
            obj.PrimitiveModifiesBackbone   = false(LVPBody.MaxPrimitivesNumber, 1);
            obj.QIdx                        = zeros(LVPBody.MaxPrimitivesNumber, 2);
            qStart                          = 1;
            for i = 1: length(Primitives)
                if isa(Primitives{i}, "LVPPrimitive")
                    obj.NPrimitives   = obj.NPrimitives + 1;
                    obj.n             = obj.n + Primitives{i}.n;
                    % Check if the primitive acts on the backbone
                    if isa(Primitives{i}, "BackbonePrimitive")
                        obj.PrimitiveModifiesBackbone(i) = true;
                    end
                    % Assign the starting and end index for the configurationv vector
                    obj.QIdx(i, 1:2) = [qStart, qStart + Primitives{i}.n - 1];
                    qStart           = obj.QIdx(i, 2) + 1;
                end
            end
            % Save the primitives
            obj.Primitives                  = Primitives;

            % Initialize all the quantities associated to the Gaussian points
            [obj.GaussPoints, obj.GaussWeights] = lgwt(obj.NGaussPoints, 0, obj.RestLength);
            obj.GaussPoints                     = [obj.GaussPoints, obj.RestLength];
            obj.GaussWeights                    = [obj.GaussWeights, 0];
            obj.NGaussPointsInt                 = obj.NGaussPoints + 1;
            obj.gGauss                          = repmat(eye(4), [1, 1, obj.NGaussPointsInt]);
            obj.EtaGauss                        = zeros(6, obj.NGaussPointsInt);
            obj.dEtaGauss                       = zeros(6, obj.NGaussPointsInt);
            obj.JEtaGauss                       = zeros(6, obj.n, obj.NGaussPointsInt);
            obj.xiGauss                         = repmat(obj.ReferenceStrain, 1, obj.NGaussPointsInt);
            obj.dxiGauss                        = zeros(6, obj.NGaussPointsInt);
            obj.ddxiGauss                       = zeros(6, obj.NGaussPointsInt);
            obj.JxiGauss                        = zeros(6, obj.n, obj.NGaussPointsInt);
            obj.JxiGaussZanna                   = zeros(6, obj.n, obj.NGaussPointsInt, 2);
            obj.xiGauss_ds                      = zeros(6, obj.NGaussPointsInt);
            obj.dxiGauss_ds                     = zeros(6, obj.NGaussPointsInt);
            obj.dJxiGauss                       = zeros(6, obj.n, obj.NGaussPointsInt);
            obj.dJxiGaussZanna                  = zeros(6, obj.n, obj.NGaussPointsInt, 2);
            obj.JEtaGaussLocal                  = zeros(6, obj.n, obj.NGaussPointsInt);
            obj.dEtaGaussLocal                  = zeros(6, obj.NGaussPointsInt);

            % Pre-evaluate the strain basis and its derivative w.r.t.s at the Gaussian points
            S = 0;
            for i = 1:obj.NGaussPointsInt
                [obj.JxiGauss(1:6, 1:obj.n, i), obj.dJxiGauss(1:6, 1:obj.n, i)] = obj.StrainBasis(obj.GaussPoints(i));
                % Compute the strain basis at the zanna points
                if i == 1
                    h                             = obj.GaussPoints(i);
                else
                    h                             = obj.GaussPoints(i)-obj.GaussPoints(i-1);
                end
                zanna1                            = S + h/2 - sqrt(3)*h/6;
                [JxiGaussZanna1, dJxiGaussZanna1] = obj.StrainBasis(zanna1);
                zanna2                            = S + h/2 + sqrt(3)*h/6;
                [JxiGaussZanna2, dJxiGaussZanna2] = obj.StrainBasis(zanna2);
                obj.JxiGaussZanna(1:6, 1:obj.n, i, 1)  = JxiGaussZanna1;
                obj.JxiGaussZanna(1:6, 1:obj.n, i, 2)  = JxiGaussZanna2;
                obj.dJxiGaussZanna(1:6, 1:obj.n, i, 1) = dJxiGaussZanna1;
                obj.dJxiGaussZanna(1:6, 1:obj.n, i, 2) = dJxiGaussZanna2;
                % Prepare for the next iteration
                S = obj.GaussPoints(i);
            end

            obj.ElongationGauss                 = zeros(1, obj.NGaussPointsInt);
            obj.JElongationGauss                = zeros(obj.n, obj.NGaussPointsInt);
            obj.dElongationGauss                = zeros(1, obj.NGaussPointsInt);
            obj.ddElongationGauss               = zeros(1, obj.NGaussPointsInt);
            

            % Assign the points along the backbone and other quantities associated to the points
            obj.setPoints(Points);
            
            % Initialize the Gaussian points and weights for the intervals [0, GaussPoints(i)].
            obj.GaussWeightsBackbonePoints  = zeros(1, obj.NGaussPointsInt);
            obj.GaussPointsBackbonePoints   = zeros(1, obj.NGaussPointsInt);
            xPrev = 0;
            for i = 1:obj.NGaussPointsInt
                [xGauss, wGauss] = lgwt(1, xPrev, obj.GaussPoints(i));
                obj.GaussPointsBackbonePoints(i)  = xGauss';
                obj.GaussWeightsBackbonePoints(i) = wGauss';
                xPrev = obj.GaussPoints(i);
            end

            % Preallocate the Body variables that depend on n and are used in by the LVPBody
            obj.v_par_                          = zeros(3, obj.n);
            obj.omega_par_                      = zeros(3, obj.n);
            
        end

        % Compute the compact mapping from the points to the Gaussian points
        function M = MapPointsToGaussianPoints(obj, Points)
            arguments (Input)
                obj     (1, 1) LVPBackbone
                Points  (1, :) double
            end

            arguments (Output)
                M       (:, 2)
            end
            % Gauss points along the backbone
            BackboneGaussPoints     = obj.GaussPoints;
            
            % Preallocate variables with repetition
            PointToBackboneExtended = zeros(1, length(Points));
            for i = 0:length(BackboneGaussPoints)-1
                
                % Get the current interval bounds
                if i == 0
                    IntvLb      = 0;
                else
                    IntvLb      = BackboneGaussPoints(i);
                end
                IntvUb      = BackboneGaussPoints(i+1);

                % Process the centroids
                idxPoints   = (Points >= IntvLb) & (Points <= IntvUb);
                if any(idxPoints)
                    PointToBackboneExtended(1, idxPoints) = i;
                end
                
            end

            % Store a compact representation of the above interval as detailed in the properties
            [Y, C]  = findRepetitions(PointToBackboneExtended);
            % Return the mapping
            M       = [C', Y'];
        end

        % Setter for the points where the backbone is evaluated
        function setPoints(obj, Points)
            arguments
                obj                 (1, 1) LVPBackbone
                Points              (:, 1) double 
            end

            % Store the new points along the backbone
            obj.Points                                              = Points;
            obj.nPoints                                             = length(obj.Points); 
            % Find the unique points among the provided points
            [obj.UniquePoints, ~, obj.UniquePointsToPoints]         = unique(Points);
            obj.nUniquePoints                                       = length(obj.UniquePoints);
            obj.UniquePointsToGaussPoints                           = obj.MapPointsToGaussianPoints(obj.UniquePoints);
            obj.UniquePointsToGaussPointsExpanded                   = repelem(obj.UniquePointsToGaussPoints(:, 1), obj.UniquePointsToGaussPoints(:, 2));
            
            
            % Initialize the variables of the backbone points
            obj.Strain                          = repmat(obj.ReferenceStrain, [1, obj.nPoints]);
            obj.dStrain                         = zeros(6, obj.nPoints);
            obj.ddStrain                        = zeros(6, obj.nPoints);
            obj.JStrain                         = zeros(6, obj.n, obj.nPoints);
            obj.dJStrain                        = zeros(6, obj.n, obj.nPoints);
            obj.Strain_ds                       = zeros(6, obj.nPoints);
            obj.dStrain_ds                      = zeros(6, obj.nPoints);
            obj.Elongation                      = zeros(obj.nPoints, 1);
            obj.JElongation                     = zeros(obj.n, obj.nPoints);
            obj.dElongation                     = zeros(obj.nPoints, 1);
            obj.ddElongation                    = zeros(obj.nPoints, 1);
            obj.Pose                            = zeros(4, 4, obj.nPoints);
            obj.Pose(1:4, 1:4, 1:obj.nPoints)   = repmat(eye(4), 1, 1, obj.nPoints);
            obj.dPose                           = zeros(6, obj.nPoints);
            obj.ddPose                          = zeros(6, obj.nPoints);
            obj.JPose                           = zeros(6, obj.n, obj.nPoints);
            
            % Initialize the variables of the unique backbone points
            obj.UniqueStrain                          = repmat(obj.ReferenceStrain, [1, obj.nUniquePoints]);
            obj.UniquedStrain                         = zeros(6, obj.nUniquePoints);
            obj.UniqueddStrain                        = zeros(6, obj.nUniquePoints);
            obj.UniqueJStrain                         = zeros(6, obj.n, obj.nUniquePoints);
            obj.UniquedJStrain                        = zeros(6, obj.n, obj.nUniquePoints);
            obj.UniqueStrain_ds                       = zeros(6, obj.nUniquePoints);
            obj.UniquedStrain_ds                      = zeros(6, obj.nUniquePoints);
            obj.UniqueElongation                      = zeros(obj.nUniquePoints, 1);
            obj.UniqueJElongation                     = zeros(obj.n, obj.nUniquePoints);
            obj.UniquedElongation                     = zeros(obj.nUniquePoints, 1);
            obj.UniqueddElongation                    = zeros(obj.nUniquePoints, 1);
            obj.UniquePose                            = repmat(eye(4), 1, 1, obj.nUniquePoints);
            obj.UniquedPose                           = zeros(6, obj.nUniquePoints);
            obj.UniqueddPose                          = zeros(6, obj.nUniquePoints);
            obj.UniqueJPose                           = zeros(6, obj.n, obj.nUniquePoints);
            obj.UniqueJStrainZanna1                   = zeros(6, obj.n, obj.nUniquePoints);
            obj.UniqueJStrainZanna2                   = zeros(6, obj.n, obj.nUniquePoints);
            obj.UniquedJStrainZanna                   = zeros(6, obj.n, obj.nUniquePoints, 2);

            % Compute the Zanna quadrature points for the unique points
            Zanna1 = zeros(1, obj.nUniquePoints);
            Zanna2 = zeros(1, obj.nUniquePoints);
            s               = obj.UniquePoints;
            sToGaussPoints  = obj.UniquePointsToGaussPoints;
            idxS            = 1;
            for i = 1:length(sToGaussPoints)
                % Get index of the current Gauss point
                idx     = sToGaussPoints(i, 1);
                % Get the ending index of the unique points between the interval [GaussPoints(idx), GaussPoints(idx+1)]
                idxE    = idxS + sToGaussPoints(i, 2) - 1;
                if idx  == 0
                    S = 0;
                else
                    S = obj.GaussPoints(idx);
                end
                % Compute the distance of the Gaussian points from S
                h                 = s(idxS:idxE) - S;
                Zanna1(idxS:idxE) =  S + h/2 - sqrt(3)*h/6;
                Zanna2(idxS:idxE) =  S + h/2 + sqrt(3)*h/6;

                % Prepare for the next iteration
                idxS    = idxE + 1;
            end
            
            % Pre-evaluate the strain basis at the unique points and at their zanna quadrature points
            for i = 1:obj.nUniquePoints
                [obj.UniqueJStrain(1:6, 1:obj.n, i), obj.UniquedJStrain(1:6, 1:obj.n, i)] = obj.StrainBasis(obj.UniquePoints(i));
                [JxiGaussZanna1, dJxiGaussZanna1]           = obj.StrainBasis(Zanna1(i));
                [JxiGaussZanna2, dJxiGaussZanna2]           = obj.StrainBasis(Zanna2(i));
                obj.UniqueJStrainZanna1(1:6, 1:obj.n, i)    = JxiGaussZanna1;
                obj.UniqueJStrainZanna2(1:6, 1:obj.n, i)    = JxiGaussZanna2;
                obj.UniquedJStrainZanna(1:6, 1:obj.n, i, 1) = dJxiGaussZanna1;
                obj.UniquedJStrainZanna(1:6, 1:obj.n, i, 2) = dJxiGaussZanna2;
            end
        end
        
        %Strain basis. By default a LVPBackbone has no DOF, thus returns and empty basis
        function [Phi, dPhi] = StrainBasis(obj, s)
            arguments (Input)
                obj (1, 1) LVPBackbone
                s   (1, 1) 
            end
            
            arguments (Output)
                Phi     (6, :)
                dPhi    (6, :)
            end

            % Output preallocation for code generation
            Phi     = zeros(6, obj.n);
            dPhi    = zeros(6, obj.n);

            % Iterate over the primitives to update the strain basis
            for i = 1:LVPBody.MaxPrimitivesNumber
                if i <= obj.NPrimitives
                    if ~isa(obj.Primitives{i}, "LVPPrimitive")
                        continue;
                    end
                    % Check for the primitives that modify the strain
                    if obj.PrimitiveModifiesBackbone(i)
                        [Phi(1:6, obj.QIdx(i, 1):obj.QIdx(i, 2)), dPhi(1:6, obj.QIdx(i, 1):obj.QIdx(i, 2))] = obj.Primitives{i}.StrainBasis(s);
                    end
                end
            end
        end


        % Compute the status of the backvone at the unique query points
        function UpdateUniquePoints(obj, q, dq, ddq)
            arguments (Input)
                obj (1, 1) LVPBackbone
                q   (:, 1) = zeros(obj.n, 1)
                dq  (:, 1) = zeros(obj.n, 1)
                ddq (:, 1) = zeros(obj.n, 1)
            end

            % Update the elongation quantities
            GaussPointsVal1         = obj.UniquePointsToGaussPoints(1:end-1, 1) + 1;
            %if GaussPointsVal1(1) == 0
            %    GaussPointsVal1 = GaussPointsVal1 + 1;
            %end
            GaussPointsVal2         = obj.UniquePointsToGaussPoints(:, 1) + 1;
            %if GaussPointsVal2(1) == 0
            %    GaussPointsVal2     = GaussPointsVal2 + 1;
            %end
            GaussPointsRep          = obj.UniquePointsToGaussPoints(:, 2);
            if obj.UniquePointsToGaussPoints(1, 1) == 0
                PointsToGaussPointsRep  = repelem(GaussPointsVal2, GaussPointsRep);
            else
                PointsToGaussPointsRep  = repelem(obj.UniquePointsToGaussPoints(:, 1), GaussPointsRep);
            end

            % Normalize the location of the points along the backbone
            s       = obj.UniquePoints;
            IntvLb  = repelem([0, obj.GaussPoints(GaussPointsVal1)], 1, GaussPointsRep);
            IntvUb  = repelem(obj.GaussPoints(GaussPointsVal2), 1, GaussPointsRep);
            s_norm  = (s-IntvLb)./(IntvUb-IntvLb);
            
            % Compute the interpolating terms
            DeltaL1     = repelem([0, obj.ElongationGauss(GaussPointsVal1)], 1, GaussPointsRep);
            DeltaL2     = repelem(obj.ElongationGauss(GaussPointsVal2), 1, GaussPointsRep);
            dDeltaL1    = repelem([0, obj.dElongationGauss(GaussPointsVal1)], 1, GaussPointsRep);
            dDeltaL2    = repelem(obj.dElongationGauss(GaussPointsVal2), 1, GaussPointsRep);
            ddDeltaL1   = repelem([0, obj.ddElongationGauss(GaussPointsVal1)], 1, GaussPointsRep);
            ddDeltaL2   = repelem(obj.ddElongationGauss(GaussPointsVal2), 1, GaussPointsRep);
            JL1Gauss    = cat(2, zeros(obj.n, 1), obj.JElongationGauss(1:obj.n, GaussPointsVal1));
            JL2Gauss    = obj.JElongationGauss(1:obj.n, GaussPointsVal2);
            JL1         = JL1Gauss(:, PointsToGaussPointsRep);
            JL2         = JL2Gauss(:, PointsToGaussPointsRep);
            
            % Update the elongation terms
            obj.UniqueElongation    = DeltaL1.*(1-s_norm) + DeltaL2.*s_norm;
            obj.UniquedElongation   = dDeltaL1.*(1-s_norm) + dDeltaL2.*s_norm;
            obj.UniqueddElongation  = ddDeltaL1.*(1-s_norm) + ddDeltaL2.*s_norm;
            obj.UniqueJElongation   = JL1.*(1-s_norm) + JL2.*s_norm;
            
            % Update the strain. 
            % TODO: We can use the Zanna quadrature points to avoid repeating these computations
            obj.UniqueStrain(1:6, :)                      = reshape(pagemtimes(obj.UniqueJStrain, q), 6, []) + obj.ReferenceStrain;
            obj.UniquedStrain(1:6, :)                     = reshape(pagemtimes(obj.UniqueJStrain, dq), 6, []);
            obj.UniqueddStrain(1:6, :)                    = reshape(pagemtimes(obj.UniqueJStrain, ddq), 6, []);
            obj.UniqueStrain_ds(1:6, :)                   = reshape(pagemtimes(obj.UniquedJStrain, q), 6, []);
            obj.UniquedStrain_ds(1:6, :)                  = reshape(pagemtimes(obj.UniquedJStrain, dq), 6, []);
            
            % Update the kinematics
            obj.UpdateKinematicsAtUniquePoints(q, dq, ddq);
            
        end

        % Compute the status of the backbone at the query points (Points)
        % Each element of the rows of sToGaussPoints has the following meaning:
        % 1) A Gaussian point x of the backbone
        % 2) The number of consecutive elements in s that belong the interval starting with x and ending with x+1, where x+1 is the following Gaussian point to x
        % In the future, this method should allow the computation NOT approximated by using the Kinematics method of the GVSBody.
        % This however requires to vectorize the strain methods of the GVSBody function and also the strain functions
        function UpdatePoints(obj, q, dq, ddq)
            % Update the pose at the points
            arguments (Input)
                obj (1, 1) LVPBackbone
                q   (:, 1) = zeros(obj.n, 1)
                dq  (:, 1) = zeros(obj.n, 1)
                ddq (:, 1) = zeros(obj.n, 1)
            end
            
            % Update the status of the unique points
            obj.UpdateUniquePoints(q, dq, ddq);

            % Map the unique points to the points specified by the user
            obj.Strain(1:6, 1:obj.nPoints)                  = obj.UniqueStrain(1:6, obj.UniquePointsToPoints);
            obj.dStrain(1:6, 1:obj.nPoints)                 = obj.UniquedStrain(1:6, obj.UniquePointsToPoints);
            obj.ddStrain(1:6, 1:obj.nPoints)                = obj.UniqueddStrain(1:6, obj.UniquePointsToPoints);
            obj.JStrain(1:6, 1:obj.n, 1:obj.nPoints)        = obj.UniqueJStrain(1:6, 1:obj.n, obj.UniquePointsToPoints);
            obj.dJStrain(1:6, 1:obj.n, 1:obj.nPoints)       = obj.UniquedJStrain(1:6, 1:obj.n, obj.UniquePointsToPoints);
            obj.Strain_ds(1:6, 1:obj.nPoints)               = obj.UniqueStrain_ds(1:6, obj.UniquePointsToPoints);
            obj.dStrain_ds(1:6, 1:obj.nPoints)              = obj.UniquedStrain_ds(1:6, obj.UniquePointsToPoints);
            obj.Elongation(1:obj.nPoints)                   = obj.UniqueElongation(obj.UniquePointsToPoints);
            obj.JElongation(1:obj.n, 1:obj.nPoints)         = obj.UniqueJElongation(1:obj.n, obj.UniquePointsToPoints);
            obj.dElongation(1:obj.nPoints)                  = obj.UniquedElongation(obj.UniquePointsToPoints);
            obj.ddElongation(1:obj.nPoints)                 = obj.UniqueddElongation(obj.UniquePointsToPoints);
            obj.Pose(1:4, 1:4, 1:obj.nPoints)               = obj.UniquePose(1:4, 1:4, obj.UniquePointsToPoints);
            obj.dPose(1:6, 1:obj.nPoints)                   = obj.UniquedPose(1:6, obj.UniquePointsToPoints);
            obj.ddPose(1:6, 1:obj.nPoints)                  = obj.UniqueddPose(1:6, obj.UniquePointsToPoints);
            obj.JPose(1:6, 1:obj.n, 1:obj.nPoints)          = obj.UniqueJPose(1:6, 1:obj.n, obj.UniquePointsToPoints);
        end

        function [JP] = PoseJacobian(obj, p)
            arguments (Input)
                obj (1, 1) LVPBackbone
                p   (1, 1) LVPPrimitive
            end

            arguments (Output)
                JP (6, :, :)
            end

            % Output preallocation
            JP       = zeros(6, obj.n, obj.nPoints);% Jacobian of the pose w.r.t the configuration variables up to the current primitive

            % Assign the Jacobian by looking among the primitives
            for i = 1:LVPBody.MaxPrimitivesNumber
                piP = obj.Primitives{i};
                if isa(piP, "LVPPrimitive")
                    if piP == p
                        idxE                            = obj.QIdx(i, 2);
                        JP(1:6, 1:idxE, 1:obj.nPoints)  = obj.JPose(1:6, 1:idxE, 1:obj.nPoints);
                        break;
                    end
                end
            end
        end

        function [Jxi, dJxi]          = StrainJacobian(obj, p)
            arguments (Input)
                obj (1, 1) LVPBackbone
                p   (1, 1) LVPPrimitive
            end

            arguments (Output)
                Jxi     (6, :, :)
                dJxi    (6, :, :)
            end

            % Output preallocation
            Jxi     = zeros(6, obj.n, obj.nPoints);% Jacobian of the pose w.r.t the configuration variables up to the current primitive
            dJxi    = zeros(6, obj.n, obj.nPoints);% Derivative w.r.t. s of the Jacobian of the pose w.r.t the configuration variables up to the current primitive

            % Assign the Jacobian by looking among the primitives
            for i = 1:LVPBody.MaxPrimitivesNumber
                piP = obj.Primitives{i};
                if isa(piP, "LVPPrimitive")
                    if piP == p
                        idxE                                = obj.QIdx(i, 2);
                        Jxi(1:6, 1:idxE, 1:obj.nPoints)     = obj.JStrain(1:6, 1:idxE, 1:obj.nPoints);
                        dJxi(1:6, 1:idxE, 1:obj.nPoints)    = obj.dJStrain(1:6, 1:idxE, 1:obj.nPoints);
                        break;
                    end
                end
            end
        end

        function [JL, JLStrain] = ElongationJacobian(obj, p)
            arguments (Input)
                obj (1, 1) LVPBackbone
                p   (1, 1) LVPPrimitive
            end

            arguments (Output)
                JL (:, :)
                JLStrain (:, :)
            end

            % Output preallocation
            JL       = zeros(obj.n, obj.nPoints);% Jacobian of the elongation
            JLStrain = zeros(obj.n, obj.nPoints);% Jacobian of the elongation strain

            % Assign the Jacobian by looking among the primitives
            for i = 1:LVPBody.MaxPrimitivesNumber
                piP = obj.Primitives{i};
                if isa(piP, "LVPPrimitive")
                    if piP == p
                        idxE = obj.QIdx(i, 2);
                        JL(1:idxE, 1:obj.nPoints)       = obj.JElongation(1:idxE, 1:obj.nPoints);
                        JLStrainQ                       = reshape(obj.JStrain(6, 1:obj.n, 1:obj.nPoints), obj.n, obj.nPoints);
                        JLStrain(1:idxE, 1:obj.nPoints) = JLStrainQ(1:idxE, 1:obj.nPoints);
                        break;
                    end
                end
            end
        end
        
        function Update(obj, q, dq, ddq)
            %Update the status of the backbone.
            %Args:
            %   obj (LVPBakcbone)       : LVPBackbone object
            %   q   ([double], [sym])   : Configuration variables
            %   dq  ([double], [sym])   : First order time derivative of the configuration variables
            %   ddq ([double], [sym])   : Second order time derivative of the configuration variables

            % Update the kinematics of the backbone
            obj.UpdateKinematicsGauss(q, dq, ddq);

            % Update the other values associated to the Gaussian points
            obj.xiGauss(1:6, :)     = pagemtimes(obj.JxiGauss, q) + obj.ReferenceStrain;
            obj.dxiGauss(1:6, :)    = pagemtimes(obj.JxiGauss, dq);
            obj.ddxiGauss(1:6, :)   = pagemtimes(obj.JxiGauss, ddq);
            obj.xiGauss_ds(1:6, :)  = pagemtimes(obj.dJxiGauss, q);
            obj.dxiGauss_ds(1:6, :) = pagemtimes(obj.dJxiGauss, dq);

            % Update also the other values associated to the Gaussian points.
            for i = 1:obj.NGaussPointsInt
                JGaussT              = obj.JxiGauss(1:6, 1:obj.n, i);
                % Update the elongation
                if i == 1
                    obj.ElongationGauss(i)              = obj.GaussWeightsBackbonePoints(i)*obj.xiGauss(6, i);
                    obj.JElongationGauss(1:obj.n, i)    = obj.GaussWeightsBackbonePoints(i)*JGaussT(6, 1:obj.n)';
                    obj.dElongationGauss(i)             = obj.GaussWeightsBackbonePoints(i)*obj.dxiGauss(6, i);
                    obj.ddElongationGauss(i)            = obj.GaussWeightsBackbonePoints(i)*obj.ddxiGauss(6, i);
                else
                    obj.ElongationGauss(i)              = obj.ElongationGauss(i-1) + obj.GaussWeightsBackbonePoints(i)*obj.xiGauss(6, i);
                    obj.JElongationGauss(1:obj.n, i)    = obj.JElongationGauss(1:obj.n, i-1) + obj.GaussWeightsBackbonePoints(i)*JGaussT(6, 1:obj.n)';
                    obj.dElongationGauss(i)             = obj.dElongationGauss(i-1) + obj.GaussWeightsBackbonePoints(i)*obj.dxiGauss(6, i);
                    obj.ddElongationGauss(i)            = obj.ddElongationGauss(i-1) + obj.GaussWeightsBackbonePoints(i)*obj.ddxiGauss(6, i);
                end
            end
            obj.ElongationGauss = obj.ElongationGauss - obj.GaussPoints;

            % Update the status of the backbone at the query Points
            obj.UpdatePoints(q, dq, ddq);
        end
        
        % Update the kinematics for the Gaussian points
        function UpdateKinematicsGauss(obj, q, dq, ddq)
            arguments (Input)
                obj (1, 1) LVPBackbone
                q   (:, 1) = zeros(obj.n, 1);
                dq  (:, 1) = zeros(obj.n, 1);
                ddq (:, 1) = zeros(obj.n, 1);
            end

            % Evaluate the strain and its time derivatives at the zanna quadrature points
            Jxi_1   = squeeze(obj.JxiGaussZanna(1:6, 1:obj.n, 1:obj.NGaussPointsInt, 1));
            Jxi_2   = squeeze(obj.JxiGaussZanna(1:6, 1:obj.n, 1:obj.NGaussPointsInt, 2));
            xi_1    = squeeze(pagemtimes(Jxi_1, q)) + obj.ReferenceStrain;
            xi_2    = squeeze(pagemtimes(Jxi_2, q)) + obj.ReferenceStrain;
            dxi_1   = squeeze(pagemtimes(Jxi_1, dq));
            dxi_2   = squeeze(pagemtimes(Jxi_2, dq));
            ddxi_1  = squeeze(pagemtimes(Jxi_1, ddq));
            ddxi_2  = squeeze(pagemtimes(Jxi_2, ddq));

            % Compute the relative distance between Gaussian points
            h       = [obj.GaussPoints(1), diff(obj.GaussPoints)];

            % Precompute useful quantities
            adxi_1   = ad(xi_1);
            adxi_2   = ad(xi_2);
            Sqrth212 = reshape((sqrt(3)*h.^2)/12, 1, 1, []);

            % Compute the Magnus expansion of the strain
            Omega                   = (h/2).*(xi_1 + xi_2) + ((sqrt(3)*h.^2)/12).*reshape(pagemtimes(adxi_1, reshape(xi_2, 6, 1, [])), 6, []);
            OmegaHat                = zeros(4, 4, obj.NGaussPointsInt, "like", q);
            OmegaHat(1:3, 1:3, :)   = skew(Omega(1:3, :));
            OmegaHat(1:3, 4,   :)   = Omega(4:6, :);

            % Compute the first order time deritvative of the Magnus expansion of the strain
            PhiOmega= pagemtimes(reshape(h/2, 1, 1, []), Jxi_1 + Jxi_2) + ...
                      pagemtimes(Sqrth212, pagemtimes(adxi_1, Jxi_2) - pagemtimes(adxi_2, Jxi_1));
            dOmega  = squeeze(pagemtimes(PhiOmega, dq));

            % Compute the second order time derivative of the Magnus expansion of the strain
            ddOmega = (h/2).*(ddxi_1 + ddxi_2) + ....
                            squeeze(pagemtimes(Sqrth212, pagemtimes(ad(ddxi_1), reshape(xi_2, 6, 1, [])) + ...
                                                         pagemtimes(2*ad(dxi_1), reshape(dxi_2, 6, 1, [])) + ...
                                                         pagemtimes(adxi_1, reshape(ddxi_2, 6, 1, []))));

            % Compute the tanget map
            TOmega  = Tanexpmat(Omega);

            % Compute the exponential mat
            ExpOmega        = expmat(OmegaHat);

            % Compute the inverse adjoint of the exponential map
            invAdExpOmega   = invAd(ExpOmega);

            % Compute useful products
            TOmegaPhiOmega = pagemtimes(TOmega, PhiOmega);
            TOmegaddOmega  = squeeze(pagemtimes(TOmega, reshape(ddOmega, 6, 1, [])));
            adTOmegadOmega = ad(squeeze(pagemtimes(TOmega, reshape(dOmega, 6, 1, []))));

            % Loop over the Gaussian points
            for i = 1:obj.NGaussPointsInt
                if i == 1
                    obj.gGauss(:, :, i)          = ExpOmega(1:4, 1:4, i);
                    obj.JEtaGaussLocal(:, :, i)  = invAdExpOmega(1:6, 1:6, i)*TOmegaPhiOmega(1:6, 1:obj.n, i);
                    obj.EtaGauss(:, i)           = obj.JEtaGaussLocal(:, :, i)*dq;
                    obj.dEtaGaussLocal(:, i)     = invAdExpOmega(1:6, 1:6, i)*(-adTOmegadOmega(1:6, 1:6, i)*obj.EtaGauss(:, i) + TOmegaddOmega(1:6, i));
                else
                    obj.gGauss(:, :, i)          = obj.gGauss(:, :, i-1)*ExpOmega(1:4, 1:4, i);
                    obj.JEtaGaussLocal(:, :, i)  = invAdExpOmega(1:6, 1:6, i)*(obj.JEtaGaussLocal(:, :, i-1) + TOmegaPhiOmega(1:6, 1:obj.n, i));
                    obj.EtaGauss(:, i)           = obj.JEtaGaussLocal(:, :, i)*dq;
                    obj.dEtaGaussLocal(:, i)     = invAdExpOmega(1:6, 1:6, i)*(obj.dEtaGaussLocal(:, i-1) -adTOmegadOmega(1:6, 1:6, i)*obj.EtaGauss(:, i) + TOmegaddOmega(1:6, i));
                end
            end

            % Rotate all the quantities in the global frame
            R                                   = obj.gGauss(1:3, 1:3, :);
            omega                               = squeeze(pagemtimes(R, reshape(obj.EtaGauss(1:3, :), 3, 1, [])));
            dR                                  = pagemtimes(skew(omega), R);
            % NOTE: We first compute the acceleration because the velocity must be in local coordinates!
            obj.dEtaGauss(1:3, :)               = squeeze(pagemtimes(dR, reshape(obj.EtaGauss(1:3, :), 3, 1, [])) + pagemtimes(R, reshape(obj.dEtaGaussLocal(1:3, :), 3, 1, [])));
            obj.dEtaGauss(4:6, :)               = squeeze(pagemtimes(dR, reshape(obj.EtaGauss(4:6, :), 3, 1, [])) + pagemtimes(R, reshape(obj.dEtaGaussLocal(4:6, :), 3, 1, [])));
            % Update the velocity
            obj.JEtaGauss(1:3, 1:obj.n, :)      = pagemtimes(R, reshape(obj.JEtaGaussLocal(1:3, 1:obj.n, :), 3, obj.n, []));
            obj.JEtaGauss(4:6, 1:obj.n, :)      = pagemtimes(R, reshape(obj.JEtaGaussLocal(4:6, 1:obj.n, :), 3, obj.n, []));
            obj.EtaGauss                        = squeeze(pagemtimes(obj.JEtaGauss, dq));

            % Update the kinematics of the tip since these are used in the other computations
            obj.T_                  = obj.gGauss(1:4, 1:4, end);
            obj.v_rel_              = obj.EtaGauss(4:6, end);
            obj.omega_rel_          = obj.EtaGauss(1:3, end);
            obj.a_rel_              = obj.dEtaGauss(4:6, end);
            obj.domega_rel_         = obj.dEtaGauss(1:3, end);
            %NOTE: The jacobians are expressed in the tip frame
            RT                      = obj.T_(1:3, 1:3);
            obj.v_par_              = RT*obj.JEtaGauss(4:6, 1:obj.n, end);
            obj.omega_par_          = RT*obj.JEtaGauss(1:3, 1:obj.n, end);
        end

    end

    methods (Access = private)

        function UpdateKinematicsAtUniquePoints(obj, q, dq, ddq)
            arguments (Input)
                obj (1, 1) LVPBackbone
                q   (:, 1) = zeros(obj.n, 1)
                dq  (:, 1) = zeros(obj.n, 1)
                ddq  (:, 1)= zeros(obj.n, 1)
            end

            % Evaluate the strain and its time derivatives at the Zanna quadrature points
            xi_1    = squeeze(pagemtimes(obj.UniqueJStrainZanna1, q)) + obj.ReferenceStrain;
            xi_2    = squeeze(pagemtimes(obj.UniqueJStrainZanna2, q)) + obj.ReferenceStrain;
            dxi_1   = squeeze(pagemtimes(obj.UniqueJStrainZanna1, dq));
            dxi_2   = squeeze(pagemtimes(obj.UniqueJStrainZanna2, dq));
            ddxi_1  = squeeze(pagemtimes(obj.UniqueJStrainZanna1, ddq));
            ddxi_2  = squeeze(pagemtimes(obj.UniqueJStrainZanna2, ddq));

            % Compute the relative distance between the points and the Gaussian point
            PointsToGaussPointsRep  = repelem(obj.UniquePointsToGaussPoints(:, 1) + 1, obj.UniquePointsToGaussPoints(:, 2));
            sGauss                  = [0, obj.GaussPoints(1:obj.NGaussPoints)];
            h                       = obj.UniquePoints-sGauss(PointsToGaussPointsRep);

            % Precompute useful quantities
            adxi_1   = ad(xi_1);
            adxi_2   = ad(xi_2);
            Sqrth212 = reshape((sqrt(3)*h.^2)/12, 1, 1, []);

            % Compute the Magnus expansion of the strain
            Omega                   = (h/2).*(xi_1 + xi_2) + ((sqrt(3)*h.^2)/12).*reshape(pagemtimes(adxi_1, reshape(xi_2, 6, 1, [])), 6, []);
            OmegaHat                = zeros(4, 4, obj.nUniquePoints, "like", q);
            OmegaHat(1:3, 1:3, :)   = skew(Omega(1:3, :));
            OmegaHat(1:3, 4,   :)   = Omega(4:6, :);

            % Compute the first order time deritvative of the Magnus expansion of the strain
            PhiOmega= pagemtimes(reshape(h/2, 1, 1, []), obj.UniqueJStrainZanna1 + obj.UniqueJStrainZanna2) + ...
                      pagemtimes(Sqrth212, pagemtimes(adxi_1, obj.UniqueJStrainZanna2) - pagemtimes(adxi_2, obj.UniqueJStrainZanna1));
            dOmega  = squeeze(pagemtimes(PhiOmega, dq));

            % Compute the second order time derivative of the Magnus expansion of the strain
            ddOmega = (h/2).*(ddxi_1 + ddxi_2) + ....
                            squeeze(pagemtimes(Sqrth212, pagemtimes(ad(ddxi_1), reshape(xi_2, 6, 1, [])) + ...
                                                         pagemtimes(2*ad(dxi_1), reshape(dxi_2, 6, 1, [])) + ...
                                                         pagemtimes(adxi_1, reshape(ddxi_2, 6, 1, []))));

            % Compute the tanget map
            TOmega  = Tanexpmat(Omega);

            % Compute the exponential mat
            ExpOmega        = expmat(OmegaHat);

            % Compute the inverse adjoint of the exponential map
            invAdExpOmega   = invAd(ExpOmega);

            % Compute useful products
            TOmegaPhiOmega = pagemtimes(TOmega, PhiOmega);
            TOmegaddOmega  = pagemtimes(TOmega, reshape(ddOmega, 6, 1, []));
            adTOmegadOmega = ad(squeeze(pagemtimes(TOmega, reshape(dOmega, 6, 1, []))));

            % Expand the quantities of the Gaussian points
            GG             = cat(3, eye(4), obj.gGauss(:, :, 1:end-1));
            gGauss         = GG(:, :, PointsToGaussPointsRep);
            %gGauss         = repelem(cat(3, eye(4), obj.gGauss(:, :, 1:end-1)), 1, 1, PointsToGaussPointsRep);
            JG             = cat(3, zeros(6, obj.n), obj.JEtaGaussLocal(:, :, 1:end-1));
            JEtaGauss      = JG(:, :, PointsToGaussPointsRep);
            EtaGauss       = squeeze(pagemtimes(JEtaGauss, dq));
            dEtaG          = cat(2, zeros(6, 1), obj.dEtaGaussLocal(:, 1:end-1));
            dEtaGauss      = dEtaG(:, PointsToGaussPointsRep);

            % Compute the kinematics terms
            obj.UniquePose    = pagemtimes(gGauss, ExpOmega);
            obj.UniqueJPose   = pagemtimes(invAdExpOmega, JEtaGauss + TOmegaPhiOmega);
            EtaPages          = pagemtimes(obj.UniqueJPose, dq);
            obj.UniquedPose   = squeeze(EtaPages);
            obj.UniqueddPose  = squeeze(pagemtimes(invAdExpOmega, reshape(dEtaGauss, 6, 1, []) - pagemtimes(adTOmegadOmega, EtaPages) + TOmegaddOmega));

            % Rotate the terms in the base frame
            R                                   = obj.UniquePose(1:3, 1:3, :);
            omega                               = squeeze(pagemtimes(R, reshape(obj.UniquedPose(1:3, :), 3, 1, [])));
            dR                                  = pagemtimes(skew(omega), R);
            % NOTE: First update the acceleration as in the below expression Eta must be expressed in the local frame
            obj.UniqueddPose(1:3, :)            = squeeze(pagemtimes(dR, reshape(obj.UniquedPose(1:3, :), 3, 1, [])) + pagemtimes(R, reshape(obj.UniqueddPose(1:3, :), 3, 1, [])));
            obj.UniqueddPose(4:6, :)            = squeeze(pagemtimes(dR, reshape(obj.UniquedPose(4:6, :), 3, 1, [])) + pagemtimes(R, reshape(obj.UniqueddPose(4:6, :), 3, 1, [])));
            obj.UniqueJPose(1:3, :, :)          = pagemtimes(R, reshape(obj.UniqueJPose(1:3, :, :), 3, obj.n, []));
            obj.UniqueJPose(4:6, :, :)          = pagemtimes(R, reshape(obj.UniqueJPose(4:6, :, :), 3, obj.n, []));
            obj.UniquedPose                     = squeeze(pagemtimes(obj.UniqueJPose, dq));
        end
        
        function [g, JEta, Eta, dEta] = UpdateKinematicsAtS(sGauss, gGauss, JEtaGauss, dEtaGauss, s, Jxi_1, Jxi_2, q, dq, ddq, ReferenceStrain)
            arguments (Input)
                sGauss          (1, 1)     double
                gGauss          (4, 4)     double
                JEtaGauss       (6, :)     double
                dEtaGauss       (6, 1)     double
                s               (1, :)     double
                Jxi_1           (6, :, :)  double
                Jxi_2           (6, :, :)  double
                q               (:, 1)     double
                dq              (:, 1)     double
                ddq             (:, 1)     double
                ReferenceStrain (6, 1)     double
            end

            arguments (Output)
                g            (4, 4, :) double
                JEta         (6, :, :) double
                Eta          (6, :)    double
                dEta         (6, :)    double
            end
            
            % Size of q
            n   = size(q, 1);
            % Get the number of points
            ns  = size(s, 2);
            
            % Evaluate the strain and its time derivatives at the Zanna quadrature points
            xi_1    = squeeze(pagemtimes(Jxi_1, q)) + ReferenceStrain;
            xi_2    = squeeze(pagemtimes(Jxi_2, q)) + ReferenceStrain;
            dxi_1   = squeeze(pagemtimes(Jxi_1, dq));
            dxi_2   = squeeze(pagemtimes(Jxi_2, dq));
            ddxi_1  = squeeze(pagemtimes(Jxi_1, ddq));
            ddxi_2  = squeeze(pagemtimes(Jxi_2, ddq));

            % Compute the relative distance between the points and the Gaussian point
            h       = s-sGauss;

            % Precompute useful quantities
            adxi_1   = ad(xi_1);
            adxi_2   = ad(xi_2);
            Sqrth212 = reshape((sqrt(3)*h.^2)/12, 1, 1, []);

            % Compute the Magnus expansion of the strain
            Omega                   = (h/2).*(xi_1 + xi_2) + ((sqrt(3)*h.^2)/12).*reshape(pagemtimes(adxi_1, reshape(xi_2, 6, 1, [])), 6, []);
            OmegaHat                = zeros(4, 4, ns, "like", q);
            OmegaHat(1:3, 1:3, :)   = skew(Omega(1:3, :));
            OmegaHat(1:3, 4,   :)   = Omega(4:6, :);

            % Compute the first order time deritvative of the Magnus expansion of the strain
            PhiOmega= pagemtimes(reshape(h/2, 1, 1, []), Jxi_1 + Jxi_2) + ...
                      pagemtimes(Sqrth212, pagemtimes(adxi_1, Jxi_2) - pagemtimes(adxi_2, Jxi_1));
            dOmega  = squeeze(pagemtimes(PhiOmega, dq));

            % Compute the second order time derivative of the Magnus expansion of the strain
            ddOmega = (h/2).*(ddxi_1 + ddxi_2) + ....
                            squeeze(pagemtimes(Sqrth212, pagemtimes(ad(ddxi_1), reshape(xi_2, 6, 1, [])) + ...
                                                         pagemtimes(2*ad(dxi_1), reshape(dxi_2, 6, 1, [])) + ...
                                                         pagemtimes(adxi_1, reshape(ddxi_2, 6, 1, []))));

            % Compute the tanget map
            TOmega  = Tanexpmat(Omega);

            % Compute the exponential mat
            ExpOmega        = expmat(OmegaHat);

            % Compute the inverse adjoint of the exponential map
            invAdExpOmega   = invAd(ExpOmega);

            % Compute useful products
            TOmegaPhiOmega = pagemtimes(TOmega, PhiOmega);
            TOmegaddOmega  = pagemtimes(TOmega, reshape(ddOmega, 6, 1, []));
            adTOmegadOmega = ad(squeeze(pagemtimes(TOmega, reshape(dOmega, 6, 1, []))));

            % Compute the kinematics terms
            g                 = pagemtimes(gGauss, ExpOmega);
            JEta              = pagemtimes(invAdExpOmega, JEtaGauss + TOmegaPhiOmega);
            EtaPages          = pagemtimes(JEta, dq);
            Eta               = squeeze(EtaPages);
            dEta              = squeeze(pagemtimes(invAdExpOmega, dEtaGauss - pagemtimes(adTOmegadOmega, EtaPages) + TOmegaddOmega));

            % Rotate the terms in the base frame
            R                                   = g(1:3, 1:3, :);
            omega                               = squeeze(pagemtimes(R, reshape(Eta(1:3, :), 3, 1, [])));
            dR                                  = pagemtimes(skew(omega), R);
            % NOTE: First update the acceleration as in the below expression Eta must be expressed in the local frame
            dEta(1:3, :)                        = squeeze(pagemtimes(dR, reshape(Eta(1:3, :), 3, 1, [])) + pagemtimes(R, reshape(dEta(1:3, :), 3, 1, [])));
            dEta(4:6, :)                        = squeeze(pagemtimes(dR, reshape(Eta(4:6, :), 3, 1, [])) + pagemtimes(R, reshape(dEta(4:6, :), 3, 1, [])));
            Eta(1:3, :)                         = omega;
            Eta(4:6, :)                         = squeeze(pagemtimes(R, reshape(Eta(4:6, :), 3, 1, [])));
            JEta(1:3, 1:n, :)                   = pagemtimes(R, reshape(JEta(1:3, 1:n, :), 3, n, []));
            JEta(4:6, 1:n, :)                   = pagemtimes(R, reshape(JEta(4:6, 1:n, :), 3, n, []));
        end
    end
end

