classdef LVPBackbone < GVSBody
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
        JStrain     (:, :) double
        % Derivative of the strain w.r.t. s (curvilinear abscissa)
        Strain_ds   (6, :) double
        % Derivative of the first order time derivative of the strain w.r.t. s (curvilinear abscissa)
        dStrain_ds  (6, :) double
        % Derivative w.r.t. s of the strain Jacobian
        dJStrain    (:, :) double
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
        GaussPointsBackbonePoints   (:, 1)
        GaussWeightsBackbonePoints  (:, 1)
        xiGauss                     (6, :)
        dxiGauss                    (6, :)
        ddxiGauss                   (6, :)
        JxiGauss                    (:, :)
        xiGauss_ds                  (6, :) % Derivative of the strain w.r.t. s
        dxiGauss_ds                 (6, :) % Derivative of the first order time derivative of the strain w.r.t. s
        dJxiGauss                   (:, :) % Derivative w.r.t. s of the Jacobian strain
        ElongationGauss             (:, 1)
        JElongationGauss            (:, :)
        dElongationGauss            (:, 1)
        ddElongationGauss           (:, 1)
        
        % Points along the backbone (curvilinear abscissa) where the quantities must be evaluated
        Points                      (1, :) double
        % Mapping from the points along the backbone to the Gauss points
        PointsToGaussPoints         (:, 2) double
        % Expanded version of PointsToGaussPoints to compute faster some terms
        PointsToGaussPointsExpanded (:, 1) double
        % NUmber of points
        nPoints                     (1, 1) = 0
    end
    
    % Implement the abstract properties inherited from the GVSBody
    properties
        % By default a LVPBackbone has no DOFs
        n    = 0;
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
            Parameters = [RestLength, 0, 0, 0, 0, 0, 0, NGaussPoints];%TODO: This has to be changed to handle cell arrays
            obj = obj@GVSBody(Parameters);

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

            % Initialize the strain and time derivatives along the backbone
            obj.xiGauss                     = repmat(obj.ReferenceStrain, 1, obj.NGaussPointsInt);
            obj.dxiGauss                    = zeros(6, obj.NGaussPointsInt);
            obj.ddxiGauss                   = zeros(6, obj.NGaussPointsInt);
            obj.JxiGauss                    = zeros(6*obj.n, obj.NGaussPointsInt);
            obj.xiGauss_ds                  = zeros(6, obj.NGaussPointsInt);
            obj.dxiGauss_ds                 = zeros(6, obj.NGaussPointsInt);
            obj.dJxiGauss                   = zeros(6*obj.n, obj.NGaussPointsInt);

            obj.ElongationGauss             = zeros(obj.NGaussPointsInt, 1);
            obj.JElongationGauss            = zeros(obj.n, obj.NGaussPointsInt);
            obj.dElongationGauss            = zeros(obj.NGaussPointsInt, 1);
            obj.ddElongationGauss           = zeros(obj.NGaussPointsInt, 1);

            % Assign the points along the backbone and other quantities associated to the points
            obj.setPoints(Points);
            
            % Initialize the Gaussian points and weights for the intervals [0, GaussPoints(i)].
            obj.GaussWeightsBackbonePoints  = zeros(obj.NGaussPointsInt, 1);
            obj.GaussPointsBackbonePoints   = zeros(obj.NGaussPointsInt, 1);
            xPrev = 0;
            for i = 1:obj.NGaussPointsInt
                [xGauss, wGauss] = lgwt(1, xPrev, obj.GaussPoints(i));
                obj.GaussPointsBackbonePoints(i)  = xGauss;
                obj.GaussWeightsBackbonePoints(i) = wGauss;
                xPrev = obj.GaussPoints(i);
            end

            % Preallocate variables of the superclass GVSBody that depend on the number of DOFs
            obj.JEtaGauss           = zeros(6, obj.n, obj.NGaussPointsInt);
            obj.JdrGauss            = zeros(3, obj.n, obj.NGaussPointsInt);
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

            % Store the new points along the backbone and change accordingly the other parameters
            obj.Points                      = Points;
            obj.PointsToGaussPoints         = obj.MapPointsToGaussianPoints(obj.Points);
            obj.PointsToGaussPointsExpanded = repelem(obj.PointsToGaussPoints(:, 1), obj.PointsToGaussPoints(:, 2));
            obj.nPoints                     = length(Points); 

            % Initialize the variables of the backbone
            obj.Strain                          = repmat(obj.ReferenceStrain, [1, obj.nPoints]);
            obj.dStrain                         = zeros(6, obj.nPoints);
            obj.ddStrain                        = zeros(6, obj.nPoints);
            obj.JStrain                         = zeros(6*obj.n, obj.nPoints);
            obj.dJStrain                        = zeros(6*obj.n, obj.nPoints);
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

        % Compute the status of the backbone at the query points (Points)
        % Each element of the rows of sToGaussPoints has the following meaning:
        % 1) A Gaussian point x of the backbone
        % 2) The number of consecutive elements in s that belong the interval starting with x and ending with x+1, where x+1 is the following Gaussian point to x
        % In the future, this method should allow the computation NOT approximated by using the Kinematics method of the GVSBody.
        % This however requires to vectorize the strain methods of the GVSBody function and also the strain functions
        function UpdatePoints(obj)
            % Update the pose at the points
            arguments (Input)
                obj             (1, 1) LVPBackbone
            end
            % Variables initialization
            s               = obj.Points;
            sToGaussPoints  = obj.PointsToGaussPoints;

            % Iterate over the rows of sToGaussPoints
            idxS = 1;
            idxE = 0;
            for i = 1:length(sToGaussPoints)
                %% Get the data for the current interval of Gauss points
                idx     = sToGaussPoints(i, 1);
                idxE    = idxS + sToGaussPoints(i, 2) - 1;
                % Extract all the quantities required for the below computation
                if idx == 0
                    IntvLb      = 0;
                    DeltaL1     = 0;
                    dDeltaL1    = 0;
                    ddDeltaL1   = 0;
                    JL1         = zeros(obj.n, 1);
                    T1          = eye(4);
                    dPose1      = zeros(6, 1);
                    ddPose1     = zeros(6, 1);
                    JPose1      = zeros(6*obj.n, 1);
                    xi1         = obj.ReferenceStrain;
                    Jxi1vec     = zeros(6*obj.n, 1);% Expressed in vectorized format
                    dJxi1vec    = zeros(6*obj.n, 1);% Expressed in vectorized format
                    dxi1        = zeros(6, 1);
                    ddxi1       = zeros(6, 1);
                    xi1_ds      = zeros(6, 1);
                    dxi1_ds     = zeros(6, 1);
                else
                    % Lower bound of the current Gauss interval
                    IntvLb      = obj.GaussPoints(idx);
                    % Elongation at the Gaussian point
                    DeltaL1     = obj.ElongationGauss(idx);
                    dDeltaL1    = obj.dElongationGauss(idx);
                    ddDeltaL1   = obj.ddElongationGauss(idx);
                    JL1         = obj.JElongationGauss(1:obj.n, idx);
                    % Transformation matrix at the Gaussian point
                    T1          = obj.gGauss(1:4, 1:4, idx);
                    dPose1      = obj.EtaGauss(1:6, idx);
                    ddPose1     = obj.dEtaGauss(1:6, idx);
                    JPose1      = reshape(obj.JEtaGauss(1:6, 1:obj.n, idx), 6*obj.n, 1);
                    % Strain at the current Gauss point
                    xi1         = obj.xiGauss(1:6, idx);
                    Jxi1vec     = obj.JxiGauss(1:end, idx);
                    dJxi1vec    = obj.dJxiGauss(1:end, idx);
                    dxi1        = obj.dxiGauss(1:6, idx);
                    ddxi1       = obj.ddxiGauss(1:6, idx);
                    xi1_ds      = obj.xiGauss_ds(1:6, idx);
                    dxi1_ds     = obj.dxiGauss_ds(1:6, idx);
                end
                IntvUb                          = obj.GaussPoints(idx + 1);
                s12                             = s(idxS:idxE); % All points of the backbone in the interval [GaussPoints(idx), GaussPoints(idx+1)]
                % Normalize the points of the interval
                s12                             = (s12-IntvLb)./(IntvUb - IntvLb);
                %% Interpolate the transformation matrices for the pose
                T2                              = obj.gGauss(1:4, 1:4, idx+1);
                dPose2                          = obj.EtaGauss(1:6, idx+1);
                ddPose2                         = obj.dEtaGauss(1:6, idx+1);
                JPose2                          = reshape(obj.JEtaGauss(1:6, 1:obj.n, idx+1), 6*obj.n, 1);
                
                % Perform the interpolation for transformation matrix and store the result
                obj.Pose(1:4, 1:4, idxS:idxE)       = Tinterpolate(T1, T2, s12);
                obj.dPose(1:6, idxS:idxE)           = dPose1.*(1-s12) + dPose2.*s12;
                obj.ddPose(1:6, idxS:idxE)          = ddPose1.*(1-s12) + ddPose2.*s12;
                obj.JPose(1:6, 1:obj.n, idxS:idxE)  = reshape(JPose1.*(1-s12) + JPose2.*s12, 6, obj.n, []);

                %% Interpolate the elongation
                DeltaL2                                     = obj.ElongationGauss(idx+1);
                JL2                                         = obj.JElongationGauss(1:obj.n, idx+1);
                dDeltaL2                                    = obj.dElongationGauss(idx+1);
                ddDeltaL2                                   = obj.ddElongationGauss(idx+1);
                obj.Elongation(idxS:idxE)                   = DeltaL1.*(1-s12) + DeltaL2.*s12;
                obj.JElongation(1:obj.n, idxS:idxE)         = JL1.*(1-s12) + JL2.*s12;
                obj.dElongation(idxS:idxE)                  = dDeltaL1.*(1-s12) + dDeltaL2.*s12;
                obj.ddElongation(idxS:idxE)                 = ddDeltaL1.*(1-s12) + ddDeltaL2.*s12;

                %% Interpolate the strain
                xi2                                     = obj.xiGauss(1:6, idx + 1);
                Jxi2vec                                 = obj.JxiGauss(1:end, idx + 1);
                dJxi2vec                                = obj.dJxiGauss(1:end, idx + 1);
                dxi2                                    = obj.dxiGauss(1:6, idx + 1);
                ddxi2                                   = obj.ddxiGauss(1:6, idx + 1);
                xi2_ds                                  = obj.xiGauss_ds(1:6, idx + 1);
                dxi2_ds                                 = obj.dxiGauss_ds(1:6, idx + 1);
                
                obj.Strain(1:6, idxS:idxE)              = xi1.*(1-s12) + xi2.*s12;
                obj.JStrain(1:6*obj.n, idxS:idxE)       = Jxi1vec.*(1-s12) + Jxi2vec.*s12;
                obj.dJStrain(1:6*obj.n, idxS:idxE)      = dJxi1vec.*(1-s12) + dJxi2vec.*s12;
                obj.dStrain(1:6, idxS:idxE)             = dxi1.*(1-s12) + dxi2.*s12;
                obj.ddStrain(1:6, idxS:idxE)            = ddxi1.*(1-s12) + ddxi2.*s12;
                obj.Strain_ds(1:6, idxS:idxE)           = xi1_ds.*(1-s12) + xi2_ds.*s12;
                obj.dStrain_ds(1:6, idxS:idxE)          = dxi1_ds.*(1-s12) + dxi2_ds.*s12;
                %% Prepare for the next iteration
                idxS                                    = idxE + 1;
            end
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
                        JStrain_                            = reshape(obj.JStrain, 6, obj.n, obj.nPoints);
                        dJStrain_                           = reshape(obj.dJStrain, 6, obj.n, obj.nPoints);
                        Jxi(1:6, 1:idxE, 1:obj.nPoints)     = JStrain_(1:6, 1:idxE, 1:obj.nPoints);
                        dJxi(1:6, 1:idxE, 1:obj.nPoints)    = dJStrain_(1:6, 1:idxE, 1:obj.nPoints);
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
                        JLStrainQ                       = obj.JStrain(6*(1:obj.n), 1:obj.nPoints);
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
            
            % Update the kinematics at the Gauss points using the superclass method
            Update@GVSBody(obj, q, dq, ddq, "EvaluateKinematicTerms", true, ...
                                            "EvaluateInertialTerms" , false, ...
                                            "EvaluateExternalForces", false);

            % Update also the other values associated to the Gaussian points.
            for i = 1:obj.NGaussPointsInt
                [JGaussT, dJGaussT]             = obj.StrainBasis(obj.GaussPoints(i));
                obj.JxiGauss(1:end, i)          = reshape(JGaussT, [], 1);
                obj.dJxiGauss(1:end, i)         = reshape(dJGaussT, [], 1);
                obj.xiGauss(1:6, i)             = JGaussT*q + obj.ReferenceStrain;
                obj.dxiGauss(1:6, i)            = JGaussT*dq;
                obj.ddxiGauss(1:6, i)           = JGaussT*ddq;
                obj.xiGauss_ds(1:6, i)          = dJGaussT*q;
                obj.dxiGauss_ds(1:6, i)         = dJGaussT*dq;

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
            obj.UpdatePoints();
        end
    end
end

