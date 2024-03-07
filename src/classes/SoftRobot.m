classdef SoftRobot < BodyTree
    %SoftRobot extends BodyTree class and offers additional methods
    %for plotting. 
    %TODO: In the future, this class and its methods should be put into a separate MATALB package.
    
    properties
        %Actuators of the soft robot organized in a cell array.
        Actuators;
        %Total number of actuators in the kinematic tree.
        N_A = 0;
        
        % Graphical properties
        SegmentRadius
        Color
        SegmentBaseColor
    end

    properties (Constant)
        % Maximum number of bodies and joints, required for code generation.
        MaxActuatorsNumber = 50;
    end

    properties (Constant, Access = private)
        DefaultColor = [0 160 219]./256;
    end
    
    %Private properties
    properties(Access = public)
        % Gaussian points for the actuators
        ActuatorGaussianPoints;
        % Gaussian weights for the actuators
        ActuatorGaussianWeights;
        % Matrix that stores the interval defining each body
        BodiesInterval = [0, 0];
        % Vector containing the number of Gaussian points for each actuator
        ActuatorGaussianPointLength = 0;
    end
    
    methods
        %Class constructor.
        %TODO: Improve the input parsing process. Now, it can lead to
        %several errors.
        %Note that we can use the input parser because this method is never
        %called during code generation in simulink.
        function obj = SoftRobot(Joints, Bodies, Actuators)
            if nargin <= 2
                Actuators = {};
            end

            % Build the body tree
            obj = obj@BodyTree(Joints, Bodies);

            % Store the rest length of each body
            if coder.target("MATLAB")
                N_B = obj.N_B;
            else
                N_B = BodyTree.MaxBodiesNumber;
            end
            obj.BodiesInterval = zeros(N_B, 2);
            
            for i = 1:N_B
                if isnumeric(obj.Bodies{i})
                    continue;
                end
                if i == 1
                    obj.BodiesInterval(i, 1:2) = [0, obj.Bodies{i}.RestLength];
                else
                    obj.BodiesInterval(i, 1:2) = [obj.BodiesInterval(i-1, 2), obj.BodiesInterval(i-1, 2) + obj.Bodies{i}.RestLength];
                end
            end

            % Compute the number of actuators
            obj.N_A = 0;
            l_B     = length(Actuators);
            for i = 1:SoftRobot.MaxActuatorsNumber
                if i <= l_B
                    if ~isnumeric(Actuators{i})
                        obj.N_A = obj.N_A + 1;
                    end
                end
            end
            % Assign the actuators
            obj.Actuators = cell(SoftRobot.MaxActuatorsNumber, 1);
            obj.Actuators = Actuators;
            
            if coder.target("MATLAB")
                N_A = obj.N_A;
            else
                N_A = SoftRobot.MaxActuatorsNumber;
            end

            % Allocate the Gaussian points and weights for each actuator. Because of the discontinuity of the strain, the integrals are expanded for actuators that span more than one body.
            obj.ActuatorGaussianPoints      = cell(N_A, 1);
            obj.ActuatorGaussianWeights     = cell(N_A, 1);
            obj.ActuatorGaussianPointLength = zeros(N_A, 1);
            % Get the rest length of the robot
            RobotLength = obj.BodiesInterval(obj.N_B, 2);
            for i = 1:N_A
                if isnumeric(Actuators{i})%Needed for code generation
                    obj.ActuatorGaussianPoints{i}  = 0;
                    obj.ActuatorGaussianWeights{i} = 0;
                    continue;
                end
                % Get the index of the body where the actuator starts
                sStart   = Actuators{i}.sStart;
                if sStart == RobotLength
                    idxBodyStart = obj.N_B;
                elseif sStart == 0
                    idxBodyStart = 1;
                else
                    idxBodyStart = find(sStart >= obj.BodiesInterval(1:obj.N_B, 1) & sStart <= obj.BodiesInterval(1:obj.N_B, 2), 1);
                end
                % Get the index of the body where the actuator ends
                sEnd     = Actuators{i}.sEnd;
                if sEnd == RobotLength
                    idxBodyEnd = obj.N_B;
                elseif sEnd == 0
                    idxBodyEnd = 1;
                else
                    idxBodyEnd = find(sEnd >= obj.BodiesInterval(1:obj.N_B, 1) & sEnd <= obj.BodiesInterval(1:obj.N_B, 2), 1);
                end
                % Compute the number of bodies trough which the actuator extends
                idxBodyInterval = idxBodyEnd(1) - idxBodyStart(1) + 1;
                
                % Preallocate the Gauss points and weights
                NGauss = obj.Actuators{i}.NGaussPoints;
                obj.ActuatorGaussianPointLength(i) = idxBodyInterval*NGauss;
                obj.ActuatorGaussianPoints{i}  = zeros(idxBodyInterval*NGauss, 1);
                obj.ActuatorGaussianWeights{i} = zeros(idxBodyInterval*NGauss, 1);
                % Assign the Gauss points and weights
                k = 1;
                for j = idxBodyStart(1):idxBodyEnd(1)
                    % Check if we are at the start body of the actuator
                    if j == idxBodyStart(1)
                        [obj.ActuatorGaussianPoints{i}((k-1)*NGauss+1:k*NGauss), obj.ActuatorGaussianWeights{i}((k-1)*NGauss+1:k*NGauss)] = lgwt(NGauss, obj.Actuators{i}.sStart, obj.BodiesInterval(j, 2));
                    elseif j == idxBodyEnd(1) % Check if we are at the end body of the actuator
                        [obj.ActuatorGaussianPoints{i}((k-1)*NGauss+1:k*NGauss), obj.ActuatorGaussianWeights{i}((k-1)*NGauss+1:k*NGauss)] = lgwt(NGauss, obj.BodiesInterval(j, 1), obj.Actuators{i}.sEnd);
                    else % Check if we are in a body through which the actuator passes completely
                        [obj.ActuatorGaussianPoints{i}((k-1)*NGauss+1:k*NGauss), obj.ActuatorGaussianWeights{i}((k-1)*NGauss+1:k*NGauss)] = lgwt(NGauss, obj.BodiesInterval(j, 1), obj.BodiesInterval(j, 2));
                    end
                    % Update iterationv variables
                    k = k + 1;
                end
            end
            
            % Store the radius of each body for plotting purposes
            SegmentRadius = cell(N_B, 1);
            for i = 1:N_B
                if isnumeric(Bodies{i})% Required for code generation
                    SegmentRadius{i} = [0, 0];
                    continue;
                end
                SegmentRadius{i, 1} = [Bodies{i}.Parameters(2), Bodies{i}.Parameters(3)];
            end
            obj.SegmentRadius = SegmentRadius;
            
            %TODO: These could be optional arguments
            obj.Color            = SoftRobot.DefaultColor;
            obj.SegmentBaseColor = repmat(obj.Color, N_B, 1);
        end

        function T = DirectKinematics(obj, q, points)
            %Evaluate the direct kinematics at specific points along the robot
            %backbone
            %Args:
            %   q   ([double], [sym]): Configuration variables
            %   points     ([double]): Ordered array of lengths indicating the backbone points for which the direct kinematics has to be evaluated
            %Return:
            %   ([double], [sym]): Homogeneous transformation matrices for the backbone points specified by idx.   

            switch nargin
                case 2
                    % If the points are not specified, evaluate the DK using the BodyTree method
                    T = DirectKinematics@BodyTree(obj, q);
                    return;
            end

            % Check that the input vectors are in column format.
            if ~iscolumn(q)
                q = q';
            end
            
            % Useful variables
            pointsLength = length(points);

            % For each point, find the index of the corresponding body in
            % the chain and the remapped curvilinear abscissa
            [idx, pointsBody] = obj.getBodyIndexFromBackboneDistance(points);

            % Remove repeated occurences of the bodies indexes
            [idxUnique, ~, idxUniqueOcc] = unique(idx);

            % T is matrix of vertically stacked 4x4 transformation matrices
            T   = repmat(eye(4, 'like', q), pointsLength, 1);

            % Compute the DK for all the precedeing bodies to which the points belong
            TPrevBodies = DirectKinematics@BodyTree(obj, q, idxUnique-1);

            % Retreive the index vector for the configuration of the bodies
            idxBody     = obj.getBodyConfigurationIndex(idx);
            
            % Compute the output
            T_i = zeros(4, 4, "like", q);% Output preallocation for code generation
            % Iterate over all the bodies because of code generation
            j = 1;
            for i = 1:obj.MaxBodiesNumber
                if i <= obj.N_B
                    if isnumeric(obj.Bodies{i})
                        continue;
                    end
                    % Check if the current body is the starting index
                    if i == idx(j)
                        for k = j:pointsLength
                            % Check if the current index corresponds to a
                            % different body and in case stop the loop
                            if i ~= idx(k)
                                % Update the j index and continue the iteration over
                                % the remaining bodies
                                j = k;
                                break;
                            end
                            % Evaluate the direct kinematics at the query point
                            T_i(1:4, 1:4)         = obj.Bodies{i}.T_s(q(idxBody(k, 1):idxBody(k, 2), 1), pointsBody(k));
                            % Update the output
                            T(4*(k-1)+1:4*k, 1:4) = TPrevBodies(4*(idxUniqueOcc(k)-1)+1:4*idxUniqueOcc(k), 1:4)*T_i;% Use the same transform of the preceideing body for the same occurrences
                        end
                    end
                end
            end
        end
        
        
        
        function J = BodyJacobian(obj, q, points)
            %Evaluate the body Jacobian at specified points along the backbone. If i is not specified, the method returns the body Jacobian of each body of the chain.
            %
            %Args:
            %   q   ([double], [sym]): Configuration variables
            %   points     ([double]): Array of backbone points indicating where the body jacobian has to be evaluated
            %Return:
            %   {[double], [sym]}: length(points)*6 x n body Jacobian with
            %   angular and linear components for each point specified by
            %   points

            switch nargin
                case 1
                    q   = zeros(obj.n, 1, 'like', q);
                    J   = BodyJacobian@BodyTree(obj, q);
                    return;
                case 2
                    J   = BodyJacobian@BodyTree(obj, q);
                    return;
            end

            % Check that the input vectors are in column format.
            if ~iscolumn(q)
                q = q';
            end

            % Store useful variables
            pointsLength = length(points);

            % Preallocate the output for code generation
            J = zeros(6*pointsLength, obj.n, "like", q);

            % For each point, find the index of the corresponding body in
            % the chain and the remapped curvilinear abscissa in the body
            [idx, sBody] = obj.getBodyIndexFromBackboneDistance(points);

            % Remove repeated occurences of the bodies indexes
            [idxUnique, ~, idxUniqueOcc] = unique(idx);

            % Get the indexes of the configuration variables corresponding
            % to the required bodies
            idxQBody = obj.getBodyConfigurationIndex(idx);

            % Evalute the body jacobian for the previous bodies 
            J_prev = BodyJacobian@BodyTree(obj, q, idxUnique-1);
            % Define useful variables for iteration and code generation
            % support
            J_i          = zeros(6, obj.n, "like", q);
            q_idx_start  = 0;
            j            = 1;
            % Update the Jacobian considering the current body
            % Iteration over all the bodies is required to allow for code
            % generation
            for i = 1:obj.MaxBodiesNumber
                if i <= obj.N_B
                    if isnumeric(obj.Bodies{i})
                        continue;
                    end

                    % Check if the current body is the starting index
                    if i == idx(j)
                        for k = j:pointsLength
                            % Check if the current index corresponds to a
                            % different body and in case stop the loop
                            if i ~= idx(k)
                                % Update the j index and continue the iteration over
                                % the remaining bodies
                                j = k;
                                break;
                            end

                            % Prepare the variables for the iteration
                            q_idx_start = idxQBody(k, 1);
                            q_idx_end   = idxQBody(k, 2);
            
                            % Compute the terms associated to the body
                            [J_omega, J_v, T]   = obj.Bodies{i}.BodyJacobian(q(q_idx_start:q_idx_end, 1), sBody(k));
                            R_i_T               = T(1:3, 1:3)';
                            t_i                 = T(1:3, 4);
            
                            % Initialize the Jacobian value using the Jacobian of the
                            % previous bodies
                            J_i = J_prev(1+6*(idxUniqueOcc(k)-1):6*idxUniqueOcc(k), 1:obj.n);
                            
                            % Update the value of the Jacobian for the bodies previous
                            % to the current one by rotating them in the body frame of
                            % the current body
                            if q_idx_start ~= 1
                                J_i(4:6, 1:q_idx_start-1) = real(R_i_T*(J_i(4:6, 1:q_idx_start-1) - skew(t_i)*J_i(1:3, 1:q_idx_start-1)));
                                J_i(1:3, 1:q_idx_start-1) = real(R_i_T*J_i(1:3, 1:q_idx_start-1));
                            end
                            J_i(1:6, q_idx_start:q_idx_end) = [J_omega; J_v];
            
                            % Assign J_i to the output
                            J(1 + 6*(k-1):6*k, 1:obj.n) = J_i(1:6, 1:obj.n);
                        end
                    end
                end
            end

        end
        
        
        function [q, converged, e] = InverseKinematics(obj, T, points, q0, N, task_flags, AngularErrorThsd, LinearErrorThsd)
            %Evaluate the inverse kinematics numerically using a Newton-Rapson iteration scheme.
            %
            %Args:
            %   T   ([double, double])      : Target transformation matrices vertically stacked
            %   q0  ([double, double])      : Initial guess, the default value is q0 = zeros(n, 1)
            %   points ([double])           : Vector of points along the robot backbone for which T has to be found
            %   N   (double)                : Maximum number of iterations
            %   task_flags ([double, bool]) : Vector of flags specifying for each indexed body what components of the task vector should be considered
            %Return:
            %   {[double], [sym]}: Homogeneous transformation matrices for each body.

            % Default values
            DefaultN                    = 4;   %Number of Newton iterations
            DefaultAngularErrorThsd     = 1e-3;%Default threshold in the Newton scheme for the angular position
            DefaultLinearErrorThsd      = 1e-2;%Default threshold in the Newton scheme for the linear position
            
            switch nargin
                case 2
                    points      = linspace(1, obj.N_B, obj.N_B);
                    q0          = zeros(obj.n, 1);
                    N           = DefaultN;
                    task_flags  = ones(obj.N_B*6, 1);
                    AngularErrorThsd = DefaultAngularErrorThsd;
                    LinearErrorThsd  = DefaultLinearErrorThsd;
                case 3
                    q0  = zeros(obj.n, 1);
                    N   = DefaultN;
                    task_flags  = ones(obj.N_B*6, 1);
                    AngularErrorThsd = DefaultAngularErrorThsd;
                    LinearErrorThsd  = DefaultLinearErrorThsd;
                case 4
                    N   = DefaultN;
                    task_flags  = ones(obj.N_B*6, 1);
                    AngularErrorThsd = DefaultAngularErrorThsd;
                    LinearErrorThsd  = DefaultLinearErrorThsd;
                case 5
                    task_flags  = ones(obj.N_B*6, 1);
                    AngularErrorThsd = DefaultAngularErrorThsd;
                    LinearErrorThsd  = DefaultLinearErrorThsd;
                case 6
                    AngularErrorThsd = DefaultAngularErrorThsd;
                    LinearErrorThsd  = DefaultLinearErrorThsd;
                case 7
                    LinearErrorThsd  = DefaultLinearErrorThsd;
            end

            % Call the superclass method that still works because of the
            % overloading of the methods DirectKinematics and BodyJacobian
            % of this class
            [q, converged, e] = InverseKinematics@BodyTree(obj, T, points, q0, N, task_flags, AngularErrorThsd, LinearErrorThsd);
        end
        
        % Evaluate the strain at a given point along the robot structure and in the given configuration
        function [xi_, J_xi_] = xi(obj, q, s)
            %
            
            % Get the rest length of the robot
            RobotLength = obj.BodiesInterval(obj.N_B, 2);
            
            % Preallocate the output
            xi_   = zeros(6, 1); 
            xi_(6)= 1;
            J_xi_ = zeros(6, obj.n);

            % Return 0 if s is outside the allowed range
            if s > RobotLength || s < 0
                return;
            end
            % % Find the index of the body to which the point belongs
            % if s == RobotLength
            %     idx = obj.N_B;
            % elseif s == 0
            %     idx = 1;
            % else
            %     idx = find(s >= obj.BodiesInterval(1:obj.N_B, 1) & s <= obj.BodiesInterval(1:obj.N_B, 2), 1);
            % end
            % idxBody = idx(1);
            % 
            % % Remap s into the body interval
            % s = s-obj.BodiesInterval(idxBody, 1);
            [idx, s] = obj.getBodyIndexFromBackboneDistance(s);
            idxBody = idx(1);

            % Get the configuration of the body
            idxQ  = obj.getBodyConfigurationIndex(idxBody);
            qBody = q(idxQ(1):idxQ(2));
            % Compute the strain and its Jacobian
            xi_                         = obj.Bodies{idxBody}.xi(qBody, s);
            J_xi_(1:6, idxQ(1):idxQ(2)) = obj.Bodies{idxBody}.Jxi(qBody, s);

            % idxQ = 1;
            % for i = 1:BodyTree.MaxBodiesNumber
            %     if isnumeric(obj.Bodies{i})
            %         continue;
            %     end
            %     idxQ = idxQ + obj.Joints{i}.n;
            %     if i == idxBody
            %         qBody = q(idxQ:idxQ+obj.Bodies{i}.n-1);
            % 
            %         % Compute the strain and its Jacobian
            %         xi_     = obj.Bodies{i}.xi(qBody, s);
            %         J_xi_(1:6, idxQ:idxQ+obj.Bodies{i}.n-1) = obj.Bodies{i}.Jxi(qBody, s);
            %         break;
            %     end
            %     idxQ = idxQ + obj.Bodies{i}.n;
            % end
        end

        % Computation of the actuation matrix and the actuator elongation
        function [A, y] = ActuationMatrix(obj, q)
            % Preallocate the output
            A = zeros(obj.n, obj.N_A);
            y = zeros(obj.N_A, 1);
            % Iterate over the actuators
            l_A = length(obj.Actuators);
            for i = 1:SoftRobot.MaxActuatorsNumber
                if i <= l_A
                    % Check that we have an actuator, used for code generation
                    if isnumeric(obj.Actuators{i})
                        continue;
                    end
                    % Integrate numerically 
                    for j = 1:obj.ActuatorGaussianPointLength(i)
                        % Retrieve the Gaussian point
                        sGauss = obj.ActuatorGaussianPoints{i}(j);
                        % Compute the distance of the actuator from the backbone and its spatial derivative
                        d   = obj.Actuators{i}.ActuatorDistance(sGauss);
                        dds = obj.Actuators{i}.dActuatorDistance(sGauss);
                        % Compute the strain and its Jacobian w.r.t. q
                        [xi, Jxi]  = obj.xi(q, sGauss);
                        % Compute the unit tangent vector to the actuator
                        ut  = skew(xi(1:3))*d + xi(4:6) + dds;
                        t   = ut/norm(ut);
                        % Compute the actuation basis
                        Phi_a = [skew(d)*t; t];
                        % Update the actuation matrix
                        A(1:obj.n, i) = A(1:obj.n, i) + obj.ActuatorGaussianWeights{i}(j)*(Jxi'*Phi_a); 
                        % Compute the reference strain
                        [xi_ref, ~] = obj.xi(zeros(obj.n, 1), sGauss);
                        % Compute the reference value for the actuation basis
                        ut_ref    = skew(xi_ref(1:3))*d + xi_ref(4:6) + dds;
                        t_ref     = ut_ref/norm(ut_ref);
                        Phi_a_ref = [skew(d)*t_ref; t_ref];
                        % Update the actuator elongation
                        y(i) = y(i) + obj.ActuatorGaussianWeights{i}(j)*((Phi_a-Phi_a_ref)'*(xi_ref + [zeros(3, 1); dds]));
                    end
                end
            end
           
            % Add to y the contribution of the current configuration
            y = y + A'*q;

        end

        
        %Plot the robot in the given configuration.
        function [ptch] = plot(obj, q, varargin)
            %Parse the optional paramters
            p = inputParser;
            validScalar = @(x) isnumeric(x);
            addRequired(p,  'q');
            addParameter(p, 'FaceAlpha', 1, validScalar);
            addParameter(p, 'Color', obj.Color);
            addParameter(p, 'SegmentBaseColor', obj.SegmentBaseColor);
            addParameter(p, 'AmbientStrength', 0.3, validScalar);
            addParameter(p, 'SegmentBaseColorExtension', 0);%Number of cylinders used for the base
            parse(p, q, varargin{:});

            %Compute the mesh of the robot in the current configuration
            [verts, faces, faces_color] = obj.robotGeometry(q, p.Results.Color, p.Results.SegmentBaseColor, p.Results.SegmentBaseColorExtension);
            ptch = patch('Faces',faces, 'Vertices',verts, ...
                  'FaceVertexCData', faces_color, ...
                  'FaceColor', 'flat', ... %Set to 'interp' if we provide the color of the vertexes instead of that of the faces
                  'LineStyle', 'none', ...
                  'AmbientStrength', p.Results.AmbientStrength, ...
                  'FaceAlpha', p.Results.FaceAlpha);
        end

        function plotBackbone(obj, q, varargin)
            %Parse the optional paramters
            p = inputParser;
            addRequired(p,  'q');
            addParameter(p, 'Color', obj.Color);
            addParameter(p, 'LineWidth', 1);
            addParameter(p, 'Frame', true);
            addParameter(p, 'FrameScale', 0.03);
            addParameter(p, 'FrameSkip', 10);
            parse(p, q, varargin{:});
            %Compute the backbone of the robot
            [x, y, z, R] = obj.robotBackbone(q);
            plot3(x, y, z, 'Color', p.Results.Color, 'LineWidth', p.Results.LineWidth);
            hold on;
            %Plot also the reference frames
            if p.Results.Frame
                for i = 1:p.Results.FrameSkip:size(R, 3)
                    n1 = norm(R(:, 1))/p.Results.FrameScale;
                    n2 = norm(R(:, 2))/p.Results.FrameScale;
                    n3 = norm(R(:, 3))/p.Results.FrameScale;
                    quiver3(x(i), y(i), z(i), R(1, 1, i)/n1, R(2, 1, i)/n1, R(3, 1, i)/n1, 'Color', 'r');
                    quiver3(x(i), y(i), z(i), R(1, 2, i)/n2, R(2, 2, i)/n2, R(3, 2, i)/n2, 'Color', 'g');
                    quiver3(x(i), y(i), z(i), R(1, 3, i)/n3, R(2, 3, i)/n3, R(3, 3, i)/n3, 'Color', 'b');
                end
            end
        end
    end

    methods(Access = private)

        % Get the index of the body for the point along the backbone and
        % also the remapped abscissa in the body frame
        function [idx, pointsBody]= getBodyIndexFromBackboneDistance(obj, points)
            % Get the rest length of the robot
            RobotLength = obj.BodiesInterval(obj.N_B, 2);
            % Store the length of the input
            pointsLength = length(points);

            % For each point, find the index of the corresponding body in the chain
            idx         = zeros(pointsLength, 1);% Vector of indexes
            pointsBody  = zeros(pointsLength, 1);% Vector containg points remapped into their body interval
            for i = 1:pointsLength
                if points(i) == RobotLength
                    idx(i) = obj.N_B;
                elseif points(i) == 0
                    idx(i) = 1;
                else
                    idx(i) = find(points(i) >= obj.BodiesInterval(1:obj.N_B, 1) & points(i) <= obj.BodiesInterval(1:obj.N_B, 2), 1);
                end
                
                % Remap the curvilinear abscissa into the body interval
                pointsBody(i) = points(i)-obj.BodiesInterval(idx(i), 1);
            end
        end


        %Compute the position of the robot backbone
        function [x, y, z, R] = robotBackbone(obj, q)
            %Function hyperparameters increase to have a better representation
            %Number of disks along each segment
            s_step   = 100;
            
            %Preallocate vertex and face matrices
            %The total number of vertexes includes also the origin and the
            %tip of the backbone (+2)
            n_points = obj.N_B*s_step;
            verts = zeros(3, n_points);
            R     = zeros(3, 3, n_points);
        
            s  = linspace(0, 1, s_step);

            T_prev = obj.T0;
            idx_verts = 1;
            idx_q     = 1;
            %Iterate over the jointa and bodies
            for i = 1:obj.N_B
                if obj.Joints{i}.n ~= 0
                    %Get the configuration variables of the joint
                    q_i      = q(idx_q:idx_q + obj.Joints{i}.n - 1);
                    obj.Joints{i}.Update(q_i, zeros(size(q_i)), zeros(size(q_i)));
                else
                    obj.Joints{i}.Update([], [], []);
                end
                %Update the transformation matrix
                T_prev   = T_prev*obj.Joints{i}.T_;
                %Update the index of the configuration variables
                idx_q    = idx_q + obj.Joints{i}.n;
                %Get the configuration variables of the current body
                q_i      = q(idx_q:idx_q + obj.Bodies{i}.n - 1);
                for j=1:s_step
                    %Get the transformation matrix from current frame to world
                    if s(j) == 0%We use the limit expression at the base to avoid numerical problems
                        T_s_j = eye(4);
                    else
                        T_s_j = obj.Bodies{i}.T_s(q_i, s(j)*obj.Bodies{i}.RestLength);
                    end
                    T_ = T_prev*T_s_j;
                    verts(:, idx_verts) = T_(1:3, 4);
                    R(:, :, idx_verts) = T_(1:3, 1:3);
                    idx_verts = idx_verts + 1;
                end
                %Prepare for the next iteration
                T_prev = T_;
                %Update the index of the configuration variables
                idx_q    = idx_q + obj.Bodies{i}.n;
            end
            x = verts(1, :);
            y = verts(2, :);
            z = verts(3, :);
        end

        function [verts, faces, faces_color] = robotGeometry(obj, q, color, segmentBaseColor, segmentBaseColorExtension)
            %Function hyperparameters increase to have a better representation
            %Number of disks along each segment
            s_step   = 100;
            %Number of points representing each cross sectional area (disk)
            phi_step = 100;
        
            %Preallocate vertex and face matrices
            %The total number of vertexes includes also the origin and the
            %tip of the backbone (+2)
            n_verts = obj.N_B*s_step*phi_step + 2;
            n_faces = obj.N_B*s_step*phi_step + phi_step;
            
            verts = zeros(n_verts, 3);
            faces = zeros(n_faces, 4);
            %Assign the colors of the segment bases
            faces_color = repmat(color, n_faces, 1);
            for i = 1:obj.N_B
                color_intv = (i-1)*(s_step*phi_step)+1:(i-1)*(s_step*phi_step) + segmentBaseColorExtension*phi_step;
                faces_color(color_intv, :) = repmat(segmentBaseColor(i, :), phi_step*segmentBaseColorExtension, 1);
            end
        
            s  = linspace(0, 1, s_step);
            phi= linspace(0, 2*pi, phi_step);

            T_prev = obj.T0;
            %Body origin
            verts(n_verts-1, :) = T_prev(1:3, 4)';

            idx_verts = 1;
            idx_q     = 1;
            %Iterate over the joints and bodies
            for i = 1:obj.N_B
                if obj.Joints{i}.n ~= 0
                    %Get the configuration variables of the joint
                    q_i      = q(idx_q:idx_q + obj.Joints{i}.n - 1);
                    obj.Joints{i}.Update(q_i, zeros(size(q_i)), zeros(size(q_i)));
                else
                    obj.Joints{i}.Update([], [], []);
                end
                %Update the transformation matrix
                T_prev   = T_prev*obj.Joints{i}.T_;
                %Update the index of the configuration variables
                idx_q    = idx_q + obj.Joints{i}.n;
                %Compute the position for the body points
                r_base   = obj.SegmentRadius{i}(1);
                r_tip    = obj.SegmentRadius{i}(2);
                %Get the configuration variables of the current body
                q_i      = q(idx_q:idx_q + obj.Bodies{i}.n - 1);
                for j=1:s_step
                    
                    %Get the transformation matrix from current frame to world
                    if s(j) == 0%We use the limit expression at the base to avoid numerical problems
                        T_s_j = eye(4);
                    else
                        T_s_j = obj.Bodies{i}.T_s(q_i, s(j)*obj.Bodies{i}.RestLength);
                    end
                    T_ = T_prev*T_s_j;
                    %Linearly interpolate the radius of the segment base
                    %and tip.
                    %rho = (1-s(j))*r_base + s(j)*r_tip;
                    for k = 1:phi_step
                        try
                            rho = obj.Bodies{i}.Radius(s(j)*obj.Bodies{i}.RestLength, phi(k), q_i);
                        catch
                            rho = obj.Bodies{i}.Radius(s(j)*obj.Bodies{i}.RestLength);
                        end
                        vertex = T_*[rho*cos(phi(k));...
                                     rho*sin(phi(k));...
                                                  0;...
                                                  1];
                        verts(idx_verts, :) = vertex(1:3)';
                        idx_verts = idx_verts + 1;
                    end
                end
                %Prepare for the next iteration
                T_prev = T_prev*obj.Bodies{i}.T(q_i);
                %Update the index of the configuration variables
                idx_q    = idx_q + obj.Bodies{i}.n;
            end

            
            %Close the base layer
            idx_verts = 1;
            for idx_faces = 1:phi_step
                if mod(idx_verts, phi_step) == 0
                    faces(idx_faces, :) = [idx_verts, n_verts-1, n_verts-1, idx_verts-phi_step+1];
                else
                    faces(idx_faces, :) = [idx_verts, n_verts-1, n_verts-1, idx_verts + 1];
                end
                idx_verts = idx_verts + 1;
            end
            
            %Reset the counter for the vertices
            idx_verts = 1;
        
            %Build the faces of the mesh
            for idx_faces = idx_faces+1:n_faces-phi_step
                if mod(idx_verts, phi_step) == 0
                    faces(idx_faces, :) = idx_verts + [0 1 phi_step -phi_step+1];
                else
                    faces(idx_faces, :) = idx_verts + [0 phi_step phi_step+1 1];
                end
                idx_verts = idx_verts + 1;
            end
        
            %Close the last layer
            %T_prev stores the position of the last point on the backbone (tip)
            verts(n_verts, :) = T_prev(1:3, 4)';
            for idx_faces = idx_faces+1:n_faces
                if mod(idx_verts, phi_step) == 0
                    faces(idx_faces, :) = [idx_verts, n_verts, n_verts, idx_verts-phi_step+1];
                else
                    faces(idx_faces, :) = [idx_verts, n_verts, n_verts, idx_verts + 1];
                end
                idx_verts = idx_verts + 1;
            end
        end
    
    end
end