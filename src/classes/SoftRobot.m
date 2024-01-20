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
        MaxActuatorsNumber = 40;
    end

    properties (Constant, Access = private)
        DefaultColor = [0 160 219]./256;
    end

    properties(Access = public)
        % Gaussian points for the actuators
        ActuatorGaussianPoints;
        % Gaussian weights for the actuators
        ActuatorGaussianWeights;
        % Vector that stores the interval defining each body
        BodiesInterval = 0;
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
            obj.BodiesInterval = zeros(obj.N_B + 1, 1);
            for i = 1:obj.N_B
                obj.BodiesInterval(i+1) = obj.BodiesInterval(i) + obj.Bodies{i}.RestLength;
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

            % Allocate the Gaussian points and weights for each actuator
            obj.ActuatorGaussianPoints  = cell(obj.N_A, 1);
            obj.ActuatorGaussianWeights = cell(obj.N_A, 1);
            for i = 1:obj.N_A
                [obj.ActuatorGaussianPoints{i}, obj.ActuatorGaussianWeights{i}] = lgwt(obj.Actuators{i}.NGaussPoints, obj.Actuators{i}.sStart, obj.Actuators{i}.sEnd);
            end
            
            % Store the radius of each body for plotting purposes
            N_B = obj.N_B;
            SegmentRadius = cell(N_B, 1);
            for i = 1:N_B
                if isnumeric(Bodies{i})
                    continue;
                end
                SegmentRadius{i, 1} = [Bodies{i}.Parameters(2), Bodies{i}.Parameters(3)];
            end
            obj.SegmentRadius = SegmentRadius;
            
            %TODO: These could be optional arguments
            obj.Color            = SoftRobot.DefaultColor;
            obj.SegmentBaseColor = repmat(obj.Color, obj.N_B, 1);
        end
        
        % Evaluate the strain at a given point along the robot structure and in the given configuration
        function [xi_, J_xi_] = xi(obj, q, s)
            %
            
            % Get the rest length of the robot
            RobotLength = obj.BodiesInterval(end);
            
            % Preallocate the output
            xi_   = zeros(6, 1); 
            xi_(6)= 1;
            J_xi_ = zeros(6, obj.n);

            % Return 0 if s is outside the allowed range
            if s > RobotLength || s < 0
                return;
            end

            % Find the index of the body to which the point belongs
            if s == RobotLength
                idxBody = obj.N_B;
            elseif s == 0
                idxBody = 1;
            else
                idxBody = find(s >= obj.BodiesInterval, 1, 'last');
            end

            % Remap s into the body interval
            s = s-obj.BodiesInterval(idxBody);

            % Get the configuration of the body
            idxQ = 1;
            for i = 1:idxBody-1
                idxQ = idxQ + obj.Joints{i}.n + obj.Bodies{i}.n;
            end
            idxQ  = idxQ + obj.Joints{idxBody}.n;
            qBody = q(idxQ:idxQ+obj.Bodies{idxBody}.n-1);

            % Compute the strain and its Jacobian
            xi_     = obj.Bodies{idxBody}.xi(qBody, s);
            J_xi_(1:6, idxQ:idxQ+obj.Bodies{idxBody}.n-1) = obj.Bodies{idxBody}.Jxi(qBody, s);
        end

        % Computation of the actuation matrix and the actuator elongation
        function [A, y] = ActuationMatrix(obj, q)
            % Preallocate the output
            A = zeros(obj.n, obj.N_A);
            y = zeros(obj.N_A, 1);
            % Iterate over the actuators
            for i = 1:obj.N_A
                % Check that we have an actuator, used for code generation
                if isnumeric(obj.Actuators{i})
                    continue;
                end
                % Integrate numerically 
                for j = 1:obj.Actuators{i}.NGaussPoints
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