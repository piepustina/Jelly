classdef SoftRobot < BodyTree
    %SoftRobot extends BodyTree class and offers additional methods
    %for plotting. 
    %SOFTROBOT Class that represents a soft robot (for the moment only to make plots)
    
    properties (Constant, Access = private)
        DefaultColor = [0 160 219]./256;
    end

    properties
        SegmentRadius
        Color
        SegmentBaseColor
    end
    
    methods
        %Class constructor.
        %TODO: Improve the input parsing process. Now, it can lead to
        %several errors.
        %Note that we can use the input parser because this method is never
        %called during code generation in simulink.
        function obj = SoftRobot(Joints, Bodies, SegmentRadius, T0)
            if nargin <= 3
                N_B = length(Bodies);
                SegmentRadius = cell(N_B, 1);
                for i = 1:N_B
                    SegmentRadius{i, 1} = [Bodies{i}.Parameters(2), Bodies{i}.Parameters(2)];
                end
            end
            if nargin <= 4
                T0 = eye(4);
            end
            obj = obj@BodyTree(Joints, Bodies)
            
            %Process the input
            obj.SegmentRadius = SegmentRadius;
            obj.T0 = T0;

            %These could be optional arguments
            obj.Color = SoftRobot.DefaultColor;
            obj.SegmentBaseColor = repmat(obj.Color, obj.N_B, 1);
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
            
            %Approximate the configuration vector to run the limits
            q = obj.ApproxQ(q);

            T_prev = obj.T0;
            idx_verts = 1;
            idx_q     = 1;
            %Iterate over the bodies
            for i = 1:obj.N_B
                %Get the configuration variables of the current body
                q_i      = q(idx_q:idx_q + obj.Bodies{i}.n - 1);
                for j=1:s_step
                    %Get the transformation matrix from current frame to world
                    if s(j) == 0%We use the limit expression at the base to avoid numerical problems
                        T_s_j = eye(4);
                    else
                        T_s_j = obj.Joints{i}.T_s(q_i, obj.Bodies{i}.Parameters, s(j)*obj.Bodies{i}.RestLength);
                    end
                    T_ = T_prev*T_s_j;
                    verts(:, idx_verts) = T_(1:3, 4);
                    R(:, :, idx_verts) = T_(1:3, 1:3);
                    idx_verts = idx_verts + 1;
                end
                %Prepare for the next iteration
                %T_prev = T_prev*obj.Joints{i}.T(q_i, obj.Joints{i}.Parameters);
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
            s_step   = 40;
            %Number of points representing each cross sectional area (disk)
            phi_step = 40;
        
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
            
            %Approximate the configuration vector to run the limits
            q = obj.ApproxQ(q);

            T_prev = obj.T0;
            %Body origin
            verts(n_verts-1, :) = T_prev(1:3, 4)';

            idx_verts = 1;
            idx_q     = 1;
            %Iterate over the bodies
            for i = 1:obj.N_B
                r_base   = obj.SegmentRadius{i}(1);
                r_tip    = obj.SegmentRadius{i}(2);
                %Get the configuration variables of the current body
                q_i      = q(idx_q:idx_q + obj.Bodies{i}.n - 1);
                for j=1:s_step
                    
                    %Get the transformation matrix from current frame to world
                    if s(j) == 0%We use the limit expression at the base to avoid numerical problems
                        T_s_j = eye(4);
                    else
                        T_s_j = obj.Joints{i}.T_s(q_i, obj.Bodies{i}.Parameters, s(j)*obj.Bodies{i}.RestLength);
                    end
                    T_ = T_prev*T_s_j;
                    %Linearly interpolate the radius of the segment base
                    %and tip.
                    rho = (1-s(j))*r_base + s(j)*r_tip;
                    for k = 1:phi_step
                        vertex = T_*[rho*cos(phi(k));...
                                     rho*sin(phi(k));...
                                                  0;...
                                                  1];
                        verts(idx_verts, :) = vertex(1:3)';
                        idx_verts = idx_verts + 1;
                    end
                end
                %Prepare for the next iteration
                T_prev = T_prev*obj.Joints{i}.T(q_i, obj.Joints{i}.Parameters);
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