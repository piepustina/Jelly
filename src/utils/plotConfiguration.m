function [f,pl, tip_pos] = plotConfiguration(q, L, plotFrames, f, args)
%PLOTCONFIGURATION Function that plots the configuration q of the robot.
%The lenght of each Constant Curvature segment of the Piecewise Constant Curvature model
%is encoded in the vector L.

%Define the arguments of the function
arguments
    q
    L
    plotFrames (1,1) = false
    f (1,1) = figure 
    args.N (1,1) double = 100%number of points to represent each CC segment
    args.clear (1,1) = false%set to true to clean the figure
    args.lineWidth (1,1) = 3
    args.Color = 'b'
    args.A0 = eye(4,4)%define the pose of the base w.r.t. the world frame (planar)
    args.tipColor = 'b'
end

%figure(f)
set(0, 'CurrentFigure', f)

if args.clear == true
    %%hold off;
    clf(f)
end

%Get the number of CC segments of the robot
n = length(L);

%Start from the base frame (rotated since angles are counted positive clockwise)
T_prev = [args.A0(1:3,1:2) args.A0(2:4,4)];%T0 is generally for 3D so we have to extract the planar information

P = zeros(3, args.N*n);
O = zeros(2, n);
U = zeros(2, n);
V = zeros(2, n);


for i = 1 : n
    %Check if the curvature is infinite
    if abs(q(i)) < 0.0017%set a threshold for the infinite curvature (below 0.1 deg)
        T = T_prev*[cos(q(i)) -sin(q(i)) L(i);
                    sin(q(i))  cos(q(i)) 0;
                    0          0         1];
    else
        %Compute the transformation matrix from the current frame to the base
        T = T_prev*[cos(q(i)) -sin(q(i)) L(i)*sin(q(i))/q(i);
                    sin(q(i))  cos(q(i)) L(i)*(1 - cos(q(i)))/q(i);
                    0          0         1];
    end
    if plotFrames == true
        %Compute the position of origin of frame i w.r.t frame i-1
        O(:,i) = T(1:2,3);
        U(:,i) = T(1:2,1);
        V(:,i) = T(1:2,2);
    end
    %Compute the representation of the segment w.r.t. frame i-1
    [X,Y] = buildSegment(q(i), L(i), args.N);
    %Rotate the points in the world frame. Since the points computed by
    %build segment are expressed in frame i-1 we have to use T_prev.
    P(:,args.N*(i-1)+1:args.N*(i-1)+args.N) = T_prev*[X;Y;ones(1,args.N)];
    %Prepare for the next iteration
    T_prev = T;
end


%Make the plots at the end to be faster
pl = plot3(P(1,:), P(2,:), zeros(args.N*n, 1),'lineWidth', args.lineWidth, 'Color', args.Color); hold on;%plot the robot

%Return the tip position
tip_pos = [P(1, end); P(2, end)];

%Set the limits for the plot
L_max = max(L);
xlim([-n*L_max n*L_max])
ylim([-n*L_max n*L_max])

xlabel('$x$[m]', 'Interpreter', 'latex', 'FontSize', 14)
ylabel('$y$[m]', 'Interpreter', 'latex', 'FontSize', 14)
zlabel('$z$[m]', 'Interpreter', 'latex', 'FontSize', 14)

%Rotate the world to display correctly. Equivalent to set T_prev as
% T_prev = [0 1 0;
%           1 0 0;
%           0 0 1];
%Setting T_prev so requires to invert the (x-y) labels

view([0 90]);

if plotFrames == true
    %Plot the frames of each segment
    quiver(O(1,:),O(2,:) ,U(1,:),U(2,:), 'Color','r', 'lineWidth', args.lineWidth), hold on;
    quiver(O(1,:),O(2,:) ,V(1,:),V(2,:), 'Color','g', 'lineWidth', args.lineWidth), hold on;
end
grid on
%Set other properties
% xticks(-n*L_max:fix(0.5/L_max)*L_max:n*L_max)
% yticks(-n*L_max:fix(0.5/L_max)*L_max:n*L_max)
pbaspect([1 1 1])

drawnow%to show the plot also in live scripts
set(f,'Visible','on')%to show the plot as an external image

end

