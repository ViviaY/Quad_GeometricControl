function plot_traj(plan_traj, track_trajs, Lp, generatedTraj, view_angle, fig_dir)
% plan_traj: 3*N matrix
% track_trajs: n_traj*3*N matrix 
% Lp: the radius of safe tube
% generatedTraj: used to generate the environment
% view_angle: the angle of views
% fig_dir: the path for saving figures

%% Plotting the normal area
f = figure;
ax1 = axes('position', [0.15 0.15 0.7 0.7]);
color = {'k','b',[0 0.5 0],'r',[0.8 0.9 0.9741],[0.8 0.98 0.9],[1 0.8 0.8],[0.7, 0.7 1]};

% plotting environment map
% Target set Xt 
Xtl = generatedTraj.Xtl; 
Xtu = generatedTraj.Xtu;

% Unsafe set Xu
XulArray = generatedTraj.XulArray; 
XuuArray = generatedTraj.XuuArray; 
R_Array = generatedTraj.R_Array;
X_Array = generatedTraj.X_Array;
Xsl = generatedTraj.Xsl;
Xsu = generatedTraj.Xsu;
xlim([Xsl(1),Xsu(1)])
ylim([Xsl(2),Xsu(2)])
zlim([Xsl(3),Xsu(3)])

sXu=size(XulArray);
Nu=sXu(2);
% plotting the target region
plotcube(transpose(Xtu-Xtl),transpose(Xtl),0.5, [0,0,1]);
hold on;

% plotting the unsafe regions
for i=1:Nu
    plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.4,[1 0 0]);
end

% number of waypoints
M = size(X_Array,2); 

% % plotting the safety regions
for k=1:M
 plotcube(transpose(2*R_Array(:,k)),transpose(X_Array(:,k)-R_Array(:,k)),0.2,[0,1,1])
end

% plotting the desried trajectory
size_traj = size(track_trajs);
n_traj = size_traj(1);
cmap = jet(n_traj);  
N = size_traj(3);
view(view_angle)
plot3(plan_traj(1,:), plan_traj(2,:), plan_traj(3,:), 'b', 'LineWidth', 2.5); % ,'LineStyle','--');
hold on;

% plotting the tube centered at the plan trajectory
for i=1:N
   center = [plan_traj(1,i),plan_traj(2,i),plan_traj(3,i)];
   if i == 1 || norm(center-center_prev) > 0.03 %0.00000003
       [X,Y,Z] = ellipsoid(center(1),center(2),center(3),Lp,Lp,Lp);
       surf(X,Y,Z,'FaceColor',color{8},'FaceAlpha',0.03,'EdgeColor','none'); %'FaceLighting','flat'
   else
       aa = 1;
   end
   center_prev = center;
end

for i = 1:n_traj
    plot3( squeeze(track_trajs(i,1,:)), ...
           squeeze(track_trajs(i,2,:)), ...
           squeeze(track_trajs(i,3,:)), ...
           'Color', cmap(i,:), ...
           'LineWidth', 1.0 );
end
hold on;


% Set labels and title
set(gca,'fontsize',15);
set(gca,'ticklabelinterpreter','latex');
xlabel('$x$ [m]','interpreter','latex');
ylabel('$y$ [m]','interpreter','latex');
zlabel('$z$ [m]','interpreter','latex');
box on;
grid on;

view(view_angle)


%% Plotting the zoom-in area (upper-right, start)
ax2 = axes('position', [0.72 0.75 0.22 0.22]);
hold on;

crop_index = floor(size_traj(3)/4);
plot3(ax2, plan_traj(1,1:crop_index), ...
    plan_traj(2,1:crop_index), ...
    plan_traj(3,1:crop_index), ...
    'black', 'LineWidth', 2,'LineStyle','-');
hold on;

for i = 1:n_traj
    plot3(ax2, squeeze(track_trajs(i,1,1:crop_index)), ...
            squeeze(track_trajs(i,2,1:crop_index)), ...
            squeeze(track_trajs(i,3,1:crop_index)), ...
           'Color', cmap(i,:), ...
           'LineWidth', 1.0,'LineStyle','--' );
end

set(gca,'ticklabelinterpreter','latex');
xlabel('$x$ [m]','interpreter','latex');
ylabel('$y$ [m]','interpreter','latex');
zlabel('$z$ [m]','interpreter','latex');
xlim(ax2, [-0.1 0.6]);  % Set x-axis limits on ax3
ylim(ax2, [-0.1 0.6]);  % Set y-axis limits on ax3
zlim(ax2, [-0.5 0.5]);  % Set z-axis limits on ax3
box on;
grid on;
hold off;
view_angle = [-138 60];
view(view_angle)

%% Plotting the zoom-in area (upper-left, target)
ax3 = axes('position', [0.08 0.75 0.22 0.22]);
hold on;
plotcube(transpose(Xtu-Xtl),transpose(Xtl), 0.2, [0,0,1]);
hold on;

% Assume track_trajs is a n_traj x 3 x N array storing multiple trajectories
n_traj = size(track_trajs, 1);  % Number of trajectories
size_traj = size(track_trajs, 3);  % Total number of points in each trajectory
end_index = size_traj;  % Total length of the trajectory
crop_index = floor(size_traj /30);  % Calculate to plot the last quarter of the trajectory
start_index = max(1, end_index - crop_index + 1);  % Compute the start index to avoid out-of-bounds error

    
% Plot the planned trajectory (assuming plan_traj is a 3xN matrix)
plot3(ax3, plan_traj(1, start_index:end_index), ...
      plan_traj(2, start_index:end_index), ...
      plan_traj(3, start_index:end_index), ...
      'black', 'LineWidth', 2, 'LineStyle', '-');

% Plot all tracking trajectories
for i = 1:n_traj
    plot3(ax3, squeeze(track_trajs(i, 1, start_index:end_index)), ...
          squeeze(track_trajs(i, 2, start_index:end_index)), ...
          squeeze(track_trajs(i, 3, start_index:end_index)), 'Color', cmap(i,:), ...
           'LineWidth', 1.0,'LineStyle','--');

    
    set(gca,'ticklabelinterpreter','latex');
    xlabel('$x$ [m]','interpreter','latex');
    ylabel('$y$ [m]','interpreter','latex');
    zlabel('$z$ [m]','interpreter','latex');
    
    xlim(ax3, [4.1 4.32]);  % Set x-axis limits on ax3
    ylim(ax3, [4.1 4.32]);  % Set y-axis limits on ax3
    zlim(ax3, [4.1 4.32]);  % Set z-axis limits on ax3
    
    box on;
    grid on;
    hold off;
    view_angle = [-123 55];
    view(view_angle)
end

% % Save PNG
% exportgraphics(gcf, fullfile(fig_dir, 'Fig5_trackingTrajs.png'), 'Resolution', 600);
% % Save EPS
% print(f, fullfile(fig_dir, 'Fig5_trackingTrajs.eps'), '-depsc');
end