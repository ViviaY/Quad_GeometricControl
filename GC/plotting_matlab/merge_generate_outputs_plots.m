function merge_generate_outputs_plots(t_span, generated_trajectory, plan_traj, track_trajs, view_angle, fig_dir)
% plan_traj: 3*N matrix
% track_trajs: n_traj*3*N matrix 
% generatedTraj: generating the environment with obstacles and safe boxes
% zoom_in: number treated as Boolean. plot the zoom-in traj to show convergence
%          0 for no zoom_in image. positive integer for first 1/zoom_in
%          section of trajectories
% view_angle: the angle of views
% fig_dir: the path for saving figures

%% setting the environments 
X_Array = generated_trajectory.X_Array;
R_Array = generated_trajectory.R_Array;
Points_Array = generated_trajectory.Points_Array;
tau = generated_trajectory.tau;
XulArray = generated_trajectory.XulArray;
XuuArray = generated_trajectory.XuuArray;
Xtl = generated_trajectory.Xtl;
Xtu = generated_trajectory.Xtu;
Xsl = generated_trajectory.Xsl;
Xsu = generated_trajectory.Xsu;

colors_setting = ColorsData();
colors = colors_setting.colors;


sXu = size(XulArray);
Nu = sXu(2); % number of obstacles

M = size(X_Array,2); % number of waypoints
N_seg = M-1;% number of segments 

N = 100; % number of points per array for plotting


T = sum(tau); % total time horizon
t = linspace(0, T, N); % time array
X = zeros(3, N); % position array
V = zeros(3, N); % velocity array
A = zeros(3, N); % acceleration array
Je = zeros(3, N); % jerk array % using a different notation as J is inertia matrix
S = zeros(3, N); % snap array

for i=1:N
    [X(:,i),V(:,i),A(:,i),Je(:,i),S(:,i)] = DesiredTrajectoryAndDerivatives(t_span(i),Points_Array,tau); 
end

%%  plotting the safety tube within the map
% f = figure;
f = figure('Units', 'inches', 'Position', [1, 1, 8, 6]);
hold on

plotcube(transpose(Xtu-Xtl),transpose(Xtl),0.05,[0,0,1])

% unsafe 
for i=1:Nu
plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.4,[1 0 0]) 
end
hold on;

size_traj = size(track_trajs);
n_traj = size_traj(1);
for i=1:n_traj
    if n_traj > 20
        colors =  jet(n_traj+1);
    end
    plot3(squeeze(track_trajs(i,1,:)), squeeze(track_trajs(i,2,:)), squeeze(track_trajs(i,3,:)), 'Color', colors(i+1, :), 'LineWidth', 1);
end
hold on;

% safe 
for k=1:M
 plotcube(transpose(2*R_Array(:,k)),transpose(X_Array(:,k)-R_Array(:,k)),0.1,[0,1,1])
end
% plot3(X(1,:),X(2,:),X(3,:),'b','linewidth',2)
% scatter3(0, 0, 0, 30, 'm', '*', 'LineWidth', 0.5);
scatter3(0, 0, 0, 80, 'Marker', '*', 'MarkerEdgeColor', '#32CD32', 'LineWidth', 1); 


hold on

% plot3(X(1,:),X(2,:),X(3,:),'b','linewidth',2)
plot3(plan_traj(1,:), plan_traj(2,:), plan_traj(3,:), 'b', 'LineWidth', 2); 
hold on
% set(gca, 'Position', [0.21, 0.22, 0.65, 0.65]);
set(gca, 'Position', [0.2, 0.1, 0.65, 0.65]);
set(gca,'fontsize',15)
set(gca,'ticklabelinterpreter','latex', fontsize=12)

hx = xlabel('$x$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
hy = ylabel('$y$ [m]', 'Interpreter', 'latex', 'FontSize', 12);
hz = zlabel('$z$ [m]', 'Interpreter', 'latex', 'FontSize', 12);

% hx.Position(2) = hx.Position(2) - 0.2; % Move x-label down
% hy.Position(1) = hy.Position(1) - 0.2; % Move y-label left
% hz.Position(2) = hz.Position(2) + 0.2; % Move z-label up

hx.HorizontalAlignment = 'center';  % or 'left' / 'right'
hy.HorizontalAlignment = 'center';  % or 'left' / 'right'
hz.HorizontalAlignment = 'center';  % or 'left' / 'right'


xlim([Xsl(1),Xsu(1)])
ylim([Xsl(2),Xsu(2)])
zlim([Xsl(3),Xsu(3)])
view(view_angle)
box on
grid on

% set(gcf,'renderer','Painters')
%Save as a vector image in PDF format
% print('Figs/Fig2_SafeTube', '-dpdf', '-r600'); 
% exportgraphics(gcf, 'Figs/Fig2_SafeTube.pdf', 'ContentType', 'vector');
% exportgraphics(gcf, 'Figs/Fig2_SafeTube.pdf', 'ContentType', 'image', 'Resolution', 600);
zoom_in  = true;
%Plotting the zoom-in area 
% 
% % subplot 1
% % ax11 = axes('position', [0.78 0.78 0.2 0.2]);
% % ax11 = axes('position', [0.78 0.01 0.2 0.2]);
% ax11 = axes('position', [0.1 0.75 0.2 0.2]);  
% hold on
% plotcube(transpose(Xtu-Xtl),transpose(Xtl),0.5,[0,0,1])
% for i=1:Nu
% plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.5,[1 0 0]) 
% end
% % scatter3(0.5,0.5,1,'*')
% % scatter3(0,0,0,'*', 'b')
% % scatter3(0, 0, 0, 30, 'm', '*', 'LineWidth', 0.5);
% scatter3(0, 0, 0, 80, 'Marker', '*', 'MarkerEdgeColor', '#228B22', 'LineWidth', 1); 
% 
% % set(gca,'fontsize',15)
% set(gca,'ticklabelinterpreter','latex')
% xlabel('$x$ [m]','interpreter','latex')
% ylabel('$y$ [m]','interpreter','latex')
% zlabel('$z$ [m]','interpreter','latex')
% % xlim([Xsl(1),Xsu(1)])
% % ylim([Xsl(2),Xsu(2)])
% xlim([-1,Xsu(1)])
% ylim([-1,Xsu(2)])
% zlim([Xsl(3),Xsu(3)])
% box on
% grid on
% view([270 0])
% 
% set(gcf,'renderer','Painters')
% 
% % subplot 2
% % ax12 = axes('position', [0.08 0.78 0.2 0.2]);
% % ax12 = axes('position', [0.0 0.0 0.2 0.2]);
% ax12 = axes('position', [0.75 0.75 0.2 0.2]);  
% hold on
% plotcube(transpose(Xtu-Xtl),transpose(Xtl),0.5,[0,0,1])
% for i=1:Nu
% % plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.5,[1 0 0]) 
% plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.4,[1 0 0]);
% end
% 
% % scatter3(0, 0, 0, 30, 'm', '*', 'LineWidth', 0.5);
% scatter3(0, 0, 0, 80, 'Marker', '*', 'MarkerEdgeColor', '#228B22', 'LineWidth', 1); 
% 
% % set(gca,'fontsize',15)
% set(gca,'ticklabelinterpreter','latex')
% xlabel('$x$ [m]','interpreter','latex')
% ylabel('$y$ [m]','interpreter','latex')
% zlabel('$z$ [m]','interpreter','latex')
% % xlim([Xsl(1),Xsu(1)])
% % ylim([Xsl(2),Xsu(2)])
% xlim([-1,Xsu(1)])
% ylim([-1,Xsu(2)])
% zlim([Xsl(3),Xsu(3)])
% 
% box on
% grid on
% view([180 0])



if zoom_in
    % ax21 = axes('position', [0.78 0.78 0.2 0.2]);
    % ax21 = axes('position', [0.78 0.08 0.2 0.2]);
    % ax21 = axes('position', [0.75 0.1 0.2 0.2]);
    ax21 = axes('position', [0.75 0.65 0.2 0.2]); 
    hold on 
    % plotting the safety regions
    % for k=1:1
    %  plotcube(transpose(2*R_Array(:,k)),transpose(X_Array(:,k)-R_Array(:,k)),0.1,[0,1,1])
    % end
    hold on
    size_traj = size(X);
    % crop_index = floor(size_traj(2)/2);
    % plot3(ax21, X(1,1:crop_index+2), ...
    %     X(2,1:crop_index+2), ...
    %     X(3,1:crop_index+2), ...
    %     'b', 'LineWidth', 2,'LineStyle','-');
    % plot3(ax21, plan_traj(1,1:crop_index+2), plan_traj(2,1:crop_index+2), plan_traj(3,1:crop_index+2), 'b', 'LineWidth', 3, 'LineStyle', '-'); 
    plot3(ax21, plan_traj(1,:), plan_traj(2,:), plan_traj(3,:), 'b', 'LineWidth', 4, 'LineStyle', '-'); 
    hold on;

    for i=1:n_traj
        % plot3(squeeze(real_traj(i,1,1:crop_index)), squeeze(real_traj(i,2,1:crop_index)), squeeze(real_traj(i,3,1:crop_index)), 'Color', colors(i, :), 'LineWidth', 1);
        if n_traj > 20
            colors =  jet(n_traj+1);
        end
        plot3(ax21,squeeze(track_trajs(i,1,:)), squeeze(track_trajs(i,2,:)), squeeze(track_trajs(i,3,:)), 'Color', colors(i+1, :), 'LineWidth', 2);
        % plot3(ax21,squeeze(real_traj(i,1,:)), squeeze(real_traj(i,2,:)), squeeze(real_traj(i,3,:)), 'Color', 'b', 'LineWidth', 1);
    end
    hold on;
    % scatter3(0, 0, 0, 30, 'm', '*', 'LineWidth', 0.5);
    % scatter3(0, 0, 0, 80, 'Marker', '*', 'MarkerEdgeColor', '#228B22', 'LineWidth', 1); 


    set(gca,'ticklabelinterpreter','latex');
    xlabel('$x$ [m]','interpreter','latex');
    ylabel('$y$ [m]','interpreter','latex');
    zlabel('$z$ [m]','interpreter','latex');
    % xlim(ax21, [-0.4 0.6]);  % Set x-axis limits on ax3
    % ylim(ax21, [-0.4 0.6]);  % Set y-axis limits on ax3
    % zlim(ax21, [-0.4 1.2]);  % Set z-axis limits on ax3
    xlim(ax21, [-0.1    0.6]);
    ylim(ax21, [-0.1816    0.58]);
    zlim(ax21, [-0.5843    0.9]);
    box on;
    grid on;
    hold off;
    % view_angle = [-130 24];
    % view_angle = [-130 41];
    % view_angle = [-65.8042 38.5717];
    % view_angle = [-74.3938 38.5717];
    % view_angle = [-139.9459, 32.5009];
    % view_angle = [-219.0605, 20.3594]; % good
    view_angle = [-142.2064, 49.7816]; % good

    % view_angle = [ -138.5897, 26.7126]; %34.9293
    view(view_angle)
end


if zoom_in
    % ax22 = axes('position', [0.0 0.78 0.2 0.2]);
    % ax22 = axes('position', [0.08 0.08 0.2 0.2]);
    % ax22 = axes('position', [0.1 0.1 0.2 0.2]); 
    ax22 = axes('position', [0.1 0.65 0.2 0.2]); 
    hold on

    plotcube(transpose(Xtu-Xtl),transpose(Xtl), 0.2, [0,0,1]);
    hold on

    % k = size(R_Array, 2);
    % for k = size(R_Array, 2):size(R_Array, 2) %size(R_Array, 2)-1: size(R_Array, 2) %
    % plotcube(transpose(2*R_Array(:,k)),transpose(X_Array(:,k)-R_Array(:,k)),0.1,[0,1,1])
    % end
    % Assume real_traj is a n_traj x 3 x N array storing multiple trajectories
    size_traj = size(X, 2);  % Total number of points in each trajectory
    end_index = size_traj;  % Total length of the trajectory
    crop_index = floor(size_traj /1);  % Calculate to plot the last quarter of the trajectory
    start_index = max(1, end_index - crop_index + 1);  % Compute the start index to avoid out-of-bounds error
   
    
    % Plot the planned trajectory (assuming plan_traj is a 3xN matrix)
    % plot3(ax22, X(1,:), ...
    %       X(2, :), ...
    %       X(3, :), ...
    %       'b', 'LineWidth', 3);
    
    plot3(ax22, plan_traj(1,:), plan_traj(2,:), plan_traj(3,:), 'b', 'LineWidth', 2, 'LineStyle', '-'); 
    hold on;
    for i=1:n_traj
        if n_traj > 20
            colors =  jet(n_traj+1);
        end
        plot3(ax22, squeeze(track_trajs(i,1,:)), squeeze(track_trajs(i,2,:)), squeeze(track_trajs(i,3,:)), 'Color', colors(i+1, :), 'LineWidth', 1);
    end
    hold on;

    set(gca,'ticklabelinterpreter','latex');
    xlabel('$x$ [m]','interpreter','latex');
    ylabel('$y$ [m]','interpreter','latex');
    zlabel('$z$ [m]','interpreter','latex');

    % Adjust properties
    % set(ax3.XLabel, 'Rotation', 25, 'VerticalAlignment', 'middle', 'HorizontalAlignment', 'right');
    % set(ax3.YLabel, 'Rotation', -25, 'VerticalAlignment', 'bottom', 'HorizontalAlignment', 'right');
    % 
    % % Customize the label padding
    % ax3.XAxis.Label.FontSize = 12;  % Adjust the font size
    % ax3.XAxis.Label.Position(2) = ax3.XAxis.Label.Position(2) - 0.1; % Move closer to the axis

    % xlim(ax22, [4.2211    4.3324]);  % Set x-axis limits on ax3
    % ylim(ax22, [4.1834    4.3042]);  % Set y-axis limits on ax3
    % zlim(ax22, [4.1992    4.3721]);  % Set z-axis limits on ax3


    xlim(ax22, [4.1    4.42]);  % Set x-axis limits on ax3
    ylim(ax22, [4.14    4.3113]);  % Set y-axis limits on ax3
    zlim(ax22, [3.98    4.3320]);  % Set z-axis limits on ax3
    % zlim(ax22, [3.9 4.4]);  % Set z-axis limits on ax3

    box on;
    grid on;
    hold off;
    % view_angle = [-130 24];
    % view_angle = [-130 41];
    % view_angle = [-195.0998, 17.9310];
    % view_angle = [-65.8042 38.5717];
    % view_angle = [-74.3938 38.5717];
    % view_angle = [-155.3170, 10.6461]; % good
     % view_angle = [-142.2064, 49.7816]; % good
     view_angle = [-117, 54];
    % view_angle = [ -129.5480,  17.1405]; 
    view(view_angle)
end

% Save PNG
exportgraphics(gcf, fullfile(fig_dir, 'Fig12_merge.png'), 'Resolution', 600);
% Save EPS
print(gcf, fullfile(fig_dir, 'Fig12_merge.eps'), '-depsc');


