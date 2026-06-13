function generate_outputs_plots(generated_trajectory, view_angle)

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


sXu=size(XulArray);
Nu=sXu(2); % number of obstacles

M=size(X_Array,2); % number of waypoints
N_seg=M-1;% number of segments 

N=100; % number of points per array for plotting


T=sum(tau); % total time horizon
t=linspace(0,T,N); % time array
X=zeros(3,N); % position array
V=zeros(3,N); % velocity array
A=zeros(3,N); % acceleration array
Je=zeros(3,N); % jerk array % I'm using a different notation as J is inertia matrix
S=zeros(3,N); % snap array

for i=1:N
    [X(:,i),V(:,i),A(:,i),Je(:,i),S(:,i)]=DesiredTrajectoryAndDerivatives(t(i),Points_Array,tau); 
    
end

% %% plotting the map only
% figure 
% ax1 = axes('position', [0.2 0.15 0.65 0.65]);
% hold on
% plotcube(transpose(Xtu-Xtl),transpose(Xtl),0.5,[0,0,1])
% for i=1:Nu
% % plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.5,[1 0 0]) 
%     plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.4,[1 0 0]);
% end
% 
% % scatter3(0.5,0.5,1,'*')
% scatter3(0,0,0,'*', 'b')
% 
% 
% set(gca,'fontsize',15)
% set(gca,'ticklabelinterpreter','latex')
% xlabel('$x$ [m]','interpreter','latex')
% ylabel('$y$ [m]','interpreter','latex')
% zlabel('$z$ [m]','interpreter','latex')
% xlim([Xsl(1),Xsu(1)])
% ylim([Xsl(2),Xsu(2)])
% zlim([Xsl(3),Xsu(3)])
% box on
% grid on
% view(view_angle)
% 
% set(gcf,'renderer','Painters')
% 
% % subplot 1
% ax2 = axes('position', [0.7 0.7 0.2 0.2]);
% hold on
% plotcube(transpose(Xtu-Xtl),transpose(Xtl),0.5,[0,0,1])
% for i=1:Nu
% plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.5,[1 0 0]) 
% end
% % scatter3(0.5,0.5,1,'*')
% scatter3(0,0,0,'*', 'b')
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
% ax3 = axes('position', [0.15 0.7 0.2 0.2]);
% hold on
% plotcube(transpose(Xtu-Xtl),transpose(Xtl),0.5,[0,0,1])
% for i=1:Nu
% % plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.5,[1 0 0]) 
% plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.4,[1 0 0]);
% end
% % scatter3(0.5,0.5,1,'*')
% scatter3(0,0,0,'*', 'b')
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
% 
% % set(gcf,'renderer','Painters')
% %Save as a vector image in PDF format
% % print('Figs/Fig1_Env', '-dpdf', '-r600'); 
% % exportgraphics(gcf, 'Figs/Fig1_Env.pdf', 'ContentType', 'vector');
% % exportgraphics(gcf, 'Figs/Fig1_Env.pdf', 'ContentType', 'image', 'Resolution', 600);
% exportgraphics(gcf, 'Figs/Fig1_Env.png', 'Resolution', 600);
% exportgraphics(gcf, 'Figs/Fig1_Env.eps', 'ContentType', 'vector');
%%  plotting the safety tube within the map
f = figure;
hold on

plotcube(transpose(Xtu-Xtl),transpose(Xtl),0.05,[0,0,1])

% unsafe 
for i=1:Nu
plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.4,[1 0 0]) 
end

% safe 
for k=1:M
 plotcube(transpose(2*R_Array(:,k)),transpose(X_Array(:,k)-R_Array(:,k)),0.2,[0,1,1])
end
% plot3(X(1,:),X(2,:),X(3,:),'b','linewidth',2)
% scatter3(0, 0, 0, 30, 'm', '*', 'LineWidth', 0.5);
scatter3(0, 0, 0, 80, 'Marker', '*', 'MarkerEdgeColor', '#32CD32', 'LineWidth', 1); 


hold on

plot3(X(1,:),X(2,:),X(3,:),'b','linewidth',2)
hold on
set(gca, 'Position', [0.21, 0.22, 0.65, 0.65]);
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

% subplot 1
% ax11 = axes('position', [0.78 0.78 0.2 0.2]);
% ax11 = axes('position', [0.78 0.01 0.2 0.2]);
ax11 = axes('position', [0.1 0.75 0.2 0.2]);  
hold on
plotcube(transpose(Xtu-Xtl),transpose(Xtl),0.5,[0,0,1])
for i=1:Nu
plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.5,[1 0 0]) 
end
% scatter3(0.5,0.5,1,'*')
% scatter3(0,0,0,'*', 'b')
% scatter3(0, 0, 0, 30, 'm', '*', 'LineWidth', 0.5);
scatter3(0, 0, 0, 80, 'Marker', '*', 'MarkerEdgeColor', '#228B22', 'LineWidth', 1); 

% set(gca,'fontsize',15)
set(gca,'ticklabelinterpreter','latex')
xlabel('$x$ [m]','interpreter','latex')
ylabel('$y$ [m]','interpreter','latex')
zlabel('$z$ [m]','interpreter','latex')
% xlim([Xsl(1),Xsu(1)])
% ylim([Xsl(2),Xsu(2)])
xlim([-1,Xsu(1)])
ylim([-1,Xsu(2)])
zlim([Xsl(3),Xsu(3)])
box on
grid on
view([270 0])

set(gcf,'renderer','Painters')

% subplot 2
% ax12 = axes('position', [0.08 0.78 0.2 0.2]);
% ax12 = axes('position', [0.0 0.0 0.2 0.2]);
ax12 = axes('position', [0.75 0.75 0.2 0.2]);  
hold on
plotcube(transpose(Xtu-Xtl),transpose(Xtl),0.5,[0,0,1])
for i=1:Nu
% plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.5,[1 0 0]) 
plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.4,[1 0 0]);
end

% scatter3(0, 0, 0, 30, 'm', '*', 'LineWidth', 0.5);
scatter3(0, 0, 0, 80, 'Marker', '*', 'MarkerEdgeColor', '#228B22', 'LineWidth', 1); 

% set(gca,'fontsize',15)
set(gca,'ticklabelinterpreter','latex')
xlabel('$x$ [m]','interpreter','latex')
ylabel('$y$ [m]','interpreter','latex')
zlabel('$z$ [m]','interpreter','latex')
% xlim([Xsl(1),Xsu(1)])
% ylim([Xsl(2),Xsu(2)])
xlim([-1,Xsu(1)])
ylim([-1,Xsu(2)])
zlim([Xsl(3),Xsu(3)])

box on
grid on
view([180 0])



if zoom_in

    % ax21 = axes('position', [0.78 0.78 0.2 0.2]);
    % ax21 = axes('position', [0.78 0.08 0.2 0.2]);
    ax21 = axes('position', [0.75 0.1 0.2 0.2]);
    hold on 
    % plotting the safety regions
    for k=1:1
     plotcube(transpose(2*R_Array(:,k)),transpose(X_Array(:,k)-R_Array(:,k)),0.8,[0,1,1])
    end
    hold on
    size_traj = size(X);
    crop_index = floor(size_traj(2)/4);
    plot3(ax21, X(1,1:crop_index+2), ...
        X(2,1:crop_index+2), ...
        X(3,1:crop_index+2), ...
        'b', 'LineWidth', 2,'LineStyle','-');
    hold on;

    % scatter3(0, 0, 0, 30, 'm', '*', 'LineWidth', 0.5);
    scatter3(0, 0, 0, 80, 'Marker', '*', 'MarkerEdgeColor', '#228B22', 'LineWidth', 1); 


    set(gca,'ticklabelinterpreter','latex');
    xlabel('$x$ [m]','interpreter','latex');
    ylabel('$y$ [m]','interpreter','latex');
    zlabel('$z$ [m]','interpreter','latex');
    xlim(ax21, [-0.2 0.2]);  % Set x-axis limits on ax3
    ylim(ax21, [-0.2 0.2]);  % Set y-axis limits on ax3
    zlim(ax21, [-0.2 0.2]);  % Set z-axis limits on ax3
    box on;
    grid on;
    hold off;
    % view_angle = [-130 24];
    view_angle = [-130 41];
    view(view_angle)
end


if zoom_in
    % ax22 = axes('position', [0.0 0.78 0.2 0.2]);
    % ax22 = axes('position', [0.08 0.08 0.2 0.2]);
    ax22 = axes('position', [0.1 0.1 0.2 0.2]); 
    hold on
    plotcube(transpose(Xtu-Xtl),transpose(Xtl),0.2, [0,0,1]);
    hold on

    % k = size(R_Array, 2);
    for k = size(R_Array, 2):size(R_Array, 2) %size(R_Array, 2)-1: size(R_Array, 2) %
    plotcube(transpose(2*R_Array(:,k)),transpose(X_Array(:,k)-R_Array(:,k)),0.2,[0,1,1])
    end
    % Assume real_traj is a n_traj x 3 x N array storing multiple trajectories
    size_traj = size(X, 2);  % Total number of points in each trajectory
    end_index = size_traj;  % Total length of the trajectory
    crop_index = floor(size_traj /30);  % Calculate to plot the last quarter of the trajectory
    start_index = max(1, end_index - crop_index + 1);  % Compute the start index to avoid out-of-bounds error
   
    
    % Plot the planned trajectory (assuming plan_traj is a 3xN matrix)
    plot3(ax22, X(1, start_index:end_index), ...
          X(2, start_index:end_index), ...
          X(3, start_index:end_index), ...
          'b', 'LineWidth', 2, 'LineStyle', '-');
    

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

    xlim(ax22, [4.22 4.32]);  % Set x-axis limits on ax3
    ylim(ax22, [4.22 4.32]);  % Set y-axis limits on ax3
    zlim(ax22, [4.22 4.32]);  % Set z-axis limits on ax3
    % zlim(ax22, [3.9 4.4]);  % Set z-axis limits on ax3

    box on;
    grid on;
    hold off;
    % view_angle = [-130 24];
    view_angle = [-130 41];
    view(view_angle)
end

% exportgraphics(f, 'Figs/Fig2_SafeTube2.png', 'Resolution', 600);
print(f, 'Figs/Fig2_SafeTube2.eps', '-depsc');
% exportgraphics(f, 'Figs/Fig2_SafeTube2.eps', 'ContentType', 'vector');
%% plotting the desired trajectoy within the safet tube
% figure 
% hold on
% for k=1:M
%  plotcube(transpose(2*R_Array(:,k)),transpose(X_Array(:,k)-R_Array(:,k)),0.1,[0,1,1])
% end
% plot3(X(1,:),X(2,:),X(3,:),'b','linewidth',2)
% hold on
% 
% 
% set(gca,'fontsize',15)
% set(gca,'ticklabelinterpreter','latex')
% xlabel('$x$ [m]','interpreter','latex')
% ylabel('$y$ [m]','interpreter','latex')
% zlabel('$z$ [m]','interpreter','latex')
% xlim([Xsl(1),Xsu(1)])
% ylim([Xsl(2),Xsu(2)])
% zlim([Xsl(3),Xsu(3)])
% box on
% grid on
% view(view_angle)
% 
% % set(gcf,'renderer','Painters')
% %Save as a vector image in PDF format
% % print('Figs/Fig3_DesiredTraj', '-dpdf', '-r600'); 
% % exportgraphics(gcf, 'Figs/Fig3_DesiredTraj.pdf', 'ContentType', 'vector');
% 
% % exportgraphics(gcf, 'Figs/Fig3_DesiredTraj.pdf', 'ContentType', 'image', 'Resolution', 600);
% 
% exportgraphics(gcf, 'Figs/Fig3_DesiredTraj.png', 'Resolution', 600);
% exportgraphics(gcf, 'Figs/Fig3_DesiredTraj.eps', 'ContentType', 'vector');
% end
% 
% 
