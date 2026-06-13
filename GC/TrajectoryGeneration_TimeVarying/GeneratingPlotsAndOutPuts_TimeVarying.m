function GeneratingPlotsAndOutPuts_TimeVarying(param, anneal_options, generatedTraj, Bounds)

g = param.g;
mass = param.m;
vm = anneal_options.vm;
am  = anneal_options.am; 

L_p = Bounds.L_p;
L_v = Bounds.L_v;
L_f = Bounds.L_f;

epsilon = generatedTraj.epsilon;
X_Array = generatedTraj.X_Array;
R_Array = generatedTraj.R_Array;
Points_Array = generatedTraj.Points_Array;

t_Array = generatedTraj.t_Array;
XulArray = generatedTraj.XulArray;
XuuArray = generatedTraj.XuuArray;
Xtl = generatedTraj.Xtl;
Xtu = generatedTraj.Xtu;
Xsl = generatedTraj.Xsl;
Xsu = generatedTraj.Xsu;


sXu = size(XulArray);
Nu = sXu(2); % number of obstacles

M = size(X_Array,2); % number of waypoints
N_seg = M-1;% number of segments 

N = 500; % number of points per array for plotting
T = t_Array(end); % total time horizon
t = linspace(0,T,N); % time array
X = zeros(3,N); % position array
V = zeros(3,N); % velocity array
A = zeros(3,N); % acceleration array
Je = zeros(3,N); % jerk array % I'm using a different notation as J is inertia matrix
S = zeros(3,N); % snap array


tau=zeros(1,length(t_Array)-1);
for i=1:length(tau)    
   tau(i)=t_Array(i+1)-t_Array(i); 
end

% function DesiredTrajectory outputs p, pdot, pddot,pdddot, and pddddot as 
% functions of time. Function DesiredTrajectory requires the arrays 
%Points_Array (control points), and tau (duration of each segment) 

for i=1:N
    [X(:,i),V(:,i),A(:,i),Je(:,i),S(:,i)] = DesiredTrajectoryAndDerivatives(t(i),Points_Array,tau); 
    
end

% 
% figure % plotting the map only
% hold on
% plotcube(transpose(Xtu-Xtl),transpose(Xtl),0.5,[0,0,1])
% for i=1:Nu
% plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.5,[1 0 0]) 
% end
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
% view(225,80)


% 
% figure %plotting the safety tube within the map
% hold on
% 
% plotcube(transpose(Xtu-Xtl),transpose(Xtl),0.3,[0,0,1])
% for i=1:Nu
% plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.3,[1 0 0]) 
% end
% 
% for k=1:M
%  plotcube(transpose(2*R_Array(:,k)),transpose(X_Array(:,k)-R_Array(:,k)),0.5,[0,1,1])
% end
% set(gca,'fontsize',15)
% set(gca,'ticklabelinterpreter','latex')
% xlabel('$x$ [m]','interpreter','latex')
% ylabel('$y$ [m]','interpreter','latex')
% zlabel('$z$ [m]','interpreter','latex')
% xlim([Xsl(1),Xsu(1)])
% ylim([Xsl(2),Xsu(2)])
% zlim([Xsl(3),Xsu(3)])
% grid on
% view(225,80)
% %set(gcf,'renderer','Painters')
% %saveas(gcf,'Tube','epsc')
% box on


figure % plotting the desired trajectoy within the safet tube
hold on

plotcube(transpose(Xtu-Xtl),transpose(Xtl),0.3,[0,0,1])
hold on

for i=1:Nu
plotcube(transpose(XuuArray(:,i)-XulArray(:,i)),transpose(XulArray(:,i)),0.3,[1 0 0]) 
end
hold on

for k=1:M
 plotcube(transpose(2*R_Array(:,k)),transpose(X_Array(:,k)-R_Array(:,k)),0.3,[0,1,1])
end
plot3(X(1,:),X(2,:),X(3,:),'b','linewidth',2)
hold on


set(gca,'fontsize',15)
set(gca,'ticklabelinterpreter','latex')
xlabel('$x$ [m]','interpreter','latex')
ylabel('$y$ [m]','interpreter','latex')
zlabel('$z$ [m]','interpreter','latex')
xlim([Xsl(1),Xsu(1)])
ylim([Xsl(2),Xsu(2)])
zlim([Xsl(3),Xsu(3)])
box on
grid on
view(225,80)

% % Save PNG
% exportgraphics(gcf, fullfile(fig_dir, 'Fig7_v.png'), 'Resolution', 600);
% % Save EPS
% print(gcf, fullfile(fig_dir, 'Fig7_v.eps'), '-depsc');


figure
plot(t,X(1,:),'displayname','x_1','linewidth',2)
hold on
plot(t,X(2,:),'displayname','x_2','linewidth',2)
plot(t,X(3,:),'displayname','x_3','linewidth',2)
xlim([0,T])
legend

figure
plot(t,abs(V(1,:)),'displayname','v_1','linewidth',2)
hold on
plot(t,abs(V(2,:)),'displayname','v_2','linewidth',2)
plot(t,abs(V(3,:)),'displayname','v_3','linewidth',2)

vbound1=t;
vbound2=t;
vbound3=t;
for i=1:length(t)
vbound1(i)=vm(1)-L_v(t(i));
vbound2(i)=vm(2)-L_v(t(i));
vbound3(i)=vm(3)-L_v(t(i));
end
plot(t,vbound1,'displayname','v_{max}-Lv','linewidth',2)
%plot(t,vbound2,'displayname','v_{max,2}-Lv','linewidth',2)
%plot(t,vbound3,'displayname','v_{max,3}-Lv','linewidth',2)
legend

xlim([0,T])
xlabel('t')

figure



abound1=t;
abound2=t;
abound3=t;
abound4=t;
abound5=t;
for i=1:length(t)
abound1(i)=am(1);
abound2(i)=am(2);
abound3(i)=am(3)-g;
abound4(i)=-am(3)-g;
abound5(i)=(1/mass)*(L_f(t(i))+epsilon)-g;
end
subplot(3,1,1)
plot(t,A(1,:),'displayname','a_1','linewidth',2)
hold on
plot(t,abound1,':k','displayname','a_{max,1}','linewidth',2)

plot(t,-abound1,'displayname','-a_{max,1}','linewidth',2)
legend
xlim([0,T])
subplot(3,1,2)
plot(t,A(2,:),'displayname','a_2','linewidth',2)
hold on
plot(t,abound2,'displayname','a_{max,2}','linewidth',2)

plot(t,-abound2,'displayname','-a_{max,2}','linewidth',2)
legend
xlim([0,T])
subplot(3,1,3)
plot(t,A(3,:),'displayname','a_3','linewidth',2)
hold on
plot(t,abound3,'displayname','a_{max,3}-g','linewidth',2)

plot(t,abound4,'displayname','-a_{max,3}-g','linewidth',2)
plot(t,abound5,'displayname','(L_f+epsilon)/m-g','linewidth',2)

legend
xlim([0,T])
xlabel('t')

figure
plot(t,Je(1,:),'displayname','J_1','linewidth',2)
hold on
plot(t,Je(2,:),'displayname','J_2','linewidth',2)
plot(t,Je(3,:),'displayname','J_3','linewidth',2)
legend
xlim([0,T])
xlabel('t')
figure
plot(t,S(1,:),'displayname','S_1','linewidth',2)
hold on
plot(t,S(2,:),'displayname','S_2','linewidth',2)
plot(t,S(3,:),'displayname','S_3','linewidth',2)
legend
xlim([0,T])

end
