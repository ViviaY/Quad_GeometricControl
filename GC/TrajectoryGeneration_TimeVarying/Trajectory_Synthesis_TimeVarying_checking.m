function generatedTraj = Trajectory_Synthesis_TimeVarying_checking(param, anneal_options, Bounds, userparam) 
rng('default');
generatedTraj = [];
mass = param.m;
vm = anneal_options.vm;
am = anneal_options.am;

L_p = Bounds.L_p;
L_v = Bounds.L_v;
L_f = Bounds.L_f;

% system dimension 
n=3;
%scaling factor for safe hyper-rectangles
if userparam.name == "alpha"
    alpha = userparam.value;
   else
    alpha = 0.99;  % Default value 
end

%scaling factor for alpha_v
if userparam.name == "alpha_s"
    alpha_s = userparam.value;
   else
    alpha_s = 0.9;  % Default value 
end


%scaling factor for rrt sampling
if userparam.name == "C_sample"
    C_sample = userparam.value;
   else
    C_sample = 0.9;  % Default value 
end


if userparam.name == "Np"
    N_pts = userparam.value;
   else
    N_pts = 10;  % Default value 
end

% maximum number of vertices in RRT plans
Nv=1000;

% speed parameter
alpha_v=2;

% paramters for ALG.2
epsilon=1e-6;

%% robustness margins
% Robustness margins adjustment to fit within a 3x3x3 box while maintaining the shape

% Scaling factor based on the original max bounds and new max bounds
scale_factor = 2 / 5;

% Safe set Xs
Xsl = [-1; -1; -1] * scale_factor;
Xsu = [5; 5; 5] * scale_factor;

Cs = 0.5 * (Xsl + Xsu);
Ds = 0.5 * (Xsu - Xsl);
Gs = diag(Ds);

% Initial point
X0 = [0; 0; 0];  % Remains the same as it is the center and within the new bounds
sX = size(X0);

%% robustness margins
% % safe set Xs
Xsl=[-1;-1;-1];
Xsu=[5;5;5];

Cs=0.5*(Xsl+Xsu);
Ds=0.5*(Xsu-Xsl);
Gs=diag(Ds);


% Initial point
X0=[0;0;0];
sX=size(X0);

% system dimension
% Target set Xt 
Xtl=[4.2;4.2;4.2];
Xtu=[4.3;4.3;4.3];
Ct=0.5*(Xtl+Xtu);
Dt=0.5*(Xtu-Xtl);
Gt=diag(Dt);


% Unsafe set Xu
XulArray=[4;5;0];
XuuArray=[4;5;3];

XulArray=[XulArray, [4.5;3.5;0]];
XuuArray=[XuuArray,[5;4;4.5]];

XulArray=[XulArray, [1.5;0;3]];
XuuArray=[XuuArray,[2.5;1;5]];

%
% XulArray=[XulArray, [0;4;0]];
% XuuArray=[XuuArray,[2;5;3]];
%

XulArray=[XulArray, [2;3;2]];
XuuArray=[XuuArray,[4;4;5]];

XulArray=[XulArray, [3;4;0]];
XuuArray=[XuuArray,[3.5;5;5]];

XulArray=[XulArray, [3.5;0;0.5]];
XuuArray=[XuuArray,[4;2.5;3]];

XulArray=[XulArray, [2.5;1;1.5]];
XuuArray=[XuuArray,[3;1.5;5]];

XulArray=[XulArray, [2.5;2.5;2]];
XuuArray=[XuuArray,[3;3;5]];

XulArray=[XulArray, [1;1;0]];
XuuArray=[XuuArray,[1.5;1.5;5]];

XulArray=[XulArray, [1;2;0]];
XuuArray=[XuuArray,[2;3;2]];
sXu=size(XulArray);
Nu=sXu(2);

CuArray=0.5*(XulArray+XuuArray);
DuArray=0.5*(XuuArray-XulArray);




%% initializing the RRT Tree and associated tree structures
tic
Tree=zeros(n,Nv);
Time_Tree=zeros(1,Nv);
Margin_Tree=zeros(n,Nv);
Safety_Radius_Tree=zeros(n,Nv);
Nodes=ones(1,Nv);
Images=ones(1,Nv);

ind=1;
ind_iter=0;
Tree(:,1)=X0; % adding initial point to the tree.
Time_Tree(1)=0;
Margin_Tree(:,1)=[L_p(Time_Tree(1));L_p(Time_Tree(1));L_p(Time_Tree(1))];
Safety_Radius_Tree(:,1)=Safety_Radius_M(Tree(:,1),Cs,Ds-Margin_Tree(:,1),CuArray,DuArray+Margin_Tree(:,1),alpha);


dist_T=300;

while  ind<=Nv 
% ind 
% generating a random sample 
   sample_safety=0;
while sample_safety==0
 % generating a random sample 
 if ind<=C_sample*Nv
    x_sample=(Xsl)+((Xsu)-(Xsl)).*rand(n,1);
 else
    x_sample=(Xtl)+((Xtu)-(Xtl)).*rand(n,1);   
 end

sample_safety = Safety_Check(x_sample,Cs,Ds,CuArray,DuArray);
end
 
 dist=100*norm(Xsu-Xsl,inf);
  for ind_dist=1:ind
      dr=norm(x_sample-ClosestPoint(x_sample,Tree(:,ind_dist),Safety_Radius_Tree(:,ind_dist)),inf); 
     if dr<dist
        ind_near=ind_dist;
        dist=dr;
     end
   
 end
 
 x_near=Tree(:,ind_near);
 R_x_near=Safety_Radius_Tree(:,ind_near);
 x_new_c=ClosestPoint(x_sample,x_near,R_x_near);
 
 ind=ind+1;
 Tree(:,ind)=x_new_c;
 Time_Tree(ind)=Time_Tree(ind_near)+norm(x_new_c-x_near)/alpha_v;
 Margin_Tree(:,ind)=[L_p(Time_Tree(ind));L_p(Time_Tree(ind));L_p(Time_Tree(ind))];
 Safety_Radius_Tree(:,ind)=Safety_Radius_M(Tree(:,ind),Cs,Ds-Margin_Tree(:,ind),CuArray,DuArray+Margin_Tree(:,ind),alpha);
 Images(ind-1)=ind;
 Nodes(ind-1)=ind_near;

 dist_T = min(Dt-Margin_Tree(:,ind)-abs(Tree(:,ind)-Ct));

if dist_T > 0
    % disp('dist_T > 0, reaching the target!');
    % fprintf('Tree(:,ind) = [%.4f, %.4f, %.4f]\n', Tree(1,ind), Tree(2,ind), Tree(3,ind));
    break;
end
 
end 


if dist_T>0
    % disp('RRT solved!')   
    G = digraph(Nodes,Images);
    Path=shortestpath(G,1,ind);
    M=length(Path);
    X_Array=zeros(n,M);
    t_c1=toc;
else
  return;
end
  
 
for k=1:M
    X_Array(:,k)=Tree(:,Path(k));
end

%% Generating Bezier Curves
tic; 
[Points_Array,R_Array,t_Array, tau, param_alphas] = BezierControlPoints_IterativeLP_TimeVarying_orig(epsilon,alpha,alpha_v,alpha_s,X_Array,Cs,Ds,Ct,Dt,CuArray,DuArray,N_pts,vm,am,L_p, L_v, L_f, mass);             
t_c2=toc;

generatedTraj.C_sample = C_sample;
generatedTraj.Ns = size(X_Array,2);
generatedTraj.Np = N_pts;
generatedTraj.Nv = Nv;
generatedTraj.epsilon = epsilon;
generatedTraj.alpha = param_alphas.alpha;
generatedTraj.alpha_s = param_alphas.alpha_s;
generatedTraj.alpha_v = param_alphas.alpha_v;

generatedTraj.X_Array = X_Array; 
generatedTraj.R_Array = R_Array;
generatedTraj.Points_Array = Points_Array;

generatedTraj.t_Array = t_Array;
generatedTraj.tau = tau;
generatedTraj.XulArray = XulArray; 
generatedTraj.XuuArray = XuuArray; 
generatedTraj.Xtl = Xtl; 
generatedTraj.Xtu = Xtu; 
generatedTraj.Xsl = Xsl; 
generatedTraj.Xsu = Xsu; 
generatedTraj.RRT_time = t_c1;
generatedTraj.Bezier_ILPtime = t_c2;


fprintf("cpu time: obtainning safe tube with RRT: %.3f, " + ...
        "obtain Bezier curve from iterative linear programming: %.3f \n", t_c1, t_c2);

end
 
