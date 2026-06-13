function Bounds = ObtainingBoundsTimeVarying(annealing_output, anneal_options, param, plotting, save_dir)
% Setting parameters
% e3 = [0;0;1];
J = param.J;
lambda_m_J = min(eig(J));
% lambda_M_J = max(eig(J));
m = param.m;
g = param.g;
vm = anneal_options.vm;
am = anneal_options.am;

psi_bar = anneal_options.psi_bar;
alpha_psi = anneal_options.alpha_psi;
V1_0 = anneal_options.V1_0;


kp = annealing_output.opt_k(1);
kv = annealing_output.opt_k(2);
kR = annealing_output.opt_k(3);
kw = annealing_output.opt_k(4);
gamma1 = annealing_output.opt_k(5);
gamma2 = annealing_output.opt_k(6);

c1 = gamma1 * min([sqrt(kp*m), 4*m*kp*kv/(kv^2+4*m*kp)]);
c2 = gamma2 * min([sqrt(kR*lambda_m_J), 4*(lambda_m_J)*kR*kw/((kw^2)+4*(lambda_m_J)*kR)]);


tic;
% this is M1
M1t=0.5*[kp*eye(3,3), c1*eye(3,3); c1*eye(3,3), m*eye(3,3)];

% this is W1
W1t=[(1/m)*c1*kp*eye(3,3) (0.5/m)*c1*kv*eye(3,3);(0.5/m)*c1*kv*eye(3,3) (kv-c1)*eye(3,3)];

% this is M2_1
M2_1t=0.5*[kR*eye(3,3) c2*eye(3,3);c2*eye(3,3) J];

% this is M2_2
M2_2t=0.5*[(2*kR/(2-psi_bar))*eye(3,3) c2*eye(3,3);c2*eye(3,3) J];

% this is W2
W2t=[c2*kR*inv(J) 0.5*c2*kw*inv(J);0.5*c2*kw*inv(J) (kw-c2)*eye(3,3)];

% bound on V2(0) in terms of psi_1 (bar{V}_{2})
V2_0=(kR+2*c2*sqrt((kR/lambda_m_J)*alpha_psi*(1-alpha_psi)))*psi_bar; 

% bound on V(0)
V_0=V1_0+V2_0; 


% computing the coefficients of the bound
x1t=min(eig(inv(sqrtm(M1t))*W1t*inv(sqrtm(M1t))));

x2t=min(eig(inv(sqrtm(M2_2t))*W2t*inv(sqrtm(M2_2t))));

% this corresponds to beta
ct=x2t/2;

% this corresponds to alpha_0
a0t=min([x1t,x2t]);

% this corresponds to alpha_1*sqrt(V_2(0))
a1t=sqrt(2/(2-psi_bar))*norm([c1/m*eye(3) eye(3)]*inv(sqrtm(M1t)))*norm([kp*eye(3) kv*eye(3)]*inv(sqrtm(M1t)))*norm([eye(3) zeros(3,3)]*inv(sqrtm(M2_1t)))*sqrt(V2_0);

% this corresponds to alpha_2*sqrt(V_2(0))
a2t=m*norm(am)*sqrt(2/(2-psi_bar))*norm([c1/m*eye(3,3) eye(3,3)]*inv(sqrtm(M1t)))*norm([eye(3,3) zeros(3,3)]*inv(sqrtm(M2_1t)))*sqrt(V2_0);


% finding tm and the uniform bound depending on 0.5*alpha_0 and beta
Bound0t=(exp(0.5*a1t/ct));

if 0.5*a0t~=ct
    X=0.5*a2t*ct/(ct-0.5*a0t);
    Y=0.5*a0t*sqrt(V_0)+(0.5*a0t*0.5*a2t/(ct-0.5*a0t));
    tbt1=log(X/Y)/(ct-0.5*a0t);
    tbt=max([tbt1,0]);
    Bound1t=0.5*a2t*((exp(-ct*tbt)-exp(-0.5*a0t*tbt))/(0.5*a0t-ct));
    L_1=@(x)(Bound0t*(sqrt(V_0)*exp(-0.5*a0t*x)+0.5*a2t*((exp(-ct*x)-exp(-0.5*a0t*x))/(0.5*a0t-ct))));
    L_p_1=@(x)norm([eye(3,3) zeros(3,3)]*inv(sqrtm(M1t)))*(Bound0t*(sqrt(V_0)*exp(-0.5*a0t*x)+0.5*a2t*((exp(-ct*x)-exp(-0.5*a0t*x))/(0.5*a0t-ct))));
    L_v_1=@(x)norm([zeros(3,3),eye(3,3)]*inv(sqrtm(M1t)))*(Bound0t*(sqrt(V_0)*exp(-0.5*a0t*x)+0.5*a2t*((exp(-ct*x)-exp(-0.5*a0t*x))/(0.5*a0t-ct))));
    L_f_1=@(x)norm([kp*eye(3,3),kv*eye(3,3)]*inv(sqrtm(M1t)))*(Bound0t*(sqrt(V_0)*exp(-0.5*a0t*x)+0.5*a2t*((exp(-ct*x)-exp(-0.5*a0t*x))/(0.5*a0t-ct))));

else
    tbt1=2*(a2t-a0t*sqrt(V_0))/(a0t*a2t);
    tbt=max([tbt1,0]);
    Bound1t=0.5*a2t*tbt*exp(-0.5*a0t*tbt);
    L_1=@(x)(Bound0t*(sqrt(V_0)*exp(-0.5*a0t*x)+0.5*a2t*x*exp(-0.5*a0t*x)));
    L_p_1=@(x)norm([eye(3,3) zeros(3,3)]*inv(sqrtm(M1t)))*(Bound0t*(sqrt(V_0)*exp(-0.5*a0t*x)+0.5*a2t*x*exp(-0.5*a0t*x)));
    L_v_1=@(x)norm([zeros(3,3),eye(3,3)]*inv(sqrtm(M1t)))*(Bound0t*(sqrt(V_0)*exp(-0.5*a0t*x)+0.5*a2t*x*exp(-0.5*a0t*x)));
    L_f_1=@(x)norm([kp*eye(3,3),kv*eye(3,3)]*inv(sqrtm(M1t)))*(Bound0t*(sqrt(V_0)*exp(-0.5*a0t*x)+0.5*a2t*x*exp(-0.5*a0t*x)));

end

% bounds on sqrt(V) 
sqrt_V_Uniform_Bound_t=Bound0t*(sqrt(V_0)*exp(-0.5*a0t*tbt)+Bound1t);

%Lp bound
L_2 = @(x)sqrt_V_Uniform_Bound_t+0*x;
L_p_2 = @(x) norm([eye(3,3) zeros(3,3)]*inv(sqrtm(M1t)))*sqrt_V_Uniform_Bound_t+0*x;
L_v_2 = @(x) norm([zeros(3,3),eye(3,3)]*inv(sqrtm(M1t)))*sqrt_V_Uniform_Bound_t+0*x;
L_f_2 = @(x) norm([kp*eye(3,3),kv*eye(3,3)]*inv(sqrtm(M1t)))*sqrt_V_Uniform_Bound_t+0*x;


L = @(x) fg_mfile(x,L_1,L_2,tbt);
L_p = @(x) fg_mfile(x,L_p_1,L_p_2,tbt);
L_v = @(x) fg_mfile(x,L_v_1,L_v_2,tbt);
L_f = @(x) fg_mfile(x,L_f_1,L_f_2,tbt);

% Bounds.sqrt_V_Uniform_Bound_t = sqrt_V_Uniform_Bound_t;
Bounds.L = L;
Bounds.L_p = L_p;
Bounds.L_v = L_v;
Bounds.L_f = L_f;

if exist('save_dir','var') && ~isempty(save_dir)
    save(fullfile(save_dir, "TimeVaryingBoundsData.mat"), ...
                  "L", "L_p","L_v","L_f","m","J","vm","am","kp","kv","kR","kw","psi_bar","alpha_psi","V1_0","V2_0","c1","c2");
end
% plot the bounds:
if plotting
    PlotBounds(3, L_p, L_p_1, L_p_2, L_v, L_v_1, L_v_2, L_f, L_f_1, L_f_2);
else
    disp("Plotting is disabled");
end

end

%% plotting the time-varying and uniform bound for t \in [0,3] as an example
function PlotBounds(T, L_p, L_p_1, L_p_2, L_v, L_v_1, L_v_2, L_f, L_f_1, L_f_2)
hold on
t = 0:0.01:T;
y0=t;
y1=t;
y2=t;

for i=1:length(t)
   y0(i)=L_p_2(t(i));
   y1(i)=L_p_1(t(i)); 
   y2(i)=L_p(t(i)); 
end

plot(t,y0,':k','linewidth',2,'displayname','Lp_u')
plot(t,y1,'r','linewidth',2,'displayname','Lp_a')
plot(t,y2,'g','linewidth',2,'displayname','Lp_m')
legend
 
figure
hold on
for i=1:length(t)
   y0(i)=L_v_2(t(i));
   y1(i)=L_v_1(t(i)); 
   y2(i)=L_v(t(i)); 
end

plot(t,y0,':k','linewidth',2,'displayname','Lv_u')
plot(t,y1,'r','linewidth',2,'displayname','Lv_a')
plot(t,y2,'g','linewidth',2,'displayname','Lv_m')
legend



figure
hold on
for i=1:length(t)
   y0(i)=L_f_2(t(i));
   y1(i)=L_f_1(t(i)); 
   y2(i)=L_f(t(i)); 
end

plot(t,y0,':k','linewidth',2,'displayname','Lf_u')
plot(t,y1,'r','linewidth',2,'displayname','Lf_a')
plot(t,y2,'g','linewidth',2,'displayname','Lf_m')
legend

end

%% put this in mfile: switches the boundary befor/after tbt.
function out = fg_mfile(x,L_p_1,L_p_2,tbt)  
    fval=L_p_1(x);
    gval=L_p_2(x);
    if x>tbt
      out=fval;
    else 
      out= gval;
    end
end

