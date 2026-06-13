%% clear all caches
clc;           % Clear the command window
clear all;     % Clear all variables, functions, and class definitions (more thorough)
close all;     % Close all open figure windows
fclose('all'); % Close all open files
rehash;        % Refresh MATLAB's search path and function cache

% add pathes
addpath('aux_functions');
addpath('control');
addpath('plotting_matlab');
addpath('gain_tuning');
addpath('TrajectoryGeneration_TimeVarying');
userparam = struct();
rng('default'); % for RRT 

p = haltonset(3,'Skip',1000,'Leap',50);
p = scramble(p,'RR2');
samples = net(p,120);   % generate 120 samples
alpha_seq    = samples(:,1);
alpha_s_seq  = samples(:,2);
Csample_seq  = samples(:,3);
% common paramters: 
e3=[0;0;1];
param.g = 9.81;

param.m = 4.34; 
param.J = diag([0.0820, 0.0845, 0.1377]);
param.J_min = min(diag(param.J));
param.J_max = max(diag(param.J)); 


anneal_options.initial_solution = [10, 10, 10, 10, 0.5, 0.5];
anneal_options.lower_bound = [0.1, 0.1, 0.1, 0.1, 0.1, 0.1];
anneal_options.upper_bound = [30, 30, 30, 30, 1, 1];
anneal_options.runtime_n = 2000;
anneal_options.vm = [2;2;2];
anneal_options.am = [1;1;10];
anneal_options.psi_bar = 0.005; 
anneal_options.alpha_psi = 0.4;
anneal_options.V1_0 = 0.4; 
anneal_options.w1 = 15;
anneal_options.w2 = 1;
anneal_options.w3 = 1;

annealing_output.opt_k = [18.3154    6.4322   22.3572    1.4732    0.5651    0.5930];
annealing_output.opt_bounds = evaluate(annealing_output.opt_k, anneal_options, param);
param.c1 = annealing_output.opt_bounds.c1;
param.c2 = annealing_output.opt_bounds.c2;
plotting = false;
Bounds = ObtainingBoundsTimeVarying(annealing_output, anneal_options, param, plotting, []);

quadparam = struct(...
    'g', 9.81,...
    'm', 4.34, ...  % kg
    'J', 1e-2*diag([8.2, 8.45, 13.77]) ...  % moment of inertia
);

%% alpha
diary('alpha_experiment_log_norand.txt');
diary on;
% %rng('default');
RRT_times  = zeros(100,1);
ILP_times  = zeros(100,1);
total_times = zeros(100,1);
Ns_list = zeros(100,1);
alpha_v_list   = zeros(100,1);
valid_count = 0;     
infeasible_num = 0;    
idx = 1;
userparam.name = "alpha";

fprintf('************* checking alpha *********** \n');
while valid_count < 100
    % %rng('default');
    % rng('shuffle');
    
    alpha_k = 0.4 + 0.6 * alpha_seq(idx);
%     alpha_k = alpha_seq(idx);
    idx = idx + 1;  % Increment index for the next alpha_s_k generation
    userparam.value = alpha_k;
%     alpha_k = rand();

    generatedTraj = Trajectory_Synthesis_TimeVarying_checking(param, anneal_options, Bounds, userparam);

    if isempty(generatedTraj) || ~isstruct(generatedTraj)
        infeasible_num = infeasible_num + 1;
        fprintf('[Infeasible count %d] infeasible → skipped, alpha_k=%.3f, \n', infeasible_num, alpha_k);
        continue;
    end

    valid_count = valid_count + 1;

    T = generatedTraj.t_Array(end);

    fprintf(['[run %d] CPU RRT = %.3fs | CPU ILP = %.3fs | Travel = %.2f | ', ...
             'Ns=%d | Nv=%d | eps=%.1e | α=%.2f | α_s=%.2f | α_v=%.2f\n'], ...
             valid_count, generatedTraj.RRT_time, generatedTraj.Bezier_ILPtime, ...
             T, generatedTraj.Ns, generatedTraj.Nv, ...
             generatedTraj.epsilon, generatedTraj.alpha, ...
             generatedTraj.alpha_s, generatedTraj.alpha_v);

    RRT_times(valid_count)    = generatedTraj.RRT_time;
    ILP_times(valid_count)    = generatedTraj.Bezier_ILPtime;
    total_times(valid_count)  = RRT_times(valid_count) + ILP_times(valid_count);
    Ns_list(valid_count)      = generatedTraj.Ns;
    alpha_v_list(valid_count)   = generatedTraj.alpha_v;
end


mu_RRT = mean(RRT_times);  sd_RRT = std(RRT_times);
mu_ILP = mean(ILP_times);  sd_ILP = std(ILP_times);
mu_total = mean(total_times);  sd_total = std(total_times);
mu_Ns = mean(Ns_list);  sd_Ns = std(Ns_list(1:100));
mu_alpha_v = mean(alpha_v_list); sd_alpha_v = std(alpha_v_list);

fprintf('\n===== Summary over %d successful runs =====\n', valid_count);
fprintf('RRT time   = %.2f ± %.2f seconds\n', mu_RRT, sd_RRT);
fprintf('ILP time   = %.2f ± %.2f seconds\n', mu_ILP, sd_ILP);
fprintf('Total time = %.2f ± %.2f seconds\n', mu_total, sd_total);
fprintf('Ns         = %.2f ± %.2f \n', mu_Ns, sd_Ns);
fprintf('alpha_v = %.2f ± %.2f\n', mu_alpha_v, sd_alpha_v);
diary off;






%% alpha_s
diary('alpha_s_experiment_log.txt');
diary on;
% %rng('default');
RRT_times  = zeros(100,1);
ILP_times  = zeros(100,1);
total_times = zeros(100,1);
Ns_list = zeros(100,1);
alpha_v_list   = zeros(100,1);

valid_count = 0;     
infeasible_num = 0;     
idx = 1;
userparam.name = "alpha_s";
fprintf("************* checking alpha_s: alpha_s_k = 0.5+0.45*alpha_s_seq(idx); *********** \n");
while valid_count < 100

    alpha_s_k = 0.5+0.45*alpha_s_seq(idx);
    idx = idx + 1;
% if alpha_s_seq(idx) < 0.95 
%     alpha_s_k = alpha_s_seq(idx);
%     idx = idx + 1;
% else
%     idx = idx +1;
%     continue
% end
%     alpha_s_k = rand();
    userparam.value = alpha_s_k;


    generatedTraj = Trajectory_Synthesis_TimeVarying_checking(param, anneal_options, Bounds, userparam);


    if isempty(generatedTraj) || ~isstruct(generatedTraj)
        infeasible_num = infeasible_num + 1;
        fprintf('[Infeasible count %d] infeasible → skipped, alpha_s_k=%.3f, \n', infeasible_num, alpha_s_k);
        continue;
    end


    valid_count = valid_count + 1;


    T = generatedTraj.t_Array(end);

    fprintf(['[run %d] CPU RRT = %.3fs | CPU ILP = %.3fs | Travel = %.2f | ', ...
             'Ns=%d | Nv=%d | eps=%.1e | α=%.2f | α_s=%.2f | α_v=%.2f\n'], ...
             valid_count, generatedTraj.RRT_time, generatedTraj.Bezier_ILPtime, ...
             T, generatedTraj.Ns, generatedTraj.Nv, ...
             generatedTraj.epsilon, generatedTraj.alpha, ...
             generatedTraj.alpha_s, generatedTraj.alpha_v);

    RRT_times(valid_count)    = generatedTraj.RRT_time;
    ILP_times(valid_count)    = generatedTraj.Bezier_ILPtime;
    total_times(valid_count)  = RRT_times(valid_count) + ILP_times(valid_count);
    Ns_list(valid_count)      = generatedTraj.Ns;
    alpha_v_list(valid_count)   = generatedTraj.alpha_v; 
end



mu_RRT = mean(RRT_times);  sd_RRT = std(RRT_times);
mu_ILP = mean(ILP_times);  sd_ILP = std(ILP_times);
mu_total = mean(total_times);  sd_total = std(total_times);
mu_Ns = mean(Ns_list);  sd_Ns = std(Ns_list);
mu_alpha_v = mean(alpha_v_list);  sd_alpha_v = std(alpha_v_list);

fprintf('\n===== Summary over %d successful runs =====\n', valid_count);
fprintf('RRT time   = %.2f ± %.2f seconds\n', mu_RRT, sd_RRT);
fprintf('ILP time   = %.2f ± %.2f seconds\n', mu_ILP, sd_ILP);
fprintf('Total time = %.2f ± %.2f seconds\n', mu_total, sd_total);
fprintf('Ns         = %.2f ± %.2f \n', mu_Ns, sd_Ns);
fprintf('alpha_v         = %.2f ± %.2f \n', mu_alpha_v, sd_alpha_v);
diary off;


%% C_sampling
diary('C_sampling_experiment_log.txt');
diary on;
% %rng('default');
RRT_times  = zeros(100,1);
ILP_times  = zeros(100,1);
total_times = zeros(100,1);
Ns_list = zeros(100,1);
alpha_v_list   = zeros(100,1);

valid_count = 0;    
infeasible_num = 0;    
idx = 1;
userparam.name = "C_sample";
fprintf("************* checking C_sample：Csampling_k = 0.7+0.3*Csample_seq(idx);*********** \n");
while valid_count < 100

    %Csampling_k = Csample_seq(idx);
    Csampling_k = 0.7+0.3*Csample_seq(idx);
    idx = idx + 1;
    userparam.value = Csampling_k;

    generatedTraj = Trajectory_Synthesis_TimeVarying_checking(param, anneal_options, Bounds, userparam);


    if isempty(generatedTraj) || ~isstruct(generatedTraj)
        infeasible_num = infeasible_num + 1;
        fprintf('[Infeasible count %d] infeasible → skipped, Csampling_k=%.3f, \n', infeasible_num, Csampling_k);
        continue;
    end

    if Csampling_k <= 0.2
        generate_outputs_plots(generatedTraj, [-130.9964, 67.4]);
    end


    valid_count = valid_count + 1;


    T = generatedTraj.t_Array(end);

    fprintf(['[run %d] CPU RRT = %.3fs | CPU ILP = %.3fs | Travel = %.2f | ', ...
             'Ns=%d | Nv=%d | eps=%.1e | α=%.2f | α_s=%.2f | α_v=%.2f | C_sampling=%0.2f\n'], ...
             valid_count, generatedTraj.RRT_time, generatedTraj.Bezier_ILPtime, ...
             T, generatedTraj.Ns, generatedTraj.Nv, ...
             generatedTraj.epsilon, generatedTraj.alpha, ...
             generatedTraj.alpha_s, generatedTraj.alpha_v, Csampling_k);


    RRT_times(valid_count)    = generatedTraj.RRT_time;
    ILP_times(valid_count)    = generatedTraj.Bezier_ILPtime;
    total_times(valid_count)  = RRT_times(valid_count) + ILP_times(valid_count);
    Ns_list(valid_count)      = generatedTraj.Ns;
    alpha_v_list(valid_count)   = generatedTraj.alpha_v; 
end


mu_RRT = mean(RRT_times);  sd_RRT = std(RRT_times);
mu_ILP = mean(ILP_times);  sd_ILP = std(ILP_times);
mu_total = mean(total_times);  sd_total = std(total_times);
mu_Ns = mean(Ns_list);  sd_Ns = std(Ns_list);
mu_alpha_v = mean(alpha_v_list);  sd_alpha_v = std(alpha_v_list);

fprintf('\n===== Summary over %d successful runs =====\n', valid_count);
fprintf('RRT time   = %.2f ± %.2f seconds\n', mu_RRT, sd_RRT);
fprintf('ILP time   = %.2f ± %.2f seconds\n', mu_ILP, sd_ILP);
fprintf('Total time = %.2f ± %.2f seconds\n', mu_total, sd_total);
fprintf('Ns         = %.2f ± %.2f \n', mu_Ns, sd_Ns);
fprintf('alpha_v         = %.2f ± %.2f \n', mu_alpha_v, sd_alpha_v);
diary off;


%% N_p
diary('N_p_experiment_log.txt');
diary on;
% rng('default');
rng(0)
RRT_times  = zeros(100,1);
ILP_times  = zeros(100,1);
total_times = zeros(100,1);
Ns_list = zeros(100,1);
alpha_v_list   = zeros(100,1);


userparam.name = "Np";

for Np = 8:25
    fprintf('************* checking Np = %d *********** \n', Np);
    userparam.value = Np;
    valid_count = 0;    
    infeasible_num = 0;    
    while valid_count < 100
    generatedTraj = Trajectory_Synthesis_TimeVarying_checking(param, anneal_options, Bounds, userparam);

    if isempty(generatedTraj) || ~isstruct(generatedTraj)
        infeasible_num = infeasible_num + 1;
        fprintf('[Infeasible count %d] infeasible → skipped, Np=%1f, \n', infeasible_num, Np);
        continue;
    end

    valid_count = valid_count + 1;

    T = generatedTraj.t_Array(end);

    fprintf(['[run %d] CPU RRT = %.3fs | CPU ILP = %.3fs | Travel = %.2f | ', ...
             'Ns=%d | Nv=%d | Np=%d | eps=%.1e | α=%.2f | α_s=%.2f | α_v=%.2f\n'], ...
             valid_count, generatedTraj.RRT_time, generatedTraj.Bezier_ILPtime, ...
             T, generatedTraj.Ns, generatedTraj.Nv, generatedTraj.Np, ...
             generatedTraj.epsilon, generatedTraj.alpha, ...
             generatedTraj.alpha_s, generatedTraj.alpha_v);

    RRT_times(valid_count)    = generatedTraj.RRT_time;
    ILP_times(valid_count)    = generatedTraj.Bezier_ILPtime;
    total_times(valid_count)  = RRT_times(valid_count) + ILP_times(valid_count);
    Ns_list(valid_count)      = generatedTraj.Ns;
    alpha_v_list(valid_count)   = generatedTraj.alpha_v; 

    end

mu_RRT  = mean(RRT_times(1:valid_count));
sd_RRT  = std(RRT_times(1:valid_count));

mu_ILP  = mean(ILP_times(1:valid_count));
sd_ILP  = std(ILP_times(1:valid_count));

mu_total = mean(total_times(1:valid_count));
sd_total = std(total_times(1:valid_count));

mu_Ns  = mean(Ns_list(1:valid_count));
sd_Ns  = std(Ns_list(1:valid_count));

mu_alpha_v = mean(alpha_v_list(1:valid_count));
sd_alpha_v = std(alpha_v_list(1:valid_count));


fprintf('\n=====  Np = %d, Summary over %d successful runs =====\n', Np, valid_count);
fprintf('RRT time   = %.2f ± %.2f seconds\n', mu_RRT, sd_RRT);
fprintf('ILP time   = %.2f ± %.2f seconds\n', mu_ILP, sd_ILP);
fprintf('Total time = %.2f ± %.2f seconds\n', mu_total, sd_total);
fprintf('Ns         = %.2f ± %.2f \n', mu_Ns, sd_Ns);
fprintf('alpha_v         = %.2f ± %.2f \n', mu_alpha_v, sd_alpha_v);

end
diary off;

%% The defined pamaters setting in the paper
userparam.name = 'Nothing';
generatedTraj = Trajectory_Synthesis_TimeVarying_checking(param, anneal_options, Bounds, userparam);
    fprintf([' CPU RRT = %.3fs | CPU ILP = %.3fs | Travel = %.2f | ', ...
             'Ns=%d | Nv=%d | Np=%d | eps=%.1e | α=%.2f | α_s=%.2f | α_v=%.2f\n'], ...
             generatedTraj.RRT_time, generatedTraj.Bezier_ILPtime, ...
             T, generatedTraj.Ns, generatedTraj.Nv, generatedTraj.Np, ...
             generatedTraj.epsilon, generatedTraj.alpha, ...
             generatedTraj.alpha_s, generatedTraj.alpha_v);