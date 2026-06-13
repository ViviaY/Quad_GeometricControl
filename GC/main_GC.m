%% clear all caches
clc;           % Clear the command window
% clear all;   % Clear all variables, functions, and class definitions (more thorough)
clearvars;     % Clear variables (safer than clear all)
close all;     % Close all open figure windows
fclose('all'); % Close all open files
rehash;        % Refresh MATLAB's search path and function cache

%% add pathes
% current_dir = fileparts(mfilename('fullpath'));   
current_dir = pwd;
addpath(fullfile(current_dir, 'aux_functions'));
addpath(fullfile(current_dir, 'control'));
addpath(fullfile(current_dir, 'plotting_matlab'));
addpath(fullfile(current_dir, 'gain_tuning'));
addpath(fullfile(current_dir, 'TrajectoryGeneration_TimeVarying'));

% set up the data save path
% root_dir = fileparts(current_dir); % get the parent directory 
save_dir = fullfile(current_dir, 'results');         
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end

% set up the figure save path
fig_dir = fullfile(save_dir, 'Figs');
if ~exist(fig_dir, 'dir')
    mkdir(fig_dir);
end

%% common paramters: 
rng('default');
e3=[0;0;1];
param.g = 9.81;

% % CDC Paper parameters
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


%% Obtaining the control gains and the related opt_bounds
disp("Finding optimal control gains now.")
annealing_output = annealing(anneal_options, param);
param.c1 = annealing_output.opt_bounds.c1;
param.c2 = annealing_output.opt_bounds.c2;
save(fullfile(save_dir, 'annealing_output.mat'), 'annealing_output');

%% obtaining the time-varying Lyapunov bounds
disp("Obtaining the time varying bounds now.")
plotting = false;
Bounds = ObtainingBoundsTimeVarying(annealing_output, anneal_options, param, plotting, save_dir);


%% refrence trajectory generation
disp("Generating a refrence trajectory now.")
generatedTraj = Trajectory_Synthesis_TimeVarying(param, anneal_options, Bounds, save_dir);

% compute the desired trajectory factorial
control_points = generatedTraj.Points_Array(:,:,1);
n_fact = size(control_points,2)-1;
param.factorial_list = [1, factorial(1:n_fact)];


%% get initial points
T = sum(generatedTraj.tau);
fprintf('total travel time of the desired trajectory: %.2f\n', T);
time_step = 1e-2;
t_span = 0:time_step:T;
N = length(t_span); 
init_n = 100; 
delta.x = 0.3;
delta.v = 0.3;
delta.R = 0.5;
delta.W = 1;
show_count = false;

fprintf("Generating %d initial points now \n", init_n);
tic
initial_satisfy = get_initial_points(anneal_options, param, annealing_output, ...
                             delta, init_n, generatedTraj, show_count, save_dir);
elapsed_time1 = toc;
fprintf('Execution time for finding initial_points: %.6f s.\n', elapsed_time1);
initial = initial_satisfy;
case_name = 'satisfy';

%% Generate some points that do not satisfy the initial conditions (for figure 4)
% fprintf("Generating %d unsatisfied initial points now \n", init_n);
% tic
% [initial_outsideBounds, initial_outsideBox] = get_initial_unsatisfypoints(anneal_options, param, annealing_output, delta, init_n, ...
%                                             generatedTraj, show_count, save_dir);
% elapsed_time2 = toc;
% fprintf('Execution time for finding unsatisfied initial_points: %.6f s.\n', elapsed_time2);
% % initial = initial_outsideBounds;
% % case_name = 'outsidebounds';
% initial = initial_outsideBox;
% case_name = 'outsidebox';


%% Generated the tracking trajectory with GC
run_time = zeros(1, init_n);

% saving desired trajectories
dx_list = zeros(3, N);
dv_list = zeros(3, N);
R_list = zeros(initial.init_n, 3, 3, N);
w_list = zeros(initial.init_n, 3, N);

% saving tracking trajectories
x_list = zeros(initial.init_n, 3, N);
v_list = zeros(initial.init_n, 3, N);
f_list = zeros(initial.init_n, N);
Fd3_list = zeros(initial.init_n,N);

% saving the error trajectories
ep_list = zeros(initial.init_n, N);
ev_list = zeros(initial.init_n, N);
eR_list = zeros(initial.init_n, N);
eW_list = zeros(initial.init_n, N);

% saving the Lyapunov results
lyap_list_V = zeros(initial.init_n,N);
lyap_list_V1 = zeros(initial.init_n,N);
lyap_list_V2 = zeros(initial.init_n,N);

for j = 1:init_n

    X0 = initial.initial_points(:, j);
    tic;
    [~, X_sol] = ode45(@(t, XR) eom(t, XR, annealing_output, param, generatedTraj.Points_Array, generatedTraj.tau), t_span, X0, odeset('RelTol', 1e-6, 'AbsTol', 1e-6));
    run_time(j) = toc;
    fprintf("Tracking Trajectories: %d/%d, runtime = %.6f s, total elapsed time = %.6f s, mean elapsed time = %.6f s.\n", ...
                j, init_n, run_time(j), sum(run_time(1:j)), sum(run_time(1:j))/j);

    x_list(j,:,:) = X_sol(:, 1:3)';
    v_list(j,:,:) = X_sol(:, 4:6)';
    w_list(j,:,:) = X_sol(:,7:9)';


    [e, d, R, f, M] = generate_output_arrays(N);
    tic;
    for i = 1:N
        
        R(:,:,i) = reshape(X_sol(i,10:18), 3, 3);
        R_list(j,:,:,i) = R(:,:,i);

        % calculate the desired state and its errors
        des = DesiredTrajectory(t_span(i), generatedTraj.Points_Array, generatedTraj.tau, param);
        [f(i), M(:,i), err, calc] = position_control(X_sol(i,:)', des, annealing_output, param);
       
        d.x(:,i) = des.x;
        d.v(:,i) = des.v;
        d.b1(:,i) = des.b1;
        d.R(:,:,i) = calc.Rd;
        Fd3_list(j,i) = calc.Fd3;
        e.x(:,i) = err.x;
        e.v(:,i) = err.v;
        e.R(:,i) = err.R;
        e.W(:,i) = err.W;

    end
    err_time = toc;
    fprintf("Errors Trajectories: %d/%d, runtime = %.4f s. \n", j, init_n, err_time);
    
    dx_list = d.x;
    dv_list = d.v;
    
    ep_list(j,:) = vecnorm(e.x);
    eW_list(j,:) = vecnorm(e.W);
    ev_list(j,:) = vecnorm(e.v);
    eR_list(j,:) = vecnorm(e.R);

    f_list(j,:) = f;
    
    lyap = get_lyapunov(annealing_output, N, param, e, R, d);
    lyap_list_V(j,:) = lyap.V;
    lyap_list_V1(j,:) = lyap.V1;
    lyap_list_V2(j,:) = lyap.V2;

    clear e d R f M des err calc lyap X_sol
end

avg_time = mean(run_time);    % mean
std_time = std(run_time);     % standard dev
fprintf('Mean runtime = %.4f ± %.4f s\n', avg_time, std_time);

S = struct(); 
for k = 1:initial.init_n
    track_p = squeeze(x_list(k, :, :))';  
    track_v = squeeze(v_list(k, :, :))';  
    
    S.(sprintf('track_p_%d', k)) = track_p;
    S.(sprintf('track_v_%d', k)) = track_v;
end


save(fullfile(save_dir, sprintf("GC_tracking_results_%s.mat", case_name)), "-struct", "S");



%% Plotting trajectories in the environment map
merge_generate_outputs_plots(t_span, generatedTraj, dx_list, x_list, [-130.9964, 67.4], fig_dir);


%% plotting the errors with the bounds, as well as L_v with its bound 
L = zeros(N, 1);
Lpt = zeros(N, 1);
Lvt = zeros(N, 1);
for i = 1:N
    L(i) = Bounds.L(t_span(i))^2;
    Lpt(i) = Bounds.L_p(t_span(i));
    Lvt(i) = Bounds.L_v(t_span(i));
end
crop = 1; 
initial_n = initial.init_n;

save(fullfile(save_dir, sprintf("plot_error_%s.mat", case_name)), ...
             "t_span", "ep_list", "ev_list", "eR_list", "eW_list", ...
             "lyap_list_V", "L", "Lpt", "Lvt", "initial_n");


plot_error(t_span, ep_list, ev_list, lyap_list_V, Lpt, Lvt, L, initial, crop, fig_dir);
plot_eReW(t_span, eR_list, eW_list, initial, crop, fig_dir);

%% tracking velocity and its bounds
plot_v(t_span, v_list, initial, anneal_options.vm, 2, fig_dir);
vm = anneal_options.vm;
initial_n = initial.init_n;
save(fullfile(save_dir, sprintf("plot_v_%s.mat", case_name)), ...
             "t_span", "v_list", "vm", "initial_n");

%% plotting |f| and |F_{d,3}| 
plot_fFd3(t_span, f_list, Fd3_list, initial, annealing_output.opt_bounds, crop, fig_dir);
opt_bounds = annealing_output.opt_bounds;
save(fullfile(save_dir, sprintf("plot_fFd3_%s.mat", case_name)), ...
             "t_span", "f_list", "Fd3_list", "opt_bounds", "initial_n");

%% plotting the data of reference trajectories.
GeneratingPlotsAndOutPuts_TimeVarying(param, anneal_options, generatedTraj, Bounds);
