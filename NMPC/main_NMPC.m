
%% clear all caches
clc;           % Clear the command window
clear all;     % Clear all variables, functions, and class definitions (more thorough)
close all;     % Close all open figure windows
fclose('all'); % Close all open files
rehash;        % Refresh MATLAB's search path and function cache

%% add pathes
import casadi.*
% current_dir = fileparts(mfilename('fullpath'));   
current_dir = pwd;
addpath(fullfile(current_dir, 'utils'));
rng('default');

save_dir = fullfile(current_dir, 'results');         
if ~exist(save_dir, 'dir')
    mkdir(save_dir);
end


%% import files
root_dir = fileparts(current_dir);
data_dir = fullfile(root_dir, 'GC/results_submit');  
trajs_file = fullfile(data_dir, 'GeneratedTrajectoryData.mat');  
Bezier_data = load(trajs_file);

initial_points_file =  fullfile(data_dir, 'InitialPoints.mat'); 
Initials_data = load(initial_points_file);

error_file =  fullfile(data_dir, 'plot_error_satisfy.mat'); 
error_data = load(error_file);

%%  set the quadotor parameters
QuadConfig = struct();
QuadConfig.m      = 4.34;                     % mass
QuadConfig.g      = 9.81;                     % gravity
QuadConfig.J      = diag([0.082, 0.0845, 0.1377]); % inertia matrix
QuadConfig.J_inv  = inv(QuadConfig.J);        % inverse of inertia
QuadConfig.e3     = [0; 0; 1];                % unit vector along z-axis
QuadConfig.hover_thrust = QuadConfig.m * QuadConfig.g;

QuadConfig.vel_max    = [2.0; 2.0; 2.0];
QuadConfig.acc_max    = [1.0; 1.0; 10.0];

QuadConfig.thrust_max = 2.0 * QuadConfig.hover_thrust;
QuadConfig.thrust_min = 0.1 * QuadConfig.hover_thrust;

%% Nonlinear MPC configuration 
NMPCConfig = struct();

% System configuration
NMPCConfig.NX = 12;    % state dimension: [p, v, R, w] with R \in se(3)
NMPCConfig.NU = 4;     % input dimension: [f, tau_x, tau_y, tau_z]

% MPC configuration
NMPCConfig.H  = 25;    % prediction horizon (shortened for faster solving)
NMPCConfig.dt = 0.01;  % time step

% Solver configuration
NMPCConfig.max_iter = 100;   % maximum iterations for the optimizer
NMPCConfig.tol      = 5e-3;  % solver tolerance

% weights
NMPCConfig.Qp  = diag([200, 200, 200]);
NMPCConfig.Qv  = diag([100, 100, 200]);
NMPCConfig.Qxi = diag([1, 1, 1]);
NMPCConfig.Qw  = diag([0.1, 0.1, 0.1]);

NMPCConfig.Q = blkdiag(NMPCConfig.Qp, NMPCConfig.Qv, NMPCConfig.Qxi, NMPCConfig.Qw);   
NMPCConfig.R = diag([1e-3, 1e-3, 1e-3, 1e-3]);   
NMPCConfig.QT = blkdiag(NMPCConfig.Qp*2, NMPCConfig.Qv*2, NMPCConfig.Qxi, NMPCConfig.Qw);

%% import the Bezier curves and transfer to the reference trajectory with track_tspan (same with CG tracking)
track_tspan = error_data.t_span;
ifPlot = false;
Bezier = sample_bezier(Bezier_data, track_tspan, ifPlot);
reference_traj = build_reference(Bezier, track_tspan, QuadConfig);
if size(reference_traj,2) ~= NMPCConfig.NX && size(reference_traj,1) == NMPCConfig.NX

elseif size(reference_traj,2) == NMPCConfig.NX
    reference_traj = reference_traj.';   % make NX×T
else
    error('reference_traj has unexpected size. Expected NX×T.');
end



%% Run NMPC for all initial conditions and save results
% get the same initial position from file
X0_data = Initials_data.initial.initial_points;
X0 = InitialPoints(X0_data);
M = size(X0, 2);                      
T = size(reference_traj, 2);
prob = nmpc_opti(NMPCConfig, QuadConfig);
opti = prob.opti;      

% save the result
trajs_path = fullfile(save_dir, sprintf('NMPC_trajs_%d_H%d.mat', M, NMPCConfig.H));
S = struct();
total_time = 0;

run_time = zeros(1, M);

for case_idx = 1:M
    fprintf('=== Case %d/%d ===\n', case_idx, M);

    x_curr = X0(:, case_idx);                 % NXx1
    X_ws   = repmat(x_curr, 1, NMPCConfig.H+1);          % NXx(H+1)
    U_ws   = zeros(NMPCConfig.NU, NMPCConfig.H);                    % NUxH
    U_ws(1,:)   = QuadConfig.hover_thrust;    % hover thrust
    U_ws(1,:)   = min(max(U_ws(1,:),   QuadConfig.thrust_min), QuadConfig.thrust_max);

    % ---- preallocate executed (closed-loop) trajectory ----
    N_iter = T-1; % - NMPCConfig.H;                           % MPC steps
    X_exec = zeros(NMPCConfig.NX, T);                    % NX x T
    U_exec = zeros(NMPCConfig.NU, T-1);                  % NU x (T-1)
    X_exec(:,1) = x_curr;
    step = 0;

    tStart = tic;

    % ---- receding horizon loop ----
    for t_idx = 1:N_iter
        % set reference window and x0
        % Xref_win = reference_traj(:, t_idx:t_idx+NMPCConfig.H);  % NX x (H+1)
        Xref_win = ref_window_pad(reference_traj, t_idx, NMPCConfig.H, NMPCConfig.dt);
        opti.set_value(prob.Xref, Xref_win);
        opti.set_value(prob.x0,   x_curr);

        % warm-start
        opti.set_initial(prob.X, X_ws);
        prob.opti.set_initial(prob.U, U_ws);

        % solve with fallback
        try
            sol   = opti.solve();
            X_sol = sol.value(prob.X);        % NX x (H+1)
            U_sol = sol.value(prob.U);        % NU x H
        catch ME
            warning('NMPC step %d failed: %s', t_idx, ME.message);
            X_sol = []; 
            U_sol = [];
            try
                D = opti.debug();        % last iterate
                X_sol = D.value(prob.X);
                U_sol = D.value(prob.U);
            catch
            end
            if isempty(X_sol) || any(~isfinite(X_sol(:))), X_sol = X_ws; end
            if isempty(U_sol) || any(~isfinite(U_sol(:))), U_sol = U_ws; end
        end

        % apply first control and log
        u_apply = U_sol(:,1);
        step = step + 1;
        U_exec(:, step) = u_apply;
        x_curr = rk4_step(prob.f_dyn, x_curr, u_apply, NMPCConfig.dt);
        X_exec(:, step+1) = x_curr;

        % shift warm-start (hold last column)
        X_ws = [X_sol(:,2:end), X_sol(:,end)];
        U_ws = [U_sol(:,2:end), U_sol(:,end)];
    end


    elapsed = toc(tStart);
    run_time(case_idx) = elapsed;
    total_time = total_time + elapsed;
    fprintf('NMPC loop: total = %.3f s, avg per step = %.3f ms, steps = %d\n', elapsed, 1e3*elapsed/step, step);

    track_p  = X_exec(1:3,   :).';           % T x 3
    track_v  = X_exec(4:6,   :).';           % T x 3
    track_xi = X_exec(7:9,   :).';           % T x 3
    track_w  = X_exec(10:12, :).';           % T x 3

    track_R = zeros(T,3,3);
    for t = 1:T
        % numeric branch; project to SO(3) for robustness
        track_R(t,:,:) = so3_exp(X_exec(7:9, t));
    end

    track_u = U_exec.';                       % (T-1) x 4


    S.(sprintf('track_p_%d', case_idx)) = track_p;
    S.(sprintf('track_v_%d', case_idx)) = track_v;
    S.(sprintf('track_R_%d', case_idx)) = track_R;
    S.(sprintf('track_w_%d', case_idx)) = track_w;
    S.(sprintf('track_u_%d', case_idx)) = track_u;
    S.(sprintf('time_%d',    case_idx)) = elapsed;

    % save after each case (safer for long runs)
    save(trajs_path, '-struct', 'S');
end

fprintf('All %d cases finished. Saved to:\n%s\n', M, trajs_path);
avg_time = mean(run_time);    % mean
std_time = std(run_time);     % standard dev
fprintf('Mean runtime = %.4f ± %.4f s\n', avg_time, std_time);


%% === Collect and save NMPC error data (same format as CBF) ===
nmpc_errors = struct();
min_len_track = Inf;

for i = 1:M

    track_p = S.(sprintf('track_p_%d', i));     
    track_v = S.(sprintf('track_v_%d', i));     
    T_i = size(track_p, 1);

    ref_p = reference_traj(1:3, 1:T_i).';       
    ref_v = reference_traj(4:6, 1:T_i).';       
    t = (0:T_i-1)' * NMPCConfig.dt;             

    pos_err = track_p - ref_p;                  
    vel_err = track_v - ref_v;                

    nmpc_errors(i).t = t;
    nmpc_errors(i).ep = pos_err;
    nmpc_errors(i).ev = vel_err;

    min_len_track = min(min_len_track, T_i);
end

max_ep_nmpc = zeros(min_len_track, 1);
max_ev_nmpc = zeros(min_len_track, 1);

for i = 1:min_len_track
    ep_at_i = zeros(M, 1);
    ev_at_i = zeros(M, 1);
    
    for j = 1:M
        if i <= length(nmpc_errors(j).ep)
            ep_at_i(j) = norm(nmpc_errors(j).ep(i,:));
            ev_at_i(j) = norm(nmpc_errors(j).ev(i,:));
        end
    end
    
    max_ep_nmpc(i) = max(ep_at_i);
    max_ev_nmpc(i) = max(ev_at_i);
end

error_data = struct();
error_data.max_ep_nmpc = max_ep_nmpc;
error_data.max_ev_nmpc = max_ev_nmpc;
error_data.t_span = (0:min_len_track-1)' * NMPCConfig.dt;
error_data.avg_runtime = mean(run_time);
error_data.std_runtime = std(run_time);


for i = 1:M
    field_prefix = sprintf('traj_%d_', i);
    error_data.(sprintf('%st', field_prefix))  = nmpc_errors(i).t(1:min_len_track);
    error_data.(sprintf('%sep', field_prefix)) = nmpc_errors(i).ep(1:min_len_track,:);
    error_data.(sprintf('%sev', field_prefix)) = nmpc_errors(i).ev(1:min_len_track,:);
end

errors_path = fullfile(save_dir, sprintf('NMPC_error_%d_H%d.mat', M, NMPCConfig.H));
save(errors_path, '-struct', 'error_data');
fprintf('\nNMPC errors have been saved\n');

%% plotting
VA = track_v(end-NMPCConfig.H:end,:) - reference_traj(4:6, end-NMPCConfig.H:end)';
verr_row = vecnorm(VA, 2, 2);

figure;
L = min(numel(track_tspan), size(track_v,1));     
t = track_tspan(1:L);
v = track_v(1:L, :);               
plot(track_tspan, v(:,1), track_tspan, v(:,2), track_tspan, v(:,3), 'LineWidth',1.2);
legend('v_x','v_y','v_z','Location','best');
title('tracking velocity vs time');

figure;
v_d = reference_traj(4:6, 1:L);
plot(track_tspan, v_d(1,:), track_tspan, v_d(2,:), track_tspan, v_d(3,:), 'LineWidth',1.2);
title('desired velocity vs time');
legend('v_x','v_y','v_z','Location','best');

figure;
e_v = abs(v' - v_d);
plot(track_tspan, e_v(1,:), track_tspan, e_v(2,:), track_tspan, e_v(3,:), 'LineWidth',1.2);
 legend('ev_x','ev_y','ev_z','Location','best');
title('desired velocity vs time');
