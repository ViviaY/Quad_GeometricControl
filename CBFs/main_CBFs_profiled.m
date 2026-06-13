%% CBFQP Control for Bezier Curve Trajectory Tracking with LQR
%% Profiled version: same logic as main_CBFs.m, with per-part timing
clc;
clear all;
close all;
fclose('all');
rehash;

%% add paths
import casadi.*
GRAV = 9.81;
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

initial_points_file = fullfile(data_dir, 'InitialPoints.mat');
Initials_data = load(initial_points_file); %#ok<NASGU>

error_file = fullfile(data_dir, 'plot_error_satisfy.mat');
error_data = load(error_file);

[p0_multi, v0_multi, R0_multi, w0_multi] = CBF_InitialPoints(initial_points_file); %#ok<ASGLU>

CBF_delta_type = 'initial';
% CBF_delta_type = 'safety';
% CBF_delta_type = 'reach';
track_tspan = error_data.t_span;

%% set the quadrotor parameters
quadparam = struct( ...
    'g', 9.81, ...
    'm', 4.34, ...
    'J', 1e-2 * diag([8.2, 8.45, 13.77]) ...
);

bezier_traj = GenerateBezierCurves(Bezier_data, track_tspan);

% Build reference from Bezier data
[p_d, v_d, a_d, R_d, f_d, w_d, tau_d, yaw_d] = CBF_build_reference(bezier_traj, quadparam, track_tspan); %#ok<ASGLU>

%% calculate all initial position errors/reaching margin/safety margin as deltas
N_init = size(p0_multi, 1);
delta_all = zeros(N_init, 1);

p0_ref = p_d(1,:);
v0_ref = v_d(1,:); %#ok<NASGU>
dp_infity = zeros(N_init, 1);
dv_infity = zeros(N_init, 1); %#ok<NASGU>

for i = 1:N_init
    pos_err_vec = abs(p0_multi(i,:) - p0_ref);
    pos_err = max(pos_err_vec);
    dp_infity(i) = pos_err;

    switch CBF_delta_type
        case 'safety'
            delta_all(i) = 0.006579;
            save_type = 'safety';
        case 'reach'
            delta_all(i) = 0.0111;
            save_type = 'reach';
        case 'initial'
            delta_all(i) = dp_infity(i);
            save_type = 'initial';
    end
end

%% CBF-QP controller for tracking
cbfqp = CBFQP();

% Set up LQR controller
K_lqr = setup_lqr();

% Simulation parameters
dt = 0.01;
N = length(track_tspan);

fprintf('\n=== CBF simulation settings ===\n');
fprintf('N_init = %d\n', N_init);
fprintf('N steps = %d\n', N);
fprintf('dt = %.6f\n', dt);
fprintf('save_type = %s\n', save_type);
fprintf('================================\n\n');

result_path = fullfile(save_dir, sprintf('CBFs_tracking_%d_%s.mat', N_init, save_type));
all_results = struct();
S = struct();
simulation_times = zeros(N_init, 1);

%% profiling arrays, one total per trajectory
profile_ref        = zeros(N_init, 1);
profile_lqr        = zeros(N_init, 1);
profile_recover    = zeros(N_init, 1);
profile_filter     = zeros(N_init, 1);
profile_virtualize = zeros(N_init, 1);
profile_rk4        = zeros(N_init, 1);
profile_logging    = zeros(N_init, 1);
profile_overhead   = zeros(N_init, 1);
profile_calls      = zeros(N_init, 1);

for ic_idx = 1:N_init

    if mod(ic_idx, 10) == 0 || ic_idx == 1
        fprintf('\n=== Case %d/%d ===\n', ic_idx, N_init);
    end

    %% per-case timing counters
    t_ref        = 0;
    t_lqr        = 0;
    t_recover    = 0;
    t_filter     = 0;
    t_virtualize = 0;
    t_rk4        = 0;
    t_logging    = 0;
    filter_calls = 0;

    %% Initialize state
    z = zeros(6,1);
    z(1:3) = p0_multi(ic_idx,:);
    z(4:6) = v0_multi(ic_idx,:);

    R = squeeze(R0_multi(ic_idx,:,:)); %#ok<NASGU>
    omega = w0_multi(ic_idx,:); %#ok<NASGU>

    %% Pre-allocate histories
    t_hist = zeros(N,1);
    pos_hist = zeros(N,3);
    ref_hist = zeros(N,3);
    err_hist = zeros(N,3);
    margin_hist = zeros(N,3);
    state_hist = zeros(N,6);
    ref_state_hist = zeros(N,9);
    u_hist = zeros(N,4);

    %% Simulation loop timing starts here: same interval as original elapsed=toc
    tStart = tic;

    for k = 1:N
        %% reference lookup + Euler conversion
        tt = tic;

        t = track_tspan(k);

        x_ref = p_d(k,:);
        v_ref_lin = v_d(k,:);

        try
            psi_ref = yaw_d(k);
        catch
            [~, ~, psi_r] = euler_from_R_zyx(R_d{k});
            psi_ref = psi_r;
        end

        x_r = zeros(9,1);
        x_r(1:3) = x_ref;

        [phi_r, theta_r, psi_r] = euler_from_R_zyx(R_d{k});
        x_r(4:6) = [phi_r; theta_r; psi_r];
        x_r(7:9) = v_ref_lin;

        t_ref = t_ref + toc(tt);

        %% LQR nominal control
        tt = tic;

        mu_nom = nominal_lqr_controller(z, x_r, K_lqr);

        t_lqr = t_lqr + toc(tt);

        %% recover: virtual acceleration -> reduced input
        tt = tic;

        v_nom = cbfqp.recover(mu_nom, psi_ref);

        t_recover = t_recover + toc(tt);

        %% CBF-QP safety filter
        tt = tic;

        v_s = cbfqp.filter(z, v_nom, x_r, v_nom, delta_all(ic_idx));
        filter_calls = filter_calls + 1;

        t_filter = t_filter + toc(tt);

        %% virtualize: reduced input -> virtual acceleration
        tt = tic;

        mu_safe = cbfqp.virtualize(v_s);

        t_virtualize = t_virtualize + toc(tt);

        %% RK4 integration
        tt = tic;

        z = rk4_step(z, mu_safe, dt);

        t_rk4 = t_rk4 + toc(tt);

        %% logging
        tt = tic;

        t_hist(k) = t;
        pos_hist(k,:) = z(1:3)';
        ref_hist(k,:) = x_ref;
        err_hist(k,:) = (z(1:3)' - x_ref);
        margin_hist(k,:) = delta_all(ic_idx) - abs(err_hist(k,:));
        state_hist(k,:) = z';
        ref_state_hist(k,:) = x_r';
        u_hist(k,:) = v_s';

        t_logging = t_logging + toc(tt);
    end

    elapsed = toc(tStart);
    simulation_times(ic_idx) = elapsed;

    %% save profiling counters for this case
    profile_sum = t_ref + t_lqr + t_recover + t_filter + t_virtualize + t_rk4 + t_logging;

    profile_ref(ic_idx)        = t_ref;
    profile_lqr(ic_idx)        = t_lqr;
    profile_recover(ic_idx)    = t_recover;
    profile_filter(ic_idx)     = t_filter;
    profile_virtualize(ic_idx) = t_virtualize;
    profile_rk4(ic_idx)        = t_rk4;
    profile_logging(ic_idx)    = t_logging;
    profile_overhead(ic_idx)   = elapsed - profile_sum;
    profile_calls(ic_idx)      = filter_calls;

    %% print profiling for selected cases
    if mod(ic_idx, 10) == 0 || ic_idx == 1 || ic_idx == N_init
        fprintf('\n--- Profiling Case %d ---\n', ic_idx);
        fprintf('Total loop elapsed     = %.6f s\n', elapsed);
        fprintf('Profiled sum           = %.6f s\n', profile_sum);
        fprintf('Unprofiled overhead    = %.6f s\n', elapsed - profile_sum);
        fprintf('Steps                  = %d\n', N);
        fprintf('Filter calls           = %d\n', filter_calls);

        fprintf('ref        total %.6f s | avg %.6f ms/step | %.2f %%\n', ...
            t_ref,        1000*t_ref/N,        100*t_ref/elapsed);
        fprintf('lqr        total %.6f s | avg %.6f ms/step | %.2f %%\n', ...
            t_lqr,        1000*t_lqr/N,        100*t_lqr/elapsed);
        fprintf('recover    total %.6f s | avg %.6f ms/step | %.2f %%\n', ...
            t_recover,    1000*t_recover/N,    100*t_recover/elapsed);
        fprintf('filter     total %.6f s | avg %.6f ms/step | %.2f %%\n', ...
            t_filter,     1000*t_filter/N,     100*t_filter/elapsed);
        fprintf('virtualize total %.6f s | avg %.6f ms/step | %.2f %%\n', ...
            t_virtualize, 1000*t_virtualize/N, 100*t_virtualize/elapsed);
        fprintf('rk4        total %.6f s | avg %.6f ms/step | %.2f %%\n', ...
            t_rk4,        1000*t_rk4/N,        100*t_rk4/elapsed);
        fprintf('logging    total %.6f s | avg %.6f ms/step | %.2f %%\n', ...
            t_logging,    1000*t_logging/N,    100*t_logging/elapsed);
        fprintf('--------------------------\n\n');
    end

    %% Save to all_results
    all_results(ic_idx).t_hist = t_hist;
    all_results(ic_idx).pos_hist = pos_hist;
    all_results(ic_idx).ref_hist = ref_hist;
    all_results(ic_idx).err_hist = err_hist;
    all_results(ic_idx).margin_hist = margin_hist;
    all_results(ic_idx).state_hist = state_hist;
    all_results(ic_idx).ref_state_hist = ref_state_hist;
    all_results(ic_idx).control_hist = u_hist;

    %% Also save to S struct for .mat output, same field names as original
    S.(sprintf('track_p_%d', ic_idx)) = pos_hist;
    S.(sprintf('track_v_%d', ic_idx)) = state_hist(:,4:6);
    S.(sprintf('track_u_%d', ic_idx)) = u_hist;
    S.(sprintf('time_%d', ic_idx)) = elapsed;

    %% keep original behavior: save after each case
    save(result_path, '-struct', 'S');
end

%% profiling summary
valid_calls = max(profile_calls, 1);
fprintf('\n=== PROFILING SUMMARY ===\n');
fprintf('mean ref        = %.6f ms/step\n', 1000 * mean(profile_ref ./ valid_calls));
fprintf('mean lqr        = %.6f ms/step\n', 1000 * mean(profile_lqr ./ valid_calls));
fprintf('mean recover    = %.6f ms/step\n', 1000 * mean(profile_recover ./ valid_calls));
fprintf('mean filter     = %.6f ms/step\n', 1000 * mean(profile_filter ./ valid_calls));
fprintf('mean virtualize = %.6f ms/step\n', 1000 * mean(profile_virtualize ./ valid_calls));
fprintf('mean rk4        = %.6f ms/step\n', 1000 * mean(profile_rk4 ./ valid_calls));
fprintf('mean logging    = %.6f ms/step\n', 1000 * mean(profile_logging ./ valid_calls));
fprintf('mean overhead   = %.6f ms/step\n', 1000 * mean(profile_overhead ./ valid_calls));
fprintf('mean total      = %.6f ms/step\n', 1000 * mean(simulation_times ./ valid_calls));
fprintf('==========================\n\n');

profile_path = fullfile(save_dir, sprintf('CBFs_profile_%d_%s.mat', N_init, save_type));
profile_data = struct();
profile_data.profile_ref        = profile_ref;
profile_data.profile_lqr        = profile_lqr;
profile_data.profile_recover    = profile_recover;
profile_data.profile_filter     = profile_filter;
profile_data.profile_virtualize = profile_virtualize;
profile_data.profile_rk4        = profile_rk4;
profile_data.profile_logging    = profile_logging;
profile_data.profile_overhead   = profile_overhead;
profile_data.profile_calls      = profile_calls;
profile_data.simulation_times   = simulation_times;
profile_data.N                  = N;
profile_data.N_init             = N_init;
save(profile_path, '-struct', 'profile_data');
fprintf('Profiling data saved to %s\n', profile_path);

%% Collect and save error data for MATLAB analysis
cbf_errors = struct();
min_len_track = Inf;

for i = 1:N_init
    result = all_results(i);
    pos_err = result.err_hist;
    vel_err = result.state_hist(:,4:6) - result.ref_state_hist(:,7:9);
    t = result.t_hist;
    min_len_track = min(min_len_track, length(t));

    cbf_errors(i).t = t;
    cbf_errors(i).ep = pos_err;
    cbf_errors(i).ev = vel_err;
end

max_ep_cbf = zeros(min_len_track, 1);
max_ev_cbf = zeros(min_len_track, 1);

for i = 1:min_len_track
    ep_at_i = zeros(N_init, 1);
    ev_at_i = zeros(N_init, 1);

    for j = 1:N_init
        if i <= length(cbf_errors(j).ep)
            ep_at_i(j) = norm(cbf_errors(j).ep(i,:));
            ev_at_i(j) = norm(cbf_errors(j).ev(i,:));
        end
    end

    max_ep_cbf(i) = max(ep_at_i);
    max_ev_cbf(i) = max(ev_at_i);
end

avg_runtime = mean(simulation_times);
std_runtime = std(simulation_times);
fprintf('Mean runtime = %.4f ± %.4f s\n', avg_runtime, std_runtime);

error_data_out = struct();
error_data_out.max_ep_cbf = max_ep_cbf;
error_data_out.max_ev_cbf = max_ev_cbf;
error_data_out.t_span = track_tspan(1:min_len_track);
error_data_out.avg_runtime = avg_runtime;
error_data_out.std_runtime = std_runtime;

for i = 1:N_init
    field_prefix = sprintf('traj_%d_', i);
    error_data_out.(sprintf('%st', field_prefix)) = cbf_errors(i).t(1:min_len_track);
    error_data_out.(sprintf('%sep', field_prefix)) = cbf_errors(i).ep(1:min_len_track,:);
    error_data_out.(sprintf('%sev', field_prefix)) = cbf_errors(i).ev(1:min_len_track,:);
end

error_path = fullfile(save_dir, sprintf('CBFs_error_%d_%s.mat', N_init, save_type));
save(error_path, '-struct', 'error_data_out');
fprintf('\nCBF error data saved to %s\n', error_path);

%% Calculate overall maximum errors
max_pos_inf_norm_overall = 0;
max_pos_2_norm_overall = 0;
max_vel_inf_norm_overall = 0;
max_vel_2_norm_overall = 0;

for i = 1:N_init
    pos_err = cbf_errors(i).ep;
    vel_err = cbf_errors(i).ev;

    max_pos_inf_temp = max(max(abs(pos_err), [], 2));
    max_pos_2_temp = max(sqrt(sum(pos_err.^2, 2)));
    max_vel_inf_temp = max(max(abs(vel_err), [], 2));
    max_vel_2_temp = max(sqrt(sum(vel_err.^2, 2)));

    max_pos_inf_norm_overall = max(max_pos_inf_norm_overall, max_pos_inf_temp);
    max_pos_2_norm_overall = max(max_pos_2_norm_overall, max_pos_2_temp);
    max_vel_inf_norm_overall = max(max_vel_inf_norm_overall, max_vel_inf_temp);
    max_vel_2_norm_overall = max(max_vel_2_norm_overall, max_vel_2_temp);
end

fprintf('\n=== OVERALL ERROR ANALYSIS ===\n');
fprintf('Maximum position errors across all trajectories:\n');
fprintf('  Infinity norm: %.4f m\n', max_pos_inf_norm_overall);
fprintf('  2-norm: %.4f m\n', max_pos_2_norm_overall);
fprintf('\nMaximum velocity errors across all trajectories:\n');
fprintf('  Infinity norm: %.4f m/s\n', max_vel_inf_norm_overall);
fprintf('  2-norm: %.4f m/s\n', max_vel_2_norm_overall);

fprintf('\n=== SIMULATION COMPLETE ===\n');
fprintf('Tested %d initial conditions\n', N_init);
fprintf('Generated plots:\n');
fprintf('- tracking_error_norms.png: Position and velocity error analysis over time\n');
fprintf('- bezier_tracking_3d_multi.png: 3D trajectory comparison\n');

fprintf('Metric, Value\n');
fprintf('max_i max_t ||e_p||_inf,%.4f\n', max_pos_inf_norm_overall);
fprintf('max_i max_t ||e_p||_2,%.4f\n', max_pos_2_norm_overall);
fprintf('max_i max_t ||e_v||_inf,%.4f\n', max_vel_inf_norm_overall);
fprintf('max_i max_t ||e_v||_2,%.4f\n', max_vel_2_norm_overall);
fprintf('mean_runtime,%.4f\n', avg_runtime);
fprintf('std_runtime,%.4f\n', std_runtime);

%% Create error plots
figure('Position', [100, 100, 1200, 900]);
cmap = colormap(jet(N_init));

% Position infinity-norm plot
subplot(2,2,1);
hold on;
safety_bound_label = sprintf('Safety Bound $\\delta = %.1f$ m', cbfqp.delta);
delta_max = max(delta_all);
plot([t_hist(1), t_hist(end)], [delta_max, delta_max], 'k--', 'LineWidth', 2, 'DisplayName', safety_bound_label);

for i = 1:N_init
    result = all_results(i);
    t = result.t_hist;
    pos_err = result.err_hist;
    pos_error_inf_norms = max(abs(pos_err), [], 2);
    plot(t, pos_error_inf_norms, 'Color', cmap(i,:), 'LineWidth', 1.5);
end
ylabel('Position Error $||{\bf e}_p||_{\infty}$ (m)', 'Interpreter', 'latex');
title('Position Error - Infinity Norm vs Time');
grid on;
grid minor;

% Position 2-norm plot
subplot(2,2,2);
hold on;
delta_max = max(delta_all);
plot([t_hist(1), t_hist(end)], [delta_max, delta_max], 'k--', 'LineWidth', 2, 'DisplayName', safety_bound_label);

for i = 1:N_init
    result = all_results(i);
    t = result.t_hist;
    pos_err = result.err_hist;
    pos_error_2_norms = sqrt(sum(pos_err.^2, 2));
    plot(t, pos_error_2_norms, 'Color', cmap(i,:), 'LineWidth', 1.5);
end
ylabel('Position Error $||{\bf e}_p||_2$ (m)', 'Interpreter', 'latex');
title('Position Error - Euclidean Norm vs Time');
grid on;
grid minor;

% Velocity infinity-norm plot
subplot(2,2,3);
hold on;
for i = 1:N_init
    result = all_results(i);
    t = result.t_hist;
    state = result.state_hist;
    ref_state = result.ref_state_hist;
    vel_err = state(:,4:6) - ref_state(:,7:9);
    vel_error_inf_norms = max(abs(vel_err), [], 2);
    plot(t, vel_error_inf_norms, 'Color', cmap(i,:), 'LineWidth', 1.5);
end
xlabel('Time (s)');
ylabel('Velocity Error $||{\bf e}_v||_{\infty}$ (m/s)', 'Interpreter', 'latex');
title('Velocity Error - Infinity Norm vs Time');
grid on;
grid minor;
ylim([0, inf]);

% Velocity 2-norm plot
subplot(2,2,4);
hold on;
for i = 1:N_init
    result = all_results(i);
    t = result.t_hist;
    state = result.state_hist;
    ref_state = result.ref_state_hist;
    vel_err = state(:,4:6) - ref_state(:,7:9);
    vel_error_2_norms = sqrt(sum(vel_err.^2, 2));
    plot(t, vel_error_2_norms, 'Color', cmap(i,:), 'LineWidth', 1.5);
end
xlabel('Time (s)');
ylabel('Velocity Error $||{\bf e}_v||_2$ (m/s)', 'Interpreter', 'latex');
title('Velocity Error - Euclidean Norm vs Time');
grid on;
grid minor;

% Save profiling/error plot for reproducibility
plot_path = fullfile(save_dir, sprintf('CBFs_tracking_error_norms_%d_%s.png', N_init, save_type));
saveas(gcf, plot_path);
fprintf('Tracking error plot saved to %s\n', plot_path);
