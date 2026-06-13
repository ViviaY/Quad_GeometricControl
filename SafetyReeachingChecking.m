%% clear all caches
clc;           % Clear the command window
clear all;     % Clear all variables, functions, and class definitions (more thorough)
close all;     % Close all open figure windows
fclose('all'); % Close all open files
rehash;        % Refresh MATLAB's search path and function cache

%% add pathes
addpath('GC/results_submit')
addpath('NMPC/results_submit')
addpath('CBFs/results_submit')

GC_results = load('GC_tracking_results_satisfy.mat');

NMPC_H10 = load('NMPC_trajs_100_H10.mat');
NMPC_H15 = load('NMPC_trajs_100_H15.mat');
NMPC_H20 = load('NMPC_trajs_100_H20.mat');
NMPC_H25 = load('NMPC_trajs_100_H25.mat');

CBF_initerrs = load('CBFs_tracking_100_initial.mat');
CBF_reachmargin = load('CBFs_tracking_100_reach.mat');
CBF_safetymargin = load('CBFs_tracking_100_safety.mat');
num = 100;
%% import files
trajs_file = 'GeneratedTrajectoryData.mat';
Bezier_data = load(trajs_file);
error_matlab_file = 'plot_error_satisfy.mat';
% Load error data for time span
error_data = load(error_matlab_file);
track_tspan = error_data.t_span;

quadparam = struct(...
    'g', 9.81,...
    'm', 4.34, ...  % kg
    'J', 1e-2*diag([8.2, 8.45, 13.77]) ...  % moment of inertia
);
bezier_traj = GenerateBezierCurves(Bezier_data, track_tspan);
% Build reference from Bezier data
[p_d, v_d, a_d, R_d, f_d, w_d, tau_d, yaw_d] = CBF_build_reference(bezier_traj, quadparam, track_tspan);


%% set up the environments
Xsl = Bezier_data.Xsl;
Xsu = Bezier_data.Xsu;
Cs=0.5*(Xsl+Xsu);
Ds=0.5*(Xsu-Xsl);

XulArray = Bezier_data.XulArray;
XuuArray = Bezier_data.XuuArray;

CuArray=0.5*(XulArray+XuuArray);
DuArray=0.5*(XuuArray-XulArray);

Xtl = Bezier_data.Xtl;
Xtu = Bezier_data.Xtu;
Ct=0.5*(Xtl+Xtu);
Dt=0.5*(Xtu-Xtl);

L_p = Bezier_data.L_p;
time_t = Bezier_data.t_Array(end);
L_pt = L_p(time_t);
Margin_target = [L_pt;L_pt;L_pt];


%% check all trajectories
[success,unsafety_num,notReaching_num] = VerifyRefTrajectories(p_d, Cs, Ds, CuArray, DuArray, Ct, Dt, Margin_target);
target_obs = [2, 6];
num_pts = size(p_d, 2);

for obs = target_obs
    center = CuArray(:, obs);    % 3×1
    d = DuArray(:, obs);         % 3×1

    pos_traj = p_d(:,1:3)';              % 3×T
    center   = CuArray(:, obs);          % 3×1
    d        = DuArray(:, obs);          % 3×1
    
    delta_abs = abs(pos_traj - center);  % 3×T
    s_inf_all = min(d - delta_abs, [], 1);   % 1×T (signed infinity-norm margin: positive inside, negative outside; zero on the boundary)
    

    [~, min_idx] = min(abs(s_inf_all));
    s_inf_closest = s_inf_all(min_idx);
    
    fprintf('Obstacle %d: Closest trajectory index = %d, s_inf = %.6f\n', ...
            obs, min_idx, s_inf_closest);
    fprintf('p_d(:,%d) = [%.6f %.6f %.6f]\n', min_idx, pos_traj(:, min_idx));
    fprintf('CuArray(:,%d) = [%.6f %.6f %.6f]\n', obs, center);
    fprintf('DuArray(:,%d) = [%.6f %.6f %.6f]\n', obs, d);
end


%%
[safety_rate,unsafety_rate,notReaching_rate] = VerifyTrajectories(GC_results, num, Cs, Ds, CuArray, DuArray, Ct, Dt, Margin_target);
fprintf("GC: successful rate: %.2f%%, Unsafety rate: %.2f%%, Not reaching rate: %.2f%%\n", ...
        safety_rate, unsafety_rate, notReaching_rate);


%%

[safety_rate,unsafety_rate,notReaching_rate] = VerifyTrajectories(NMPC_H10, num, Cs, Ds, CuArray, DuArray, Ct, Dt, Margin_target);
fprintf("NMPC H=10: successful rate: %.2f%%, Unsafety rate: %.2f%%, Not reaching rate: %.2f%%\n", ...
        safety_rate, unsafety_rate, notReaching_rate);

[safety_rate,unsafety_rate,notReaching_rate] = VerifyTrajectories(NMPC_H15, num, Cs, Ds, CuArray, DuArray, Ct, Dt, Margin_target);
fprintf("NMPC H=15: successful rate: %.2f%%, Unsafety rate: %.2f%%, Not reaching rate: %.2f%%\n", ...
        safety_rate, unsafety_rate, notReaching_rate);

[safety_rate,unsafety_rate,notReaching_rate] = VerifyTrajectories(NMPC_H20, num, Cs, Ds, CuArray, DuArray, Ct, Dt, Margin_target);
fprintf("NMPC H=20: successful rate: %.2f%%, Unsafety rate: %.2f%%, not reaching rate: %.2f%%\n", ...
        safety_rate, unsafety_rate, notReaching_rate);

[safety_rate,unsafety_rate,notReaching_rate] = VerifyTrajectories(NMPC_H25, num, Cs, Ds, CuArray, DuArray, Ct, Dt, Margin_target);
fprintf("NMPC H=25: successful rate: %.2f%%, Unsafety rate: %.2f%%, not reaching rate: %.2f%%\n", ...
        safety_rate, unsafety_rate, notReaching_rate);

[safety_rate,unsafety_rate,notReaching_rate] = VerifyTrajectories(GC_results, num, Cs, Ds, CuArray, DuArray, Ct, Dt, Margin_target);
fprintf("GC: successful rate: %.2f%%, Unsafety rate: %.2f%%, Not reaching rate: %.2f%%\n", ...
        safety_rate, unsafety_rate, notReaching_rate);

[safety_rate,unsafety_rate,notReaching_rate] = VerifyTrajectories(CBF_initerrs, num, Cs, Ds, CuArray, DuArray, Ct, Dt, Margin_target);
fprintf("CBFs_initialerrors_results: successful rate: %.2f%%, Unsafety rate: %.2f%%, Not reaching rate: %.2f%%\n", ...
        safety_rate, unsafety_rate, notReaching_rate);

[safety_rate,unsafety_rate,notReaching_rate] = VerifyTrajectories(CBF_reachmargin, num, Cs, Ds, CuArray, DuArray, Ct, Dt, Margin_target);
fprintf("CBFs_reachingmargin_results: successful rate: %.2f%%, Unsafety rate: %.2f%%, Not reaching rate: %.2f%%\n", ...
        safety_rate, unsafety_rate, notReaching_rate);

[safety_rate,unsafety_rate,notReaching_rate] = VerifyTrajectories(CBF_safetymargin, num, Cs, Ds, CuArray, DuArray, Ct, Dt, Margin_target);
fprintf("CBFs_safetymargin_results: successful rate: %.2f%%, Unsafety rate: %.2f%%, Not reaching rate: %.2f%%\n", ...
        safety_rate, unsafety_rate, notReaching_rate);
%%
function [success,unsafety_num,notReaching_num] = VerifyRefTrajectories(p_d, Cs, Ds, CuArray, DuArray, Ct, Dt, Margin_target)
    success = 0;
    unsafety_num = 0;
    notReaching_num = 0;
    % safety_flags_all = cell(1, num);

    T = size(p_d, 1);  
    safety_flags = zeros(1, T); 

    % check the safety
    for k = 1:T
        safety_flags(k) = SafetyChecking(p_d(k,:)', Cs, Ds, CuArray, DuArray, k, false);
    end
    all_safe = all(safety_flags);

    % check if reaching the target
    dist_T = min(Dt - Margin_target - abs(p_d(end,:)' - Ct));
    fprintf("target margin: %.4f\n", dist_T);


    if all_safe && dist_T > 0
        success = success + 1;
        % fprintf('Trajectory %d is safe and reached the target. \n', j);
    end
    if ~all_safe
        unsafe_idx = find(safety_flags==0);
        unsafety_num = unsafety_num +1;
        % fprintf('the referencre trajectory is unsafe at indices: %s\n', mat2str(unsafe_idx));
    end

    if dist_T <= 0
        notReaching_num = notReaching_num + 1;
        % fprintf('Trajectory %d did not reach the target! \n', j);
    end

    % safety_flags_all{j} = safety_flags;
end 

%% 
function [success_rate,unsafety_rate,notReaching_rate] = VerifyTrajectories(trajs, num, Cs, Ds, CuArray, DuArray, Ct, Dt, Margin_target)
    success = 0;
    unsafety_num = 0;
    notReaching_num = 0;
    safety_flags_all = cell(1, num);

    for j = 1:num
        fieldname = sprintf('track_p_%d', j);   
        X_traj = trajs.(fieldname);      

        T = size(X_traj, 1);  
        safety_flags = zeros(1, T); 

        % check the safety
        for k = 1:T
            safety_flags(k) = SafetyChecking(X_traj(k, :)', Cs, Ds, CuArray, DuArray, k, false);
        end
        all_safe = all(safety_flags);

        % check if reaching the target
        dist_T = min(Dt - Margin_target - abs(X_traj(end,:)' - Ct));

        if all_safe && dist_T > 0
            success = success + 1;
            % fprintf('Trajectory %d is safe and reached the target. \n', j);
        end
        if ~all_safe
            unsafe_idx = find(safety_flags==0);
            unsafety_num = unsafety_num +1;
            % fprintf('Trajectory %d is unsafe at indices: %s\n', j, mat2str(unsafe_idx));
            for i = 1:length(unsafe_idx)
                x_i = X_traj(unsafe_idx(i), :)';
                % fprintf('Checking unsafe index %d\n', unsafe_idx(i));
                SafetyChecking(x_i, Cs, Ds, CuArray, DuArray, unsafe_idx(i), false);
            end
        end

        if dist_T <= 0
            notReaching_num = notReaching_num + 1;
            % fprintf('Trajectory %d did not reach the target! \n', j);
        end

        safety_flags_all{j} = safety_flags;
    end 
    success_rate = 100 * success / num;
    unsafety_rate = 100 * unsafety_num / num;
    notReaching_rate = 100 * notReaching_num / num;
end


%%
function [safetycheck] = SafetyChecking(Xc,Cs,Ds,CuArray,DuArray, index, verb)
Nu = size(CuArray,2);
safetycheck = 1;

[val1, idx1] = min(Ds - abs(Xc - Cs));
tol = 1e-3;

if val1 < 0
    if verb
    fprintf('[Debug] min(Ds - abs(Xc - Cs)) = %g (at index %d)\n', val1, index);
    end
    safetycheck = 0;
else
    for ind_unsafe = 1:Nu
        [val2, idx2] = min(DuArray(:,ind_unsafe) - abs(Xc - CuArray(:,ind_unsafe)));

        if val2 >= tol
            if verb
                fprintf('[Debug] min(DuArray - abs(Xc - CuArray(:,%d))) = %g (at index %d)\n', ind_unsafe, val2, index);
                disp(['    DuArray(:,', num2str(ind_unsafe), ') = ']);
                disp(DuArray(:, ind_unsafe));

                disp(['    CuArray(:,', num2str(ind_unsafe), ') = ']);
                disp(CuArray(:, ind_unsafe));
            end
            % fprintf('[Debug] min(DuArray - abs(Xc - CuArray(:,%d))) = %g (at index %d)\n', ind_unsafe, val2, index);
            safetycheck = 0;
            break;
        end
    end
end

end
