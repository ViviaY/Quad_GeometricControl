function [initial_unsatisfy, initial_outsideBox] = get_initial_unsatisfypoints( ...
            anneal_options, param, annealing_output, delta, init_n, ...
            generatedTraj, show_count, save_dir)

rng("default");

Points_Array = generatedTraj.Points_Array;
tau = generatedTraj.tau;
t = 0;

% Controller gains
k.x = annealing_output.opt_k(1);  
k.v = annealing_output.opt_k(2);  
k.R = annealing_output.opt_k(3);  
k.W = annealing_output.opt_k(4);


% 1. unsatisfy bounds 
initial_unsatisfy.init_n = init_n;
initial_unsatisfy.initial_points = zeros(18, init_n);
initial_unsatisfy.initial_eul   = zeros(3, init_n);
initial_unsatisfy.psi_0_list    = zeros(1, init_n);
initial_unsatisfy.V1_0_list     = zeros(1, init_n);
initial_unsatisfy.eW0_norm_list = zeros(1, init_n);
unsatisfy_count = 0;

% 2. unsatisfy outsideBox
initial_outsideBox.init_n = init_n;
initial_outsideBox.initial_points = zeros(18, init_n);
initial_outsideBox.initial_eul   = zeros(3, init_n);
initial_outsideBox.psi_0_list    = zeros(1, init_n);
initial_outsideBox.V1_0_list     = zeros(1, init_n);
initial_outsideBox.eW0_norm_list = zeros(1, init_n);
outsideBox_count = 0;

xd_0 = DesiredTrajectory(t, Points_Array, tau, param);

% set up the scene box that is larger than the initial safe set. 
lb = -0.8;  
ub = 2;
scene_center = (lb + ub)/2 * ones(3,1);
scene_half   = (ub - lb)/2 * ones(3,1);
safe_half = generatedTraj.R_Array(:,1);

MAX_ITERS = 1e3;
iters = 0;

while (unsatisfy_count < init_n) || (outsideBox_count < init_n)

    iters = iters + 1;
    if iters > MAX_ITERS
        warning("Max iterations reached.");
        break;
    end

    % Sample state
    x0 = scene_center + scene_half .* (2*rand(3,1) - 1); % x0 within a scene box
    v0 = 0.4 + (0.7 - 0.4) * rand(3,1); % v0 in [0.4, 0.7]

    eul0 = delta.R*(2*rand(3,1)-1);
    R0 = expm(hat(eul0));

    W0 = delta.W*(2*rand(3,1)-1);
    X0 = [x0; v0; W0; reshape(R0,9,1)];

    % Compute errors 
    [~, ~, error0, calculated0] = position_control(X0, xd_0, annealing_output, param);
    ex_0 = error0.x; 
    ev_0 = error0.v;
    eR_0 = error0.R;
    eW_0 = error0.W;

    Rd_0 = calculated0.Rd;
    psi_0 = get_psi(R0, Rd_0);

    [V1_0, ~] = lyapunov(param, k, ex_0, ev_0, eR_0, eW_0, psi_0);

    % Conditions for satisfy 
    condition1 = (psi_0 <= anneal_options.alpha_psi * anneal_options.psi_bar);
    condition2 = (0.5*eW_0'*param.J*eW_0 <= k.R * (1 - anneal_options.alpha_psi) * anneal_options.psi_bar);
    condition3 = (V1_0 <= anneal_options.V1_0);

    is_satisfy = (condition1 && condition2 && condition3);
    is_outsidebounds = norm(ex_0) >= 0.35 && norm(ev_0) >= 0.7;
    is_outsideBox = (any(abs(x0 - xd_0.x) > safe_half)) && (norm(ev_0) >= 0.7);

    % 1: unsatisfy and outside the bounds
    if ~is_satisfy && is_outsidebounds && unsatisfy_count < init_n

        unsatisfy_count = unsatisfy_count + 1;

        initial_unsatisfy.initial_points(:, unsatisfy_count) = X0;
        initial_unsatisfy.initial_eul(:, unsatisfy_count) = eul0;
        initial_unsatisfy.psi_0_list(unsatisfy_count) = psi_0;
        initial_unsatisfy.V1_0_list(unsatisfy_count) = V1_0;
        initial_unsatisfy.eW0_norm_list(unsatisfy_count) = 0.5*eW_0'*param.J*eW_0;

        if show_count
            fprintf("Saved unsatisfy: %d/%d\n", unsatisfy_count, init_n);
        end
    end


    % 2: outsideBox
    if ~is_satisfy && is_outsideBox && outsideBox_count < init_n

        outsideBox_count = outsideBox_count + 1;

        initial_outsideBox.initial_points(:, outsideBox_count) = X0;
        initial_outsideBox.initial_eul(:, outsideBox_count) = eul0;
        initial_outsideBox.psi_0_list(outsideBox_count) = psi_0;
        initial_outsideBox.V1_0_list(outsideBox_count) = V1_0;
        initial_outsideBox.eW0_norm_list(outsideBox_count) = 0.5*eW_0'*param.J*eW_0;

        if show_count
            fprintf("Saved outsideBox: %d/%d\n", outsideBox_count, init_n);
        end
    end

end

% Save files
save(fullfile(save_dir, "InitialPoints_unsatisfy_outsideBouds.mat"), "initial_unsatisfy");
save(fullfile(save_dir, "InitialPoints_unsatisfy_outsideBox.mat"), "initial_outsideBox");

end
