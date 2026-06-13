function initial = get_initial_points(anneal_options, param, annealing_output, delta, ...
                                      init_n, generatedTraj, show_count, save_dir)
rng("default");

Points_Array = generatedTraj.Points_Array;
tau = generatedTraj.tau;
t = 0;

% Controller gains
k.x = annealing_output.opt_k(1);  
k.v = annealing_output.opt_k(2);  
k.R = annealing_output.opt_k(3);  
k.W = annealing_output.opt_k(4); 
c1 = annealing_output.opt_bounds.c1; 
c2 = annealing_output.opt_bounds.c2; 

initial.init_n = init_n;
initial.initial_points = zeros(18, init_n);
initial.initial_eul = zeros(3, init_n);
initial.psi_0_list = zeros(1,init_n);
initial.V1_0_list = zeros(1,init_n);
initial.eW0_norm_list = zeros(1,init_n);


initial_point_count = 1;
runtimes = 0;


while initial_point_count ~= init_n+1

    runtimes = runtimes + 1;
    
    % obtainning the initial point of the desired trajectory
    xd_0 = DesiredTrajectory(t, Points_Array, tau, param);

    % create the other initial points with delta
    x0 = xd_0.x + delta.x*(2*rand(3,1)-1);  % [0, 0, 0]';
    v0 = [0, 0, 0]' + delta.v*(2*rand(3,1)-1);  %  x1.v;  %  [0, 0, 0]';
    eul0 = delta.R*(2*rand(3,1)-1);
    R0 = expm(hat(eul0));
    W0 = [0, 0, 0]' + delta.W*(2*rand(3,1)-1);
    X0 = [x0; v0; W0; reshape(R0,9,1)];
    
    % calculate the initial errors
    [~, ~, error0, calculated0] = position_control(X0, xd_0, annealing_output, param);
    ex_0 = error0.x;
    ev_0 = error0.v;
    eR_0 = error0.R;
    eW_0 = error0.W;
    Rd_0 = calculated0.Rd;
    psi_0 = get_psi(R0, Rd_0);
    [V1_0, ~] = lyapunov(param, k, ex_0, ev_0, eR_0, eW_0, psi_0);

    
    % testing if the initial errors satisfy the initial conditions
    condition1 = (psi_0 <= anneal_options.alpha_psi * anneal_options.psi_bar);
    condition2 = (0.5*eW_0'*param.J*eW_0 <= k.R * (1 - anneal_options.alpha_psi) * anneal_options.psi_bar);
    condition3 = (V1_0 <= anneal_options.V1_0);
    

    if condition1 && condition2 && condition3
        
        initial.psi_0_list(initial_point_count) = psi_0;
        initial.V1_0_list(initial_point_count) = V1_0;
        initial.eW0_norm_list(initial_point_count) = 0.5*eW_0'*param.J*eW_0;
        initial.initial_points(:, initial_point_count) = X0;
        initial.initial_eul(:, initial_point_count) = eul0;
        initial_point_count = initial_point_count + 1;

        if show_count
            fprintf('Succeeded %d/%d. \n', initial_point_count - 1, init_n);
        end
    else
        % if not, show which conditions is not satisfied
        if show_count
            if ~condition1
                fprintf('Runtimes = %d, condition 1 fails. \n', runtimes);
            end
            if ~condition2
                fprintf('Runtimes = %d, condition 2 fails. \n', runtimes);
            end
            if ~condition3
                fprintf('Runtimes = %d, condition 3 fails. \n', runtimes);
            end
        end
    end

end

save(fullfile(save_dir, "InitialPoints.mat"), "initial");

end