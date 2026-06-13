function prob = nmpc_opti(NMPCConfig, QuadConfig)
    import casadi.*

    % unpack
    NX = NMPCConfig.NX; 
    NU = NMPCConfig.NU; 
    H  = NMPCConfig.H;  
    dt = NMPCConfig.dt;
    m  = QuadConfig.m;  
    g  = QuadConfig.g;  
    e3 = QuadConfig.e3;

    % decision variables
    opti = Opti();
    X = opti.variable(NX, H+1);
    U = opti.variable(NU, H);

    % parameters 
    Xref = opti.parameter(NX, H+1);
    x0   = opti.parameter(NX, 1);

    % dynamics 
    [x_sym, u_sym, f_sym] = eom_so3(QuadConfig);   
    f_func = casadi.Function('f_dyn', {x_sym, u_sym}, {f_sym});

    % initial condition
    opti.subject_to(X(:,1) == x0);

    % cost
    cost = 0;

    % constraints + RK4 dynamics
    for k = 1:H
        xk = X(:,k);  
        uk = U(:,k);  
        xk1 = X(:,k+1);

        % RK4
        k1 = f_func(xk,                 uk);
        k2 = f_func(xk + 0.5*dt*k1,     uk);
        k3 = f_func(xk + 0.5*dt*k2,     uk);
        k4 = f_func(xk + dt*k3,         uk);
        x_next_rk4 = xk + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
        opti.subject_to(xk1 == x_next_rk4);

        % input constraints
        opti.subject_to(QuadConfig.thrust_min <= uk(1) <= QuadConfig.thrust_max);
        % opti.subject_to(QuadConfig.torque_min(:) <= uk(2:4) <= QuadConfig.torque_max(:));

        % velocity constraints
        opti.subject_to(-QuadConfig.vel_max(:) <= xk(4:6) <= QuadConfig.vel_max(:));

        % acceleration constraints (nonlinear, optional)
        xi_k = xk(7:9);
        Rk   = so3_exp(xi_k);      
        a_k  = (uk(1)/m) * (Rk*e3) - g*e3;
        opti.subject_to(-QuadConfig.acc_max(:) <= a_k <= QuadConfig.acc_max(:));

        % 
        % % % stage cost
        xi_ref_k = Xref(7:9, k);
        % Rref_k   = so3_exp(xi_ref_k);
        % xi_e     = so3_log(Rref_k' * Rk);
        xi_e = xi_k - xi_ref_k;
        x_err    = [ X(1:6, k) - Xref(1:6, k);
                     xi_e;
                     X(10:12, k) - Xref(10:12, k) ];
        cost = cost + x_err' * NMPCConfig.Q * x_err + uk' * NMPCConfig.R * uk;

    end

    % terminal cost
    x_T      = X(7:9, H+1);
    x_ref_T  = Xref(7:9, H+1);
    % RT       = so3_exp(x_T);
    % Rref_T   = so3_exp(x_ref_T);
    % xi_e     = so3_log(Rref_T' * RT);
    xi_e = x_T - x_ref_T;
    x_err    = [X(1:6, H+1) - Xref(1:6, H+1);
                 xi_e;
                 X(10:12, H+1) - Xref(10:12, H+1) ];
    cost = cost + x_err'* NMPCConfig.QT *x_err;

    opti.minimize(cost);

    % solver options (Opti plugin opts in 2nd arg; Ipopt opts flat in 3rd arg)
    p_opts = struct( ...
        'print_time', false, ...
        'expand', true ...            % optional: can improve speed/stability
    );
    
    s_opts = struct( ...
        'max_iter', NMPCConfig.max_iter, ...
        'tol',      NMPCConfig.tol, ...
        'warm_start_init_point', 'yes', ...
        'print_level', 0, ...         % quiet Ipopt
        'sb', 'yes' ...               % short banner
        ... % 'linear_solver','mumps' % (optional) if needed
    );
    
    opti.solver('ipopt', p_opts, s_opts);
    
    % pack
    prob.opti  = opti;
    prob.X     = X;
    prob.U     = U;
    prob.Xref  = Xref;
    prob.x0    = x0;
    prob.f_dyn = f_func;   % for rollout
    prob.dt    = dt;
end
